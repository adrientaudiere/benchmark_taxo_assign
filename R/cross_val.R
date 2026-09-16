# k-fold cross-validation of taxonomic-assignment algorithms on a reference
# fasta. Not part of the main pipeline — source manually when you want to
# benchmark a method against a database (no mock community required).
#
# Requires dplyr / tibble / tidyr / Biostrings / MiscMetabar / dada2 on the
# search path. In the benchmark project these are loaded by the pipeline
# scripts; if running interactively, load them first.

# Wrapper: run cross_val() over a vector of min_bootstrap values for methods
# that do not natively accept a vector min_bootstrap (i.e. anything other than
# "sintax" / "dada2"). For sintax/dada2 prefer passing a vector directly to
# cross_val() — it is faster because the heavy assignment runs once per fold.
cross_val_param <- function(..., min_bootstrap = c(0.4, 0.5, 0.6)) {
  res <- list()
  for (p in min_bootstrap) {
    res[[as.character(p)]] <- cross_val(..., min_bootstrap = p)
  }
  res
}

# TODO
# réfléchir à l'inclusion de fake pour faire un TRUE negative et donc voire
# émerger un trade-off 
#
# adapté aux autres algos d'assignation
# TODO? + ajouter un argument nperm pour faire un certain nombre de permutations en
# complément du k-fold -> déjà un peu trop long

#  lca_res <- assign_vsearch_lca(fake_pq, ref_fasta= tmp_fasta)
# blast_res_tophit <- assign_blastn(fake_pq, ref_fasta= tmp_fasta,  keep_blast_metrics = TRUE, method="top-hit")
#  blast_res_vote <- assign_blastn(fake_pq, ref_fasta= tmp_fasta, method = "vote", vote_algorithm= "consensus",  keep_blast_metrics = TRUE)
#  dada2_res <- assign_dada(fake_pq, ref_fasta= tmp_fasta)

#' Cross validation of taxonomic assignation algorithm on a given fasta database
#'
#' @param ref_fasta
#' @param fold_number (int) Number of fold define the number of cut in the database
#' @param fold_tested (int) Set the number of fold for which we test the assignation
#'   method. Default = fold_number
#' @param patterns_NA
#' @param ignore.case
#' @param seed
#' @param remove_tested_sequences
#' @param min_bootstrap (Float or vector of Float, default 0.5). A value of minimum
#'   bootstrap for method sintax or dada2. If min_bootstrap is a vector, all
#'   value of minimum_bootstrap are used and return in the output data frames.
#' @param compute_by_tax_level
#' @param verbose
#' @param max_seq (int) Size of the random subsample of the database (NULL:
#'   every record). With trimmed queries, the number of queries.
#' @param reduce_reference (logical, default TRUE) TRUE reduces the reference
#'   to the drawn pool, so a run searches about `oversample * max_seq` records
#'   whatever the database (the behaviour of the runs made before 2026-09-16,
#'   ROADMAP B24). FALSE keeps the whole database and removes only the tested
#'   queries of each fold, as Bokulich et al. 2018; the queries are drawn the
#'   same way in both cases.
#' @param primer_fw,primer_rev (character) Primers of the amplicon. When both
#'   are set, the queries are trimmed with trim_cv_queries() (R/cv_queries.R)
#'   and the records without the reverse-primer site are never queried: they
#'   stay in the training part, which keeps full-length records.
#' @param primer_min_overlap (int) cutadapt minimum overlap (-O) for the
#'   trimming.
#' @param oversample (numeric, default 2) Records drawn per wanted query when
#'   the queries are trimmed.
#' @param cutadapt_prelude (character) Shell prelude activating cutadapt (NULL:
#'   the dbpq default).
#' @param query_fasta (character, default NULL) Fasta from which the queries
#'   are drawn and trimmed when the records of ref_fasta no longer hold the
#'   primer sites (the _Fungi source of a _Fungi_cut database). Queries are
#'   paired by name with the ref_fasta records, and a query whose amplicon
#'   differs from its paired record is dropped. NULL or ref_fasta: the queries
#'   come from ref_fasta.
#' @param id_fasta,query_id_fasta (character, default NULL) sintax-format
#'   fasta holding the same records as ref_fasta (query_fasta) in the same
#'   order. When it differs from ref_fasta (a dada2-format file, whose headers
#'   are the taxonomy only), the records are named by its identifiers, so that
#'   records sharing a taxonomy are not deduplicated; the dada2 headers are
#'   kept for the truth table and the training fasta (ROADMAP B22).
#' #TODO complete documentation 
#'
#' @returns
#' @export
#'
#' @examples
cross_val <- function(ref_fasta,
                      fold_number = 10,
                      fold_tested = fold_number,
                      method = c("sintax", "lca", "blastn", "dada2_2steps", "dada2"),
                      patterns_NA = NULL,
                      min_bootstrap = 0.5,
                      ignore.case = TRUE,
                      seed = NULL,
                      remove_tested_sequences = TRUE,
                      compute_by_tax_level = FALSE,
                      verbose = FALSE,
                      nproc = 1,
                      max_seq = NULL,
                      min_seq_length = 50L,
                      reduce_reference = TRUE,
                      primer_fw = NULL,
                      primer_rev = NULL,
                      primer_min_overlap = NULL,
                      oversample = 2,
                      cutadapt_prelude = NULL,
                      query_fasta = NULL,
                      id_fasta = NULL,
                      query_id_fasta = NULL,
                      ...) {
  dna <- Biostrings::readDNAStringSet(ref_fasta)
  method <- match.arg(method)

  # dada2-format headers are the taxonomy only: without identifiers, the name
  # deduplication below kept one record per taxonomy string for dada2 until
  # 2026-09-15 (ROADMAP B22). The records are named by the identifiers of
  # id_fasta and the dada2 headers are kept aside in `headers`.
  headers <- NULL
  if (!is.null(id_fasta) &&
      normalizePath(id_fasta) != normalizePath(ref_fasta)) {
    ids <- record_ids(id_fasta, widths = Biostrings::width(dna))
    headers <- stats::setNames(names(dna), ids)
    names(dna) <- ids
  }

  if (!is.null(min_seq_length)) {
    dna <- dna[Biostrings::width(dna) >= min_seq_length]
  }

  dna <- dna[!duplicated(names(dna))]

  if (!method %in% c("sintax", "dada2") && length(min_bootstrap) > 1) {
    stop("min_bootstrap must be set to one value (not a vector) exept if method
         is set to 'sintax' or 'dada2'")
  }
  # The seed is set before the subsample is drawn (before 2026-09-15 the
  # subsample followed the targets seed, not `seed`; ROADMAP S8.6).
  if (!is.null(seed)) {
    set.seed(seed)
    newseed <- round(runif(1, 1, 1e+09))
    on.exit(set.seed(newseed))
  }

  # Trimmed queries (ROADMAP D1a, 2026-09-15): the drawn records are cut with
  # cutadapt and the subsample stops at the max_seq-th record holding the
  # reverse-primer site. The training part always gets the ref_fasta records.
  # A _Fungi_cut database no longer holds the primer sites: its queries are
  # drawn and trimmed from query_fasta (its _Fungi source) and paired by name
  # with the ref_fasta records; the trimmed query must equal its paired record.
  trim_queries <- !is.null(primer_fw) && !is.null(primer_rev)
  if (!is.null(query_fasta) &&
      normalizePath(query_fasta) == normalizePath(ref_fasta)) {
    query_fasta <- NULL
  }
  if (!is.null(query_fasta) && !trim_queries) {
    stop("query_fasta is only used when primer_fw and primer_rev are set.")
  }
  if (is.null(query_fasta)) {
    pool_names <- names(dna)
  } else {
    query_index <- Biostrings::fasta.index(query_fasta)
    query_index$id <- query_index$desc
    if (!is.null(headers)) {
      query_index$id <- record_ids(query_id_fasta, widths = query_index$seqlength)
    }
    query_index <- query_index[
      !duplicated(query_index$id) & query_index$id %in% names(dna),
    ]
    pool_names <- query_index$id
  }
  if (!is.null(max_seq) && length(pool_names) > max_seq) {
    n_draw <- if (trim_queries) {
      min(length(pool_names), ceiling(oversample * max_seq))
    } else {
      max_seq
    }
    pool_names <- pool_names[sample(length(pool_names), n_draw)]
  }
  if (is.null(query_fasta)) {
    source_dna <- dna[match(pool_names, names(dna))]
  } else {
    # Rows read in file order, then named by identifier and put back in the
    # drawn order (dada2 headers are not unique, so names cannot be matched).
    rows <- sort(match(pool_names, query_index$id))
    source_dna <- Biostrings::readDNAStringSet(query_index[rows, ])
    names(source_dna) <- query_index$id[rows]
    source_dna <- source_dna[match(pool_names, names(source_dna))]
  }
  if (trim_queries) {
    queries <- trim_cv_queries(
      source_dna,
      primer_fw = primer_fw,
      primer_rev = primer_rev,
      min_overlap = primer_min_overlap,
      nproc = nproc,
      prelude = cutadapt_prelude
    )
    if (!is.null(min_seq_length)) {
      queries <- queries[Biostrings::width(queries) >= min_seq_length]
    }
    if (!is.null(query_fasta)) {
      paired <- dna[match(names(queries), names(dna))]
      queries <- queries[as.character(queries) == as.character(paired)]
    }
    kept <- cv_select_queries(pool_names, names(queries), max_seq)
    dna <- cv_reference_records(dna, kept$pool, reduce_reference)
    queries <- queries[match(kept$queries, names(queries))]
  } else {
    dna <- cv_reference_records(dna, names(source_dna), reduce_reference)
    queries <- source_dna
  }

  if (length(queries) < fold_number) {
    stop(
      "Only ", length(queries), " query sequences remain after length filtering ",
      "(min_seq_length=", min_seq_length, ") and primer trimming; need at least ",
      "fold_number=", fold_number, ". The reference database may be unsuitable ",
      "for this method."
    )
  }

  queries_shuffled <- queries[sample(length(queries))]
  folds <- cut(seq_along(queries_shuffled), breaks = fold_number, labels = FALSE)

  res <- list()
  for (f in 1:fold_tested) {
    if (verbose) {
      print(paste0(f, "/", fold_tested))
    }
    tested_data <- queries_shuffled[folds == f]
    tested_names <- names(tested_data)
    tested_data <- tested_data[!duplicated(as.character(tested_data))]
    tested_data <- tested_data[!duplicated(names(tested_data))]
    n_tested <- length(tested_data)
    # One reference fasta per fold; tempdir() is cleaned when the session ends.
    tmp_fasta <- tempfile(pattern = "cv_refseq_", fileext = ".fasta")

    fake_pq <- create_fake_pq_from_refseq(
      tested_data,
      headers = if (is.null(headers)) NULL else unname(headers[names(tested_data)])
    )
    if (!is.null(patterns_NA)) {
      fake_pq <- taxtab_replace_pattern_by_NA(fake_pq, patterns = patterns_NA, ignore.case = ignore.case)
    }

    # Training part: the records of ref_fasta (with or without primer site),
    # under their original headers (dada2 reads the taxonomy from them).
    training <- if (remove_tested_sequences) {
      dna[!names(dna) %in% tested_names]
    } else {
      dna
    }
    if (!is.null(headers)) {
      names(training) <- unname(headers[names(training)])
    }
    Biostrings::writeXStringSet(training, tmp_fasta)

    if (method == "sintax") {
      # assign_sintax() applies its own min_bootstrap (default 0.5) to
      # taxo_value: pass ours, or 0 when several values are filtered below.
      assign_res <- assign_sintax(
        fake_pq,
        ref_fasta = tmp_fasta,
        nproc = nproc,
        behavior = "return_matrix",
        min_bootstrap = if (length(min_bootstrap) > 1) 0 else min_bootstrap,
        ...
      )
      # vsearch --sintax with several threads does not keep the query order;
      # the metrics below compare rows by position (ROADMAP B23).
      assign_res$taxo_value <- cv_align_rows(assign_res$taxo_value, taxa_names(fake_pq))
      assign_res$taxo_bootstrap <- cv_align_rows(assign_res$taxo_bootstrap, taxa_names(fake_pq))

    } else if (method == "lca") {
      assign_res <- list()
      lca_raw <- assign_vsearch_lca(
        fake_pq,
        ref_fasta = tmp_fasta,
        nproc = nproc,
        behavior = "return_matrix",
        ...
      )
      if (is.null(lca_raw)) {
        stop("assign_vsearch_lca returned NULL — no LCA output produced.")
      }
      # Query order, unmatched sequences as NA rows (not missing rows).
      assign_res$taxo_value <- cv_align_rows(lca_raw, taxa_names(fake_pq)) |>
        # assign_vsearch_lca() names its ranks "<rank>_sintax"; use the plain
        # rank names so lca rows line up with the other methods (ROADMAP S1.4).
        dplyr::rename_with(\(x) sub("_sintax$", "", x), -taxa_names)
    } else if (method == "blastn") {
      assign_res <- list()
      assign_res$taxo_value <- assign_blastn(fake_pq,
                                             ref_fasta = tmp_fasta,
                                             behavior = "add_to_phyloseq",
                                             nproc = nproc,
                                             ...) |>
        tidypq::tax_table_to_df(convert = FALSE) |>
        select(-taxon) |>
        select(ends_with("_blastn")) |>
        select(-any_of("Taxa_name_db_blastn"))

      if (ncol(assign_res$taxo_value) == 0) {
        # assign_blastn returned physeq unchanged (no BLAST hits at filter thresholds)
        assign_res$taxo_value <- tibble::as_tibble(matrix(
          NA_character_,
          nrow = n_tested,
          ncol = length(phyloseq::rank_names(fake_pq)),
          dimnames = list(NULL, phyloseq::rank_names(fake_pq))
        ))
      } else {
        colnames(assign_res$taxo_value) <- phyloseq::rank_names(fake_pq)
      }


    } else if (method == "dada2_2steps") {
      # Renamed from "dada2_steps" so the branch is reachable from match.arg.
      # Body still stops — wire up assign_dada2 in MiscMetabar before enabling.
      stop("method dada2_2steps is not working for the moment")
      assign_res_pq <- assign_dada2(
        fake_pq,
        ref_fasta = tmp_fasta,
        nproc = nproc,
        ...
      )
      assign_res <- list()
      assign_res$taxo_value <- assign_res_pq@tax_table
    } else if (method == "dada2") {
      # ref_fasta is already in dada2 format (db_path = dada2_format/…), so the
      # subset written to tmp_fasta requires no conversion. dada2:: because
      # crew workers do not attach dada2 (an Import of MiscMetabar).
      # minBoot (0-100) applies min_bootstrap; with several values, every
      # call is kept and filtered below on the rescaled bootstraps.
      assign_res_dada <- dada2::assignTaxonomy(
        fake_pq@refseq,
        refFasta = tmp_fasta,
        outputBootstraps = TRUE,
        minBoot = if (length(min_bootstrap) > 1) 0 else 100 * min_bootstrap,
        multithread = nproc,
        ...
      )
      assign_res <- list()
      assign_res$taxo_value <- assign_res_dada$tax |>
        as.data.frame() |>
        tibble::rownames_to_column("taxa_names") |>
        tibble::as_tibble() |>
        # UNITE dada2 references return "k__Fungi"-style values; the truth
        # table has prefixes removed by simplify_taxo() (ROADMAP S1.4).
        dplyr::mutate(dplyr::across(-taxa_names, \(x) sub("^[a-z]__", "", x)))
      assign_res$taxo_bootstrap <- assign_res_dada$boot |>
        as.data.frame() |>
        tibble::rownames_to_column("taxa_names") |>
        tibble::as_tibble() |>
        # assignTaxonomy() bootstraps are 0-100, min_bootstrap is 0-1.
        dplyr::mutate(dplyr::across(-taxa_names, \(x) x / 100))
    }

    if (nrow(assign_res$taxo_value) != n_tested) {
      stop(
        method, " returned ", nrow(assign_res$taxo_value), " rows for ",
        n_tested, " queries."
      )
    }

    if (length(min_bootstrap) > 1) {
      assign_res$taxo <- select(assign_res$taxo_value, -taxa_names)
      assign_res$taxo_bootstrap  <- select(assign_res$taxo_bootstrap, -taxa_names)

      res_assign_NA <- matrix(
        nrow = length(min_bootstrap),
        ncol = ncol(assign_res$taxo_bootstrap)
      )
      res_assign_NA_classif <- matrix(
        nrow = length(min_bootstrap),
        ncol = ncol(assign_res$taxo_bootstrap)
      )
      res_assign_NA_database <- matrix(
        nrow = length(min_bootstrap),
        ncol = ncol(assign_res$taxo_bootstrap)
      )
      res_assign_good_classification <- matrix(
        nrow = length(min_bootstrap),
        ncol = ncol(assign_res$taxo_bootstrap)
      )
      res_assign_bad_classification <- matrix(
        nrow = length(min_bootstrap),
        ncol = ncol(assign_res$taxo_bootstrap)
      )

      if (compute_by_tax_level) {
        tib_by_tax_level <- tibble(.rows = 5)
      }
      for (i in seq_along(min_bootstrap)) {
        tax_tib <- tidypq::tax_table_to_df(fake_pq, convert = FALSE) |> select(-taxon)

        if (!is.null(patterns_NA)) {
          for (pat in patterns_NA) {
            assign_res$taxo <- assign_res$taxo |>
              mutate(across(
                everything(),
                gsub,
                pattern = pat,
                replacement = NA
              ))
          }
        }

        assign_res$taxo[assign_res$taxo_bootstrap < min_bootstrap[i]] <- NA

        NA_matrix <- is.na(assign_res$taxo == tax_tib)
        res_assign_NA[i, ] <- colSums(NA_matrix) / n_tested

        res_assign_NA_classif[i, ]  <- colSums(is.na(assign_res$taxo)) / n_tested
        res_assign_NA_database[i, ] <- colSums(is.na(tax_tib)) / n_tested

        good_classification_matrix <- assign_res$taxo == tax_tib
        res_assign_good_classification[i, ] <-
          colSums(good_classification_matrix, na.rm =   TRUE) / n_tested

        bad_classification_matrix <- assign_res$taxo != tax_tib
        res_assign_bad_classification[i, ] <-
          colSums(bad_classification_matrix, na.rm =  TRUE) / n_tested

        if (compute_by_tax_level) {
          for (taxlev in colnames(assign_res$taxo[, -1]))
          {
            val_rank <- unique(tax_tib[[taxlev]])
            for (tax_rank in val_rank) {
              cond <- tax_tib[, taxlev] == tax_rank
              TP <- sum(good_classification_matrix[cond, taxlev], na.rm = TRUE)
              FP <- sum(!good_classification_matrix[cond, taxlev], na.rm = TRUE)

              cond_fn <- assign_res$taxo[, taxlev] == tax_rank
              FN <- sum(!good_classification_matrix[cond_fn, taxlev], na.rm = TRUE)

              tib_by_tax_level <- rbind(tib_by_tax_level,
                                        c(min_bootstrap[[i]], taxlev, tax_rank, f, "TP", TP))

              tib_by_tax_level <- rbind(tib_by_tax_level,
                                        c(min_bootstrap[[i]], taxlev, tax_rank, f, "FP", FP))
              tib_by_tax_level <- rbind(tib_by_tax_level,
                                        c(min_bootstrap[[i]], taxlev, tax_rank, f, "FN", FN))

              tib_by_tax_level <- rbind(
                tib_by_tax_level,
                c(
                  min_bootstrap[[i]],
                  taxlev,
                  tax_rank,
                  f,
                  "F1_score",
                  2 * TP / (2 * TP + FP + FN)
                )
              )

            }
          }
          colnames(tib_by_tax_level) <- c(
            "bootstrap",
            "taxonomic_rank",
            "taxonomic_value",
            "fold",
            "metric",
            "value"
          )
        }
        if (verbose) {
          print(paste0(round(
            100 * i / length(min_bootstrap), 2
          ), "%"))
        }
      }

      if (sum(
        res_assign_NA + res_assign_bad_classification + res_assign_good_classification
      )
      != (ncol(res_assign_NA) * nrow(res_assign_NA))) {
        stop("The proportion of NA, good classification and bad classification must sum to 1")
      }

      colnames(res_assign_NA) <- colnames(assign_res$taxo)

      res_assign_NA <- res_assign_NA |>
        as_tibble() |>
        mutate(bootstrap = min_bootstrap) |>
        mutate(fold = as.character(f))

      colnames(res_assign_good_classification) <- colnames(assign_res$taxo)
      res_assign_good_classification <- res_assign_good_classification |>
        as_tibble() |>
        mutate(bootstrap = min_bootstrap) |>
        mutate(fold = as.character(f))

      colnames(res_assign_bad_classification) <- colnames(assign_res$taxo)
      res_assign_bad_classification <- res_assign_bad_classification |>
        as_tibble() |>
        mutate(bootstrap = min_bootstrap) |>
        mutate(fold = as.character(f))

      if (compute_by_tax_level) {
        res <- list(
          "good_classifications" =  rbind(
            res$good_classifications,
            res_assign_good_classification
          ),
          "wrong_classifications" = rbind(
            res$wrong_classifications,
            res_assign_bad_classification
          ),
          "prop_NA" = rbind(res$prop_NA, res_assign_NA),
          "metrics_by_taxonomic_level" = rbind(res$metrics_by_taxonomic_level, tib_by_tax_level)
        )

      } else {
        res <- list(
          "good_classifications" = rbind(
            res$good_classifications,
            res_assign_good_classification
          ),
          "wrong_classifications" = rbind(
            res$wrong_classifications,
            res_assign_bad_classification
          ),
          "prop_NA" = rbind(res$prop_NA, res_assign_NA)
        )
      }
    } else {
      assign_res$taxo <- assign_res$taxo_value |>
        select(-any_of(c("taxa_names", "Taxa_name_db_blastn")))

      res_assign_NA <- matrix(nrow = 1,
                              ncol = ncol(assign_res$taxo))
      res_assign_NA_classif <- matrix(nrow = 1,
                                      ncol = ncol(assign_res$taxo))
      res_assign_NA_database <- matrix(nrow = 1,
                                       ncol = ncol(assign_res$taxo))
      res_assign_good_classification <- matrix(nrow = 1,
                                               ncol = ncol(assign_res$taxo))
      res_assign_bad_classification <- matrix(nrow = 1,
                                              ncol = ncol(assign_res$taxo))

      if (compute_by_tax_level) {
        tib_by_tax_level <- tibble(.rows = 5)
      }
      tax_tib <- tidypq::tax_table_to_df(fake_pq, convert = FALSE) |> select(-taxon)

      if (!is.null(patterns_NA)) {
        for (pat in patterns_NA) {
          assign_res$taxo <- assign_res$taxo |>
            mutate(across(
              everything(),
              gsub,
              pattern = pat,
              replacement = NA
            ))
        }
      }

      NA_matrix <- is.na(assign_res$taxo == tax_tib)
      res_assign_NA[1, ] <- colSums(NA_matrix) / n_tested

      res_assign_NA_classif[1, ]  <- colSums(is.na(assign_res$taxo)) / n_tested
      res_assign_NA_database[1, ] <- colSums(is.na(tax_tib)) / n_tested

      good_classification_matrix <- assign_res$taxo == tax_tib
      res_assign_good_classification[1, ] <-
        colSums(good_classification_matrix, na.rm =   TRUE) / n_tested

      bad_classification_matrix <- assign_res$taxo != tax_tib
      res_assign_bad_classification[1, ] <-
        colSums(bad_classification_matrix, na.rm =  TRUE) / n_tested

      if (compute_by_tax_level) {
        for (taxlev in colnames(assign_res$taxo[, -1]))
        {
          val_rank <- unique(tax_tib[[taxlev]])
          for (tax_rank in val_rank) {
            cond <- tax_tib[, taxlev] == tax_rank
            TP <- sum(good_classification_matrix[cond, taxlev], na.rm = TRUE)
            FP <- sum(!good_classification_matrix[cond, taxlev], na.rm = TRUE)

            cond_fn <- assign_res$taxo[, taxlev] == tax_rank
            FN <- sum(!good_classification_matrix[cond_fn, taxlev], na.rm = TRUE)

            tib_by_tax_level <- rbind(tib_by_tax_level, c(taxlev, tax_rank, f, "TP", TP))

            tib_by_tax_level <- rbind(tib_by_tax_level, c(taxlev, tax_rank, f, "FP", FP))
            tib_by_tax_level <- rbind(tib_by_tax_level, c(taxlev, tax_rank, f, "FN", FN))

            tib_by_tax_level <- rbind(tib_by_tax_level,
                                      c(taxlev, tax_rank, f, "F1_score", 2 * TP / (2 * TP + FP + FN)))

          }
        }
        colnames(tib_by_tax_level) <- c("taxonomic_rank",
                                        "taxonomic_value",
                                        "fold",
                                        "metric",
                                        "value")
      }

      if (sum(
        res_assign_NA + res_assign_bad_classification + res_assign_good_classification
      )
      != (ncol(res_assign_NA) * nrow(res_assign_NA))) {
        stop("The proportion of NA, good classification and bad classification must sum to 1")
      }

      colnames(res_assign_NA) <- colnames(assign_res$taxo)
      res_assign_NA <- res_assign_NA |>
        as_tibble() |>
        mutate(fold = as.character(f))

      colnames(res_assign_good_classification) <- colnames(assign_res$taxo)
      res_assign_good_classification <- res_assign_good_classification |>
        as_tibble() |>
        mutate(fold = as.character(f))

      colnames(res_assign_bad_classification) <- colnames(assign_res$taxo)
      res_assign_bad_classification <- res_assign_bad_classification |>
        as_tibble() |>
        mutate(fold = as.character(f))

      if (compute_by_tax_level) {
        res <- list(
          "good_classifications" =  rbind(
            res$good_classifications,
            res_assign_good_classification
          ),
          "wrong_classifications" = rbind(
            res$wrong_classifications,
            res_assign_bad_classification
          ),
          "prop_NA" = rbind(res$prop_NA, res_assign_NA),
          "metrics_by_taxonomic_level" = rbind(res$metrics_by_taxonomic_level, tib_by_tax_level)
        )
      } else {
        res <- list(
          "good_classifications" = rbind(
            res$good_classifications,
            res_assign_good_classification
          ),
          "wrong_classifications" = rbind(
            res$wrong_classifications,
            res_assign_bad_classification
          ),
          "prop_NA" = rbind(res$prop_NA, res_assign_NA)
        )
      }
    }
  }

  final_res <- list()

  if (length(min_bootstrap) > 1) {
    final_res[["good_classifications"]] <- res$good_classifications |>
      pivot_longer(-c(fold, bootstrap)) |>
      group_by(bootstrap, name) |>
      summarise(mean = mean(value), sd = sd(value))

    final_res[["wrong_classifications"]] <- res$wrong_classifications |>
      pivot_longer(-c(fold, bootstrap)) |>
      group_by(bootstrap, name) |>
      summarise(mean = mean(value), sd = sd(value))

    final_res[["prop_NA"]] <- res$prop_NA |>
      pivot_longer(-c(fold, bootstrap)) |>
      group_by(bootstrap, name) |>
      summarise(mean = mean(value), sd = sd(value))

    if (compute_by_tax_level) {
      final_res[["metrics_by_taxonomic_level"]] <- res$metrics_by_taxonomic_level |>
        group_by(bootstrap, taxonomic_rank, taxonomic_value, metric) |>
        summarise(mean = mean(as.numeric(value)),
                  sd = sd(as.numeric(value)))
    }
  } else {
    final_res[["good_classifications"]] <- res$good_classifications |>
      pivot_longer(-c(fold)) |>
      group_by(name) |>
      summarise(mean = mean(value), sd = sd(value))

    final_res[["wrong_classifications"]] <- res$wrong_classifications |>
      pivot_longer(-c(fold)) |>
      group_by(name) |>
      summarise(mean = mean(value), sd = sd(value))

    final_res[["prop_NA"]] <- res$prop_NA |>
      pivot_longer(-c(fold)) |>
      group_by(name) |>
      summarise(mean = mean(value), sd = sd(value))

    final_res[["metrics"]] <-
      rbind(final_res[["good_classifications"]],
            final_res[["wrong_classifications"]],
            final_res[["prop_NA"]])

    final_res[["metrics"]]$metric <- c(
      rep("good_classifications", nrow(final_res[["good_classifications"]])),
      rep("wrong_classifications", nrow(final_res[["wrong_classifications"]])),
      rep("prop_NA", nrow(final_res[["prop_NA"]]))
    )

    if (compute_by_tax_level) {
      final_res[["metrics_by_taxonomic_level"]] <- res$metrics_by_taxonomic_level |>
        group_by(taxonomic_rank, taxonomic_value, metric) |>
        summarise(mean = mean(as.numeric(value)),
                  sd = sd(as.numeric(value)))
    }
  }
  return(final_res)
}
