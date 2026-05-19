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
# émerger un trade-off -> pas si urgent
#
# adapté aux autres algos d'assignation
# + ajouter un argument nperm pour faire un certain nombre de permutations en
# complément du k-fold -> déjà un peu trop long

#  lca_res <- assign_vsearch_lca(fake_pq, ref_fasta= paste0(tempdir(), "/", "test_refseq.fasta"))
# blast_res_tophit <- assign_blastn(fake_pq, ref_fasta= paste0(tempdir(), "/", "test_refseq.fasta"),  keep_blast_metrics = TRUE, method="top-hit")
#  blast_res_vote <- assign_blastn(fake_pq, ref_fasta= paste0(tempdir(), "/", "test_refseq.fasta"), method = "vote", vote_algorithm= "consensus",  keep_blast_metrics = TRUE)
#  dada2_res <- assign_dada(fake_pq, ref_fasta= paste0(tempdir(), "/", "test_refseq.fasta"))

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
                      ...) {
  dna <- Biostrings::readDNAStringSet(ref_fasta)
  method <- match.arg(method)

  if (!method %in% c("sintax", "dada2") && length(min_bootstrap) > 1) {
    stop("min_bootstrap must be set to one value (not a vector) exept if method
         is set to 'sintax' or 'dada2'")
  }
  if (!is.null(seed)) {
    set.seed(seed)
    newseed <- round(runif(1, 1, 1e+09))
    on.exit(set.seed(newseed))
  }

  dna_shuffled <- dna[sample(length(dna)), ]
  folds <- cut(seq(1, length(dna_shuffled)), breaks = fold_number, labels = FALSE)

  res <- list()
  for (f in 1:fold_tested) {
    if (verbose) {
      print(paste0(f, "/", fold_tested))
    }
    index_tested <- which(folds == f, arr.ind = TRUE)
    tested_data <- dna_shuffled[index_tested, ]

    fake_pq <- create_fake_pq_from_refseq(tested_data)
    if (!is.null(patterns_NA)) {
      fake_pq <- taxtab_replace_pattern_by_NA(fake_pq, patterns = patterns_NA, ignore.case = ignore.case)
    }

    if (remove_tested_sequences) {
      Biostrings::writeXStringSet(dna_shuffled[-index_tested],
                                  paste0(tempdir(), "/", "test_refseq.fasta"))
    } else {
      Biostrings::writeXStringSet(dna_shuffled,
                                  paste0(tempdir(), "/", "test_refseq.fasta"))
    }

    if (method == "sintax") {
      assign_res <- assign_sintax(
        fake_pq,
        ref_fasta = paste0(tempdir(), "/", "test_refseq.fasta"),
        nproc = nproc,
        behavior = "return_matrix",
        ...
      )

    } else if (method == "lca") {
      assign_res <- assign_vsearch_lca(
        fake_pq,
        ref_fasta = paste0(tempdir(), "/", "test_refseq.fasta"),
        nproc = nproc,
        behavior = "return_matrix",
        ...
      )
    } else if (method == "blastn") {
      assign_res <- list()
      assign_res$taxo_value <- assign_blastn(fake_pq,
                                             ref_fasta = paste0(tempdir(), "/", "test_refseq.fasta"),
                                             behavior = "add_to_phyloseq",
                                             ...)@tax_table

      assign_res$taxo_value <- assign_res$taxo_value |>
        data.frame() |>
        tibble() |>
        select(ends_with("_blastn")) |>
        select(-Taxa_name_db_blastn)

      colnames(assign_res$taxo_value) <- colnames(fake_pq@tax_table)


    } else if (method == "dada2_2steps") {
      # Renamed from "dada2_steps" so the branch is reachable from match.arg.
      # Body still stops — wire up assign_dada2 in MiscMetabar before enabling.
      stop("method dada2_2steps is not working for the moment")
      assign_res_pq <- assign_dada2(
        fake_pq,
        ref_fasta = paste0(tempdir(), "/", "test_refseq.fasta"),
        nproc = nproc,
        ...
      )
      assign_res <- list()
      assign_res$taxo_value <- assign_res_pq@tax_table
    } else if (method == "dada2") {
      dbpq::format2dada2(
        fasta_db = paste0(tempdir(), "/", "test_refseq.fasta"),
        output_path = paste0(tempdir(), "/", "test_refseq_dada.fasta"),
        pattern_to_remove = "\\|rep.*"
      )
      assign_res_dada <- assignTaxonomy(
        fake_pq@refseq,
        refFasta = paste0(tempdir(), "/", "test_refseq_dada.fasta"),
        outputBootstraps = TRUE,
        minBoot = 0,
        ...
      )
      assign_res <- list()
      assign_res$taxo_value <- as_tibble(assign_res_dada$tax, .name_repair	= "minimal")
      colnames(assign_res$taxo_value)  <-  c(colnames(assign_res_dada$tax)[!is.na(colnames(assign_res_dada$tax))], "taxa_names")
      assign_res$taxo_bootstrap <- as_tibble(assign_res_dada$boot, .name_repair	= "minimal")
      colnames(assign_res$taxo_bootstrap)  <-  c(colnames(assign_res_dada$boot)[!is.na(colnames(assign_res_dada$boot))], "taxa_names")
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
        tax_tib <- as_tibble(as.matrix(unclass(fake_pq@tax_table)))

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
        res_assign_NA[i, ] <- colSums(NA_matrix) / length(index_tested)

        res_assign_NA_classif[i, ]  <- colSums(is.na(assign_res$taxo)) / length(index_tested)
        res_assign_NA_database[i, ] <- colSums(is.na(tax_tib)) / length(index_tested)

        good_classification_matrix <- assign_res$taxo == tax_tib
        res_assign_good_classification[i, ] <-
          colSums(good_classification_matrix, na.rm =   TRUE) / length(index_tested)

        bad_classification_matrix <- assign_res$taxo != tax_tib
        res_assign_bad_classification[i, ] <-
          colSums(bad_classification_matrix, na.rm =  TRUE) / length(index_tested)

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
        mutate(bootstrap = min_bootstrap[[i]]) |>
        mutate(fold = as.character(f))

      colnames(res_assign_good_classification) <- colnames(assign_res$taxo)
      res_assign_good_classification <- res_assign_good_classification |>
        as_tibble() |>
        mutate(bootstrap = min_bootstrap[[i]]) |>
        mutate(fold = as.character(f))

      colnames(res_assign_bad_classification) <- colnames(assign_res$taxo)
      res_assign_bad_classification <- res_assign_bad_classification |>
        as_tibble() |>
        mutate(bootstrap = min_bootstrap[[i]]) |>
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
      tax_tib <- as_tibble(as.matrix(unclass(fake_pq@tax_table)))

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
      res_assign_NA[1, ] <- colSums(NA_matrix) / length(index_tested)

      res_assign_NA_classif[1, ]  <- colSums(is.na(assign_res$taxo)) / length(index_tested)
      res_assign_NA_database[1, ] <- colSums(is.na(tax_tib)) / length(index_tested)

      good_classification_matrix <- assign_res$taxo == tax_tib
      res_assign_good_classification[1, ] <-
        colSums(good_classification_matrix, na.rm =   TRUE) / length(index_tested)

      bad_classification_matrix <- assign_res$taxo != tax_tib
      res_assign_bad_classification[1, ] <-
        colSums(bad_classification_matrix, na.rm =  TRUE) / length(index_tested)

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
