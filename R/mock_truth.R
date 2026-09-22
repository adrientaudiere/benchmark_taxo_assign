# Per-unit truth of a mock community (ROADMAP 0.6, docs/hleap_2021_metrics.md
# §6). Hleap et al. 2021 give every sequence its own label before scoring; the
# mocks here only came with a set of expected taxa, so each ASV (or OTU) is
# matched against the Sanger sequences of the strains of its mock.
#
# The rule, confirmed by the developer (Q2a): a unit gets the truth of a Sanger
# sequence when their identity is at least `min_identity` AND the shorter of
# the two sequences is covered at `min_cover` — the inclusion goes both ways
# (the Pauvert ASVs sit inside the Sanger reads, the Tedersoo Sanger reads
# inside the ASVs). A unit matching several strains keeps their lowest common
# rank (Q2c) and leaves the matrix below it; a unit matching none stays in the
# matrix with no truth (Q2b).
#
# Names (hleap §6): the submitted name of the strain is the truth, the GBIF
# accepted name counts as correct too at Genus and Species, and the ranks above
# Genus come from the GBIF lineage — for both mocks, the Table S1 lineage of
# Pauvert serving as a cross-check (lineage_disagreements()). The lineages are
# read from a cache written by refresh_gna_lineage_cache(), so building the
# truth needs no network.
#
# Requires dplyr, tibble and vsearch; refresh_gna_lineage_cache() also needs
# the taxinfo checkout and a network connection.
#
# Measured on the production stores (2026-09-22), which is what the developer
# confirmed in hleap §6: 149 of the 196 Pauvert ASVs get a truth (12 ties) and
# 99 of the 238 Tedersoo ASVs (3 ties); on the OTUs, 126 of 143 and 74 of 180.
# Two Pauvert units stop above Genus for a reason worth knowing: one matches a
# strain named only to class (Sordariomycetes sp) and one ties with a
# genus-only name (Fusarium sp) that GBIF places in another class. Ten more
# stop at Genus because their strain is named "Genus sp": species_name()
# gives them no species truth (2026-09-22), so the Pauvert truth_depth reads
# 127 Species, 20 Genus, 1 Class, 1 Phylum and 47 without truth. The
# The GBIF accepted name differs from the submitted one for 28 Pauvert and 7
# Tedersoo units, and is missing for the names GBIF does not resolve. The
# disagreements between the Table S1 lineage and GBIF are listed in
# data/data_final/mock_truth/pauvert_lineage_disagreements.csv (61 rows: 33
# families, 18 orders, 8 classes, 1 genus, 1 species) for the hand check the
# decision asks for.

truth_ranks <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Species"
)

# Sanger sequences of the strains of a mock, from one row of
# `config.R::mock_sanger_tables`. Returns a tibble: `strain` (the identifier
# used in the truth table), `name` (the taxon name as the authors wrote it) and
# `sequence`. Rows without a sequence are dropped, with a message.
read_mock_sanger <- function(cfg, base_dir = here::here()) {
  raw <- utils::read.csv(
    file.path(base_dir, cfg$path),
    sep = cfg$sep,
    check.names = FALSE
  )
  missing_cols <- setdiff(
    c(cfg$unit_col, cfg$name_col, cfg$seq_col),
    names(raw)
  )
  if (length(missing_cols) > 0) {
    stop(
      "read_mock_sanger(): ",
      cfg$path,
      " has no column ",
      paste(missing_cols, collapse = ", ")
    )
  }
  out <- tibble::tibble(
    strain = as.character(raw[[cfg$unit_col]]),
    name = as.character(raw[[cfg$name_col]]),
    sequence = toupper(gsub("[^A-Za-z]", "", raw[[cfg$seq_col]]))
  )
  empty <- !nzchar(out$sequence)
  if (any(empty)) {
    message(
      "- mock truth: ",
      sum(empty),
      " of ",
      nrow(out),
      " strains of ",
      basename(cfg$path),
      " have no Sanger sequence; dropped"
    )
  }
  out[!empty, ]
}

# Binomial of a taxon name, underscored as the `Gen_sp_*` columns of
# analysis/01_load_and_clean.qmd are ("Agaricus_essettei"). Vectorised.
binomial_name <- function(name) {
  out <- rep(NA_character_, length(name))
  known <- !is.na(name)
  parts <- strsplit(
    trimws(gsub("[_[:space:]]+", " ", name[known])),
    " ",
    fixed = TRUE
  )
  out[known] <- vapply(
    parts,
    \(p) paste(utils::head(p[nzchar(p)], 2), collapse = "_"),
    character(1)
  )
  out
}

# Genus of a taxon name (its first word). Vectorised.
genus_name <- function(name) {
  sub("_.*$", "", binomial_name(name))
}

# Epithets that name no species: "Phoma sp" identifies a strain to the genus
# and no further. Such a name must not become a species truth, because
# analysis/01_load_and_clean.qmd turns every "<Genus>_sp" assignment into NA
# (it is not an identification either), so the unit would be unmatchable at the
# species rank and score FN or FP whatever the method answered. The unit stops
# at Genus instead, exactly as a tie between two species of the same genus
# does.
placeholder_epithets <- c("sp", "spp", "cf", "aff", "indet", "incertae")

# Species of a taxon name: its binomial, or NA when the name carries no
# epithet or a placeholder one. Vectorised.
species_name <- function(name) {
  out <- binomial_name(name)
  epithet <- rep(NA_character_, length(out))
  known <- !is.na(out)
  epithet[known] <- vapply(
    strsplit(out[known], "_", fixed = TRUE),
    \(p) if (length(p) >= 2) p[[2]] else NA_character_,
    character(1)
  )
  epithet <- tolower(sub("\\.$", "", epithet))
  out[is.na(epithet) | epithet %in% placeholder_epithets] <- NA_character_
  out
}

# Query the GNA verifier for the lineage of `names` and cache it as a csv
# (`submitted_name`, `accepted_name`, one column per rank). Names are queried
# with spaces, whatever the authors' table uses ("Agaricus_essettei"). Run it again with
# force = TRUE when the strain list changes; the truth builder only reads the
# cache, so it never calls the network.
refresh_gna_lineage_cache <- function(
  names,
  cache_path,
  data_sources = 11,
  ranks = tolower(utils::head(truth_ranks, -1)),
  force = FALSE
) {
  if (!force && file.exists(cache_path)) {
    message("- mock truth: ", cache_path, " already exists; skipping")
    return(invisible(normalizePath(cache_path)))
  }
  res <- taxinfo::gna_verifier_pq(
    taxnames = unique(gsub("_", " ", names)),
    data_sources = data_sources,
    add_to_phyloseq = FALSE,
    classification_col = TRUE,
    classification_ranks = ranks,
    verbose = FALSE
  )
  rank_names <- paste0(toupper(substring(ranks, 1, 1)), substring(ranks, 2))
  rank_cols <- paste0("classification", rank_names)
  out <- tibble::tibble(
    submitted_name = res$submittedName,
    accepted_name = res$currentCanonicalSimple
  )
  out[rank_names] <- res[rank_cols]
  dir.create(dirname(cache_path), showWarnings = FALSE, recursive = TRUE)
  utils::write.csv(out, cache_path, row.names = FALSE)
  message(
    "- mock truth: ",
    nrow(out),
    " lineages written to ",
    cache_path,
    " (",
    sum(!is.na(out$Genus)),
    " with a genus)"
  )
  invisible(normalizePath(cache_path))
}

# Lineage of every strain of `sanger`, read from the cache of
# refresh_gna_lineage_cache(). The Genus and Species of the truth are the
# submitted name; `Genus_accepted` and `Species_accepted` hold the GBIF
# accepted name, which counts as correct too (hleap §6). Stops on a name
# missing from the cache rather than returning an empty lineage.
mock_lineage <- function(sanger, cache_path) {
  if (!file.exists(cache_path)) {
    stop(
      "mock_lineage(): ",
      cache_path,
      " is missing; run refresh_gna_lineage_cache() first."
    )
  }
  cache <- utils::read.csv(cache_path, check.names = FALSE)
  missing_names <- setdiff(
    unique(gsub("_", " ", sanger$name)),
    cache$submitted_name
  )
  if (length(missing_names) > 0) {
    stop(
      "mock_lineage(): ",
      length(missing_names),
      " names are missing from ",
      basename(cache_path),
      " (e.g. ",
      paste(utils::head(missing_names, 3), collapse = ", "),
      "); run refresh_gna_lineage_cache(force = TRUE)."
    )
  }
  rows <- match(gsub("_", " ", sanger$name), cache$submitted_name)
  above_genus <- utils::head(truth_ranks, -2)
  out <- tibble::tibble(strain = sanger$strain, name = sanger$name)
  out[above_genus] <- cache[rows, above_genus]
  out$Genus <- genus_name(sanger$name)
  out$Species <- species_name(sanger$name)
  out$Genus_accepted <- genus_name(cache$accepted_name[rows])
  out$Species_accepted <- species_name(cache$accepted_name[rows])
  out
}

# Ranks where the lineage of the authors' table (columns of `sanger_raw`)
# disagrees with the GBIF lineage of `lineage`, one row per disagreement. The
# decision of 2026-09-17 keeps the GBIF lineage and asks for these to be
# checked by hand.
lineage_disagreements <- function(lineage, sanger_raw, ranks = truth_ranks) {
  shared <- intersect(ranks, names(sanger_raw))
  shared <- intersect(shared, names(lineage))
  do.call(
    rbind,
    lapply(shared, function(rank) {
      published <- as.character(sanger_raw[[rank]])
      gbif <- as.character(lineage[[rank]])
      differs <- !is.na(published) &
        nzchar(published) &
        (is.na(gbif) | published != gbif)
      tibble::tibble(
        strain = lineage$strain[differs],
        rank = rank,
        published = published[differs],
        gbif = gbif[differs]
      )
    })
  )
}

# Strains whose Sanger sequence matches each unit of `seqs` (a named
# DNAStringSet: the refseq slot of the input being scored). vsearch keeps every
# hit at `search_id` or more, and a hit is a match when the identity reaches
# `min_identity` and the shorter of the two sequences is covered at
# `min_cover`; only the hits of highest identity of a unit are kept. Returns a
# tibble: `taxon`, `strain`, `identity`, `qcov`, `tcov`.
sanger_matches <- function(
  seqs,
  sanger,
  min_identity = 0.995,
  min_cover = 0.995,
  search_id = 0.9,
  nproc = 1
) {
  if (!dbpq::is_vsearch_installed()) {
    stop("sanger_matches() needs vsearch.")
  }
  query_file <- tempfile(fileext = ".fasta")
  db_file <- tempfile(fileext = ".fasta")
  out_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(query_file, db_file, out_file)), add = TRUE)
  Biostrings::writeXStringSet(seqs, query_file, width = 20000)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(stats::setNames(sanger$sequence, sanger$strain)),
    db_file,
    width = 20000
  )
  status <- system2(
    dbpq::find_vsearch(),
    c(
      "--usearch_global",
      shQuote(query_file),
      "--db",
      shQuote(db_file),
      "--id",
      format(search_id),
      "--maxaccepts",
      "0",
      "--maxrejects",
      "0",
      "--userfields",
      "query+target+id+qcov+tcov",
      "--userout",
      shQuote(out_file),
      "--threads",
      nproc,
      "--quiet"
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  if (status != 0) {
    stop("vsearch failed (exit ", status, ") while matching the Sanger reads.")
  }
  if (file.size(out_file) == 0) {
    return(tibble::tibble(
      taxon = character(),
      strain = character(),
      identity = numeric(),
      qcov = numeric(),
      tcov = numeric()
    ))
  }
  hits <- utils::read.delim(
    out_file,
    header = FALSE,
    quote = "",
    col.names = c("taxon", "strain", "identity", "qcov", "tcov")
  )
  hits |>
    tibble::as_tibble() |>
    dplyr::filter(
      identity >= 100 * min_identity,
      pmax(qcov, tcov) >= 100 * min_cover
    ) |>
    dplyr::group_by(taxon) |>
    dplyr::filter(identity == max(identity)) |>
    dplyr::ungroup()
}

# Truth of every unit of `seqs`, one row per unit (ROADMAP 0.6):
#
# - the rank columns hold the truth of the unit, NA below `truth_depth`;
# - `truth_depth` is the deepest rank with a truth: the strain's Species when
#   one strain matched, their lowest common rank when several did (their
#   identical lineages are collapsed first, so two strains of the same species
#   are not a tie), and NA when none did. A unit is scored down to
#   `truth_depth` and left out of the matrix below it;
# - `Genus_accepted` and `Species_accepted` hold the GBIF accepted name, which
#   counts as correct too (kept only when every matched strain agrees);
# - `n_strains` and `strains` record what matched.
mock_truth_table <- function(
  seqs,
  sanger,
  lineage,
  ranks = truth_ranks,
  ...
) {
  matches <- sanger_matches(seqs, sanger, ...)
  lineage <- lineage[match(sanger$strain, lineage$strain), ]
  rows <- lapply(names(seqs), function(taxon) {
    strains <- matches$strain[matches$taxon == taxon]
    truth <- lineage[lineage$strain %in% strains, , drop = FALSE]
    values <- stats::setNames(rep(NA_character_, length(ranks)), ranks)
    accepted <- c(
      Genus_accepted = NA_character_,
      Species_accepted = NA_character_
    )
    depth <- NA_character_
    if (nrow(truth) > 0) {
      for (rank in ranks) {
        shared <- unique(truth[[rank]])
        if (length(shared) != 1 || is.na(shared)) {
          break
        }
        values[[rank]] <- shared
        depth <- rank
      }
      for (col in names(accepted)) {
        shared <- unique(truth[[col]])
        rank <- sub("_accepted$", "", col)
        if (
          length(shared) == 1 &&
            !is.na(shared) &&
            !is.na(values[[rank]])
        ) {
          accepted[[col]] <- shared
        }
      }
    }
    tibble::tibble(
      taxon = taxon,
      !!!as.list(values),
      !!!as.list(accepted),
      truth_depth = depth,
      n_strains = length(strains),
      strains = paste(sort(strains), collapse = "|")
    )
  })
  dplyr::bind_rows(rows)
}

# Truth of the units of `physeq` (its refseq slot), from one row of
# `config.R::mock_sanger_tables`. Returns the truth table with the lineage
# disagreements of the authors' table attached as the attribute
# "disagreements" (empty when the table carries no lineage).
build_mock_truth <- function(physeq, cfg, base_dir = here::here(), ...) {
  sanger <- read_mock_sanger(cfg, base_dir = base_dir)
  lineage <- mock_lineage(sanger, file.path(base_dir, cfg$lineage_cache))
  truth <- mock_truth_table(physeq@refseq, sanger, lineage, ...)
  sanger_raw <- utils::read.csv(
    file.path(base_dir, cfg$path),
    sep = cfg$sep,
    check.names = FALSE
  )
  sanger_raw <- sanger_raw[
    match(sanger$strain, as.character(sanger_raw[[cfg$unit_col]])),
    ,
    drop = FALSE
  ]
  attr(truth, "disagreements") <- lineage_disagreements(lineage, sanger_raw)
  truth
}

# True lineage of the external negative controls, parsed from their taxon
# names. The records of the fake reference keep their UNITE general headers
# ("Name|accession|SH|reps|k__Alveolata;p__Ciliophora;...;s__Oxytrichidae_sp")
# and comparpq::add_external_seq_pq() only prefixes them with "external_", so
# the truth of a control is written in its name. Returns a tibble with a
# `taxon` column and one column per rank, ready for
# comparpq::tc_metrics_unit(external_truth = ). Species values are the
# binomials of the `Gen_sp_*` columns.
external_control_truth <- function(
  taxa,
  ranks = truth_ranks,
  pattern = "^external_"
) {
  taxa <- taxa[grepl(pattern, taxa)]
  letters_of <- c(
    Kingdom = "k",
    Phylum = "p",
    Class = "c",
    Order = "o",
    Family = "f",
    Genus = "g",
    Species = "s"
  )
  out <- tibble::tibble(taxon = taxa)
  for (rank in ranks) {
    field <- paste0(letters_of[[rank]], "__")
    value <- sub(paste0("^.*", field, "([^;|]*).*$"), "\\1", taxa)
    value[!grepl(field, taxa)] <- NA_character_
    out[[rank]] <- if (rank == "Species") species_name(value) else value
  }
  out
}
