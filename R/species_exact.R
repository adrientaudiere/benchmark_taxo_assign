# Species by exact match (OD 43): the species section names a species only
# when a query is found at 100 % in a reference record. Three methods, each
# kept as raw output so that the rule turning hits into a species call is a
# scoring choice, not a rerun:
#
# - dada2::assignSpecies() on the species_format/ files
#   (make_databases.R::derive_species()): a query is an exact substring of a
#   record. Its own output rules (dada2:::mapHits()): Genus is kept only when
#   every hit has the same genus, else NA; Species holds the EPITHETS of all
#   hits, sorted and joined by "/" (allowMultiple = TRUE), without their genus.
# - vsearch --usearch_global --id 1 --query_cov 1 --maxaccepts 0
#   --maxrejects 32 on the sintax files: every (query, record) pair, then
#   kept only when the query is an exact substring of the record, because
#   vsearch counts a base matching an ambiguity code (N, Y) of the record as
#   identical. --search_exact is not used: it reports full-length identity
#   only (analysis/sandbox/check_vsearch_exact_species.R).
# - blastn: a search of its own (blastn_near_hits(), since 2026-10-07; the
#   grid's raw tables before), the pairs with no difference over the query.
#
# vsearch and blastn search every record of a database, also the ones that
# name no species ("<Genus>_sp", unclassified): their binomial is NA, so a
# scoring rule can decide whether such a hit blocks a species call.
#
# Requires Biostrings, dplyr and stringi; vsearch through dbpq::find_vsearch(),
# blastn and makeblastdb on PATH.

# "Genus_epithet" of sintax headers ("...;tax=k:...,g:Genus,s:..."), NA when
# the record names no species. UNITE writes "s:Genus_epithet", EUKARYOME the
# epithet alone ("s:epithet", its genus in "g:"). Vectorised.
binomial_of_sintax_header <- function(headers) {
  genus <- sub("^.*[;,]g:([^,;]*).*$", "\\1", headers)
  genus[!grepl("[;,]g:", headers)] <- NA_character_
  species <- sub("^.*,s:([^,;]*).*$", "\\1", headers)
  species[!grepl(",s:", headers)] <- NA_character_
  binomial <- ifelse(
    grepl("_", species),
    species,
    paste0(genus, "_", species)
  )
  named <- !is.na(species) &
    grepl("^[A-Z][a-z]+_[a-z][a-z-]+$", binomial) &
    !grepl("_sp$|unclassified|incertae", binomial, ignore.case = TRUE)
  binomial[!named] <- NA_character_
  binomial
}

# The queries of several inputs pooled into one set of distinct sequences,
# named q1..qN, and the table that maps each back to every (dataset, input,
# taxon) holding it. `inputs` is a named list of DNAStringSet (names
# "<dataset>|<input>"). Sequences with a base other than A, C, G or T are left
# out (dada2::assignSpecies() refuses them) and counted in the `dropped`
# attribute of the map, so that the three methods see the same queries.
pool_queries <- function(inputs) {
  map <- dplyr::bind_rows(lapply(names(inputs), function(key) {
    s <- inputs[[key]]
    tibble::tibble(
      dataset = sub("\\|.*$", "", key),
      input = sub("^.*\\|", "", key),
      taxon = names(s),
      sequence = as.character(s)
    )
  }))
  acgt <- !grepl("[^ACGT]", map$sequence)
  dropped <- map[!acgt, c("dataset", "input", "taxon")]
  map <- map[acgt, ]
  distinct_seqs <- unique(map$sequence)
  map$query <- paste0("q", match(map$sequence, distinct_seqs))
  seqs <- Biostrings::DNAStringSet(distinct_seqs)
  names(seqs) <- paste0("q", seq_along(seqs))
  map$sequence <- NULL
  attr(map, "dropped") <- dropped
  list(seqs = seqs, map = map)
}

# dada2::assignSpecies() as a tibble: one row per query of `seqs`, its Genus
# (NA when the hits span several genera) and Species (epithets joined by "/").
assign_species_exact <- function(seqs, species_fasta) {
  out <- dada2::assignSpecies(
    as.character(seqs),
    species_fasta,
    allowMultiple = TRUE,
    tryRC = FALSE
  )
  tibble::tibble(
    query = names(seqs),
    genus = unname(out[, "Genus"]),
    epithets = unname(out[, "Species"])
  )
}

# Every record of `db_fasta` (sintax headers) holding a query of `seqs` at
# 100 %: one row per (query, record) with the record's binomial (NA when it
# names no species). --notrunclabels keeps the whole header (vsearch cuts
# labels at the first space otherwise; 3 EUKARYOME headers hold one).
vsearch_exact_hits <- function(seqs, db_fasta, nproc = 1) {
  query_file <- tempfile(fileext = ".fasta")
  out_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(query_file, out_file)), add = TRUE)
  Biostrings::writeXStringSet(seqs, query_file, width = 20000)
  status <- system2(
    dbpq::find_vsearch(),
    c(
      "--usearch_global",
      shQuote(query_file),
      "--db",
      shQuote(db_fasta),
      "--id",
      "1",
      "--query_cov",
      "1",
      "--maxaccepts",
      "0",
      "--maxrejects",
      "32",
      "--strand",
      "plus",
      "--notrunclabels",
      "--threads",
      nproc,
      "--userout",
      shQuote(out_file),
      "--userfields",
      "query+target",
      "--quiet"
    )
  )
  if (status != 0) {
    stop("vsearch failed (exit ", status, ") on ", basename(db_fasta))
  }
  hits <- if (file.size(out_file) > 0) {
    utils::read.delim(
      out_file,
      header = FALSE,
      col.names = c("query", "record"),
      quote = "",
      colClasses = "character"
    )
  } else {
    data.frame(query = character(), record = character())
  }
  if (nrow(hits) == 0) {
    return(tibble::tibble(
      query = character(),
      record = character(),
      binomial = character()
    ))
  }
  ref <- Biostrings::readDNAStringSet(db_fasta)
  ref <- ref[unique(hits$record)]
  exact <- mapply(
    \(q, r) grepl(q, r, fixed = TRUE),
    as.character(seqs[hits$query]),
    as.character(ref[hits$record]),
    USE.NAMES = FALSE
  )
  hits <- hits[exact, ]
  tibble::tibble(
    query = hits$query,
    record = hits$record,
    binomial = binomial_of_sintax_header(hits$record)
  )
}

# ---- Near matches (OD 43, developer 2026-10-07) ------------------------------
# Two relaxations of the exact match, read by vsearch and blastn
# (assignSpecies() has no such mode): at most one difference over the whole
# query, and 100 % identity over at least 90 % of the query, the rest lying
# beyond an end of the record (the record stops early). Each search only proposes (query, record) pairs; whether a pair
# passes is then decided from the two sequences, identically for both
# methods, so the methods differ only by the pairs their search finds.

# Differences over the whole query between each `query` and its `record`,
# given the stretches a search aligned (qlo..qhi on the query, tlo..thi on the
# record): the edit distance between the query and the record stretch
# extended by the unaligned ends of the query, plus the query bases that fall
# beyond an end of the record. 0 when the query is an exact substring of the
# record; an ambiguity code of the record (N, Y) counts as a difference.
# `approx` is the count the search reports; pairs above `max_diffs` + 2 keep
# it, so the edit distance is computed only where it can matter. Vectorised.
query_diffs <- function(query, record, qlo, qhi, tlo, thi, approx, max_diffs) {
  ql <- nchar(query)
  rl <- nchar(record)
  start <- tlo - (qlo - 1)
  end <- thi + (ql - qhi)
  before <- pmax(0, 1 - start)
  after <- pmax(0, end - rl)
  q_part <- substring(query, 1 + before, ql - after)
  r_part <- substring(record, pmax(1, start), pmin(rl, end))
  out <- as.numeric(approx)
  exact <- stringi::stri_detect_fixed(record, query)
  todo <- which(!exact & approx <= max_diffs + 2)
  out[todo] <- vapply(
    todo,
    \(i) drop(utils::adist(q_part[i], r_part[i])),
    numeric(1)
  ) +
    before[todo] +
    after[todo]
  out[exact] <- 0
  as.integer(out)
}

# TRUE when the record holds the query at 100 % wherever it covers it, covers
# at least `min_cover` of it, and only stops early: the rest of the query lies
# beyond an end of the record, not against different bases (developer,
# 2026-10-07; the record of a species that ends before the ITSx region does).
# Four cases: the query is inside the record; a prefix of the query ends the
# record; a suffix of the query starts the record; the record lies inside the
# query. Vectorised.
exact_cover <- function(query, record, min_cover) {
  ql <- nchar(query)
  rl <- nchar(record)
  k_min <- ceiling(min_cover * ql)
  found <- stringi::stri_detect_fixed(record, query) |
    (rl >= k_min & rl < ql & stringi::stri_detect_fixed(query, record))
  if (length(query) == 0) {
    return(found)
  }
  for (d in seq_len(max(c(0, ql - k_min)))) {
    todo <- !found & d <= ql - k_min
    if (!any(todo)) {
      break
    }
    q <- query[todo]
    n <- ql[todo]
    found[todo] <- stringi::stri_endswith_fixed(
      record[todo],
      substring(q, 1, n - d)
    ) |
      stringi::stri_startswith_fixed(record[todo], substring(q, d + 1, n))
  }
  found
}

# The near-match columns of candidate pairs: `diffs` (query_diffs()) and
# `exact_cover` (exact_cover()), keeping the pairs within `max_diffs` or
# covered at 100 % over `min_cover` of the query. `pairs` holds query, record, qlo, qhi, tlo, thi, approx.
near_hit_columns <- function(pairs, seqs, ref, max_diffs, min_cover) {
  query_seq <- as.character(seqs[pairs$query])
  record_seq <- as.character(ref[pairs$record])
  pairs$diffs <- query_diffs(
    query_seq,
    record_seq,
    pairs$qlo,
    pairs$qhi,
    pairs$tlo,
    pairs$thi,
    pairs$approx,
    max_diffs
  )
  pairs$exact_cover <- exact_cover(query_seq, record_seq, min_cover)
  pairs <- pairs[pairs$diffs <= max_diffs | pairs$exact_cover, ]
  tibble::tibble(
    query = pairs$query,
    record = pairs$record,
    binomial = binomial_of_sintax_header(pairs$record),
    diffs = pairs$diffs,
    exact_cover = pairs$exact_cover
  )
}

empty_near_hits <- function() {
  tibble::tibble(
    query = character(),
    record = character(),
    binomial = character(),
    diffs = integer(),
    exact_cover = logical()
  )
}

# vsearch candidates for the near matches: --usearch_global at `id` over the
# aligned part (terminal gaps excluded, vsearch's default) with `query_cov`,
# every accepted record (--maxaccepts 0, --maxrejects 32 as for the exact
# search). id 0.9 also proposes a query whose last 10 % differ inside the
# record; on the Rolling queries against full EUK_ITS_v2.1 it gives about
# 10 000 pairs in 3 min (2026-10-07).
vsearch_near_hits <- function(
  seqs,
  db_fasta,
  nproc = 1,
  id = 0.9,
  query_cov = 0.9,
  max_diffs = 1L,
  min_cover = 0.9
) {
  query_file <- tempfile(fileext = ".fasta")
  out_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(query_file, out_file)), add = TRUE)
  Biostrings::writeXStringSet(seqs, query_file, width = 20000)
  status <- system2(
    dbpq::find_vsearch(),
    c(
      "--usearch_global",
      shQuote(query_file),
      "--db",
      shQuote(db_fasta),
      "--id",
      id,
      "--query_cov",
      query_cov,
      "--maxaccepts",
      "0",
      "--maxrejects",
      "32",
      "--strand",
      "plus",
      "--notrunclabels",
      "--threads",
      nproc,
      "--userout",
      shQuote(out_file),
      "--userfields",
      "query+target+ql+qilo+qihi+tilo+tihi+mism+gaps",
      "--quiet"
    )
  )
  if (status != 0) {
    stop("vsearch failed (exit ", status, ") on ", basename(db_fasta))
  }
  if (file.size(out_file) == 0) {
    return(empty_near_hits())
  }
  hits <- utils::read.delim(
    out_file,
    header = FALSE,
    col.names = c(
      "query",
      "record",
      "ql",
      "qlo",
      "qhi",
      "tlo",
      "thi",
      "mism",
      "gaps"
    ),
    quote = "",
    colClasses = c("character", "character", rep("integer", 7))
  )
  # vsearch counts no terminal gap in `gaps`: the query bases beyond the
  # aligned part are added here.
  hits$approx <- hits$mism + hits$gaps + (hits$qlo - 1) + (hits$ql - hits$qhi)
  ref <- Biostrings::readDNAStringSet(db_fasta)
  near_hit_columns(hits, seqs, ref, max_diffs, min_cover)
}

# blastn candidates for the exact and near matches, from a search of the
# species layer itself rather than the grid's tables, which keep 500 hits per
# query and no alignment coordinates (developer 2026-10-07): megablast
# (blastn_search_args of config.R for the rest), HSPs at >= `perc_identity` over
# >= `qcov_hsp` % of the query, up to `max_target_seqs` records per query.
# One row per (query, record) with the fewest differences over its HSPs, and
# `n_targets`, the records blastn returned for the query: a query at
# `max_target_seqs` may lack some.
blastn_near_hits <- function(
  seqs,
  db_fasta,
  nproc = 1,
  perc_identity = 90,
  qcov_hsp = 89,
  max_target_seqs = 2000L,
  max_diffs = 1L,
  min_cover = 0.9
) {
  dir <- tempfile("blastn_near_")
  dir.create(dir)
  on.exit(unlink(dir, recursive = TRUE), add = TRUE)
  query_file <- file.path(dir, "queries.fasta")
  out_file <- file.path(dir, "hits.tsv")
  Biostrings::writeXStringSet(seqs, query_file, width = 20000)
  # The records are indexed under their position (r1, r2...): BLAST+ rewrites
  # some ids (it drops the trailing "." of a EUKARYOME header ending in
  # "nov.inval."), so its ids cannot be matched back to the headers.
  ref <- Biostrings::readDNAStringSet(db_fasta)
  indexed <- ref
  names(indexed) <- paste0("r", seq_along(ref))
  ref_file <- file.path(dir, "ref.fasta")
  Biostrings::writeXStringSet(indexed, ref_file, width = 20000)
  rm(indexed)
  status <- system2(
    "makeblastdb",
    c(
      "-dbtype nucl -in",
      shQuote(ref_file),
      "-out",
      shQuote(file.path(dir, "db"))
    ),
    stdout = FALSE
  )
  if (status != 0) {
    stop("makeblastdb failed (exit ", status, ") on ", basename(db_fasta))
  }
  status <- system2(
    "blastn",
    c(
      "-query",
      shQuote(query_file),
      "-db",
      shQuote(file.path(dir, "db")),
      "-out",
      shQuote(out_file),
      "-strand plus",
      "-perc_identity",
      perc_identity,
      "-qcov_hsp_perc",
      qcov_hsp,
      "-num_threads",
      nproc,
      blastn_search_args,
      "-max_target_seqs",
      max_target_seqs,
      "-outfmt",
      shQuote("6 qseqid sseqid qlen length nident qstart qend sstart send")
    )
  )
  if (status != 0) {
    stop("blastn failed (exit ", status, ") on ", basename(db_fasta))
  }
  if (file.size(out_file) == 0) {
    return(dplyr::mutate(empty_near_hits(), n_targets = integer()))
  }
  hits <- utils::read.delim(
    out_file,
    header = FALSE,
    col.names = c(
      "query",
      "sseqid",
      "ql",
      "length",
      "nident",
      "qlo",
      "qhi",
      "tlo",
      "thi"
    ),
    quote = "",
    colClasses = c("character", "character", rep("integer", 7))
  )
  # Back to the whole header, as vsearch reports it.
  hits$record <- names(ref)[as.integer(sub("^r", "", hits$sseqid))]
  if (anyNA(hits$record)) {
    stop(
      "blastn returned ids that are not record positions: ",
      basename(db_fasta)
    )
  }
  hits$approx <- (hits$length - hits$nident) +
    (hits$qlo - 1) +
    (hits$ql - hits$qhi)
  n_targets <- tapply(hits$sseqid, hits$query, \(x) length(unique(x)))
  near_hit_columns(hits, seqs, ref, max_diffs, min_cover) |>
    dplyr::group_by(.data$query, .data$record, .data$binomial) |>
    dplyr::summarise(
      diffs = min(.data$diffs),
      exact_cover = any(.data$exact_cover),
      .groups = "drop"
    ) |>
    dplyr::mutate(n_targets = as.integer(n_targets[.data$query]))
}
