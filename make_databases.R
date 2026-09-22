# Reproducible download and derivation of every reference database the
# benchmark consumes. Full documentation: docs/reference_databases.md.
#
# Two idempotent steps (existing outputs are kept unless force = TRUE):
#
#   1. download_reference_sources()
#      Fetches the general FASTA release of every row of `reference_sources`
#      (config.R) with dbpq::download_unite_db() or
#      dbpq::download_eukaryome_db(), stores it as
#      data/data_raw/refseq/sources/<source>.fasta and records its provenance
#      (URL, DOI, original file name, size, md5, date) in
#      data/data_raw/refseq/sources/manifest.csv.
#
#   2. derive_all_variants()
#      Converts each general FASTA to the dada2 and sintax header formats
#      (derive_dada2(), derive_sintax()), builds the simplified variants listed
#      in `benchmark_dbs` and `zoom_dbs` (kingdom Fungi only, Fungi plus
#      non-fungal representatives, primer-trimmed), the fake reference set and
#      the mini_* smoke-test subsets.
#
#   source("make_databases.R")
#   download_reference_sources()
#   derive_all_variants()
#
# To benchmark a new release: edit `reference_sources` and `benchmark_dbs` in
# config.R (new versioned names), run both steps, then the pipelines.
#
# External tools: bash, sed, awk, grep; cutadapt in the conda env of
# config.R::cutadapt_conda_prelude (primer-trimmed variants); vsearch for
# derive_fungi_rep() (controls kept away from the representatives).
#
# derive_clustered() and derive_no_parens() were deleted on 2026-09-22
# (ROADMAP C22-S2): both were defined but never called. Why, and what to do if
# a clustered reference ever returns to the design: proposals_for_dbpq.md
# sections 3 and 5.

library("here")
here::i_am("make_databases.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
load_pqverse("dbpq")

refseq_dir <- here("data/data_raw/refseq")
sources_dir <- file.path(refseq_dir, "sources")
dada2_dir <- file.path(refseq_dir, "dada2_format")
sintax_dir <- file.path(refseq_dir, "sintax_format")

# ---- helpers ----------------------------------------------------------------

# Returns TRUE (skip) when output already exists AND force is FALSE.
# When force is TRUE and the file exists, removes it and returns FALSE so the
# caller rebuilds.
skip_if_exists <- function(output, label = NULL, force = FALSE) {
  if (force && file.exists(output)) {
    file.remove(output)
    return(FALSE)
  }
  if (file.exists(output)) {
    msg <- if (is.null(label)) {
      sprintf("- %s already exists; skipping", output)
    } else {
      sprintf("- %s: %s already exists; skipping", label, output)
    }
    message(msg)
    return(TRUE)
  }
  FALSE
}

run_bash <- function(cmd) {
  status <- system(paste("bash -c", shQuote(cmd)))
  if (status != 0) {
    stop(sprintf("bash command failed (exit %d): %s", status, cmd))
  }
  invisible(status)
}

# Path of the reference file `name` (a source or a benchmarked database) in
# the general (downloaded), dada2 or sintax header format. Vectorised on name.
refseq_path <- function(
  name,
  format = c("dada2", "sintax", "general"),
  mini = FALSE
) {
  format <- match.arg(format)
  dir <- switch(
    format,
    dada2 = dada2_dir,
    sintax = sintax_dir,
    general = sources_dir
  )
  file.path(dir, paste0(if (mini) "mini_" else "", name, ".fasta"))
}

# ---- 1. download ------------------------------------------------------------

# Download the general FASTA release of one row of `reference_sources`, store
# it as <dir>/<source>.fasta and write <dir>/<source>.provenance.csv.
download_reference_source <- function(src, dir = sources_dir, force = FALSE) {
  output <- file.path(dir, paste0(src$source, ".fasta"))
  if (skip_if_exists(output, src$source, force = force)) {
    return(normalizePath(output))
  }
  work_dir <- file.path(dir, paste0(src$source, "_download"))
  dir.create(work_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)

  fasta <- switch(
    src$provider,
    unite = dbpq::download_unite_db(
      dest_dir = work_dir,
      url = src$url,
      extract = TRUE
    ),
    eukaryome = extract_single_fasta(
      dbpq::download_eukaryome_db(dest_dir = work_dir, url = src$url),
      exdir = work_dir
    ),
    stop("Unknown provider '", src$provider, "' for source ", src$source)
  )
  if (!file.rename(fasta, output)) {
    stop("Could not move ", fasta, " to ", output)
  }

  provenance <- data.frame(
    source = src$source,
    provider = src$provider,
    release = src$release,
    doi = src$doi,
    url = src$url,
    original_file = basename(fasta),
    bytes = file.size(output),
    md5 = unname(tools::md5sum(output)),
    downloaded_on = format(Sys.Date())
  )
  utils::write.csv(
    provenance,
    file.path(dir, paste0(src$source, ".provenance.csv")),
    row.names = FALSE
  )
  normalizePath(output)
}

# Unzip `archive` into `exdir` and return its single FASTA file. EUKARYOME
# v2.1 zips hold a .7z archive that holds the FASTA; nested .7z files are
# extracted with the 7z command-line tool (p7zip) and then deleted.
extract_single_fasta <- function(archive, exdir) {
  files <- utils::unzip(archive, exdir = exdir)
  for (nested in files[grepl("\\.7z$", files, ignore.case = TRUE)]) {
    seven_zip <- Sys.which(c("7z", "7za", "7zr"))
    seven_zip <- seven_zip[seven_zip != ""]
    if (length(seven_zip) == 0) {
      stop(
        "Extracting ",
        basename(nested),
        " needs the 7z command-line tool (p7zip) on PATH."
      )
    }
    status <- system2(
      seven_zip[[1]],
      c("x", "-y", paste0("-o", shQuote(exdir)), shQuote(nested)),
      stdout = FALSE
    )
    if (status != 0) {
      stop("7z failed (exit ", status, ") on ", basename(nested))
    }
    unlink(nested)
  }
  files <- list.files(exdir, full.names = TRUE, recursive = TRUE)
  fasta <- files[grepl("\\.(fasta|fas|fa|fna)$", files, ignore.case = TRUE)]
  if (length(fasta) != 1) {
    stop(
      "Expected one FASTA file in ",
      basename(archive),
      ", found: ",
      paste(basename(files), collapse = ", ")
    )
  }
  fasta
}

# Step 1 for every row of `sources`; (re)writes <dir>/manifest.csv from the
# per-source provenance files. Returns the general FASTA paths invisibly.
download_reference_sources <- function(
  sources = reference_sources,
  dir = sources_dir,
  force = FALSE
) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  paths <- vapply(
    seq_len(nrow(sources)),
    \(i) download_reference_source(sources[i, ], dir = dir, force = force),
    character(1)
  )
  provenance_files <- file.path(dir, paste0(sources$source, ".provenance.csv"))
  provenance_files <- provenance_files[file.exists(provenance_files)]
  if (length(provenance_files) > 0) {
    manifest <- do.call(rbind, lapply(provenance_files, utils::read.csv))
    utils::write.csv(
      manifest,
      file.path(dir, "manifest.csv"),
      row.names = FALSE
    )
  }
  invisible(stats::setNames(paths, sources$source))
}

# ---- 2. format conversion (delegated to dbpq) -------------------------------

# Header rewrites applied to a general FASTA before the dbpq conversion, on
# header lines only (`sed -E` rules below).
#
# - UNITE general-release headers ("Name|Acc|SH|type|k__Kingdom;p__...") put a
#   "|" right before the taxonomy. dbpq::format2sintax() and format2dada2()
#   read "...|k__Kingdom" as the identifier and drop the kingdom, which shifts
#   every rank by one (ROADMAP S6.2). "|k__" becomes ";k__".
# - EUKARYOME qualifies some genus, family and order names (ROADMAP S7.8):
#   homonyms "g__Lactarius(Fungi)", "f__Gonostomatidae(Sporadotrichida)";
#   bracketed genera "g__(Candida)" or "g__(Candida]"; suffixes ".s.str",
#   ".s.str.", ".nom.prov", ".nom.provis". Only the name is kept
#   ("g__Lactarius", "g__Candida", "g__Mortierella"). Left as is, the second
#   "(" of the sintax prediction "g:Lactarius(Fungi)(0.97)" makes
#   MiscMetabar::assign_sintax() read the bootstrap as NA (min_bootstrap then
#   never filters the genus), and the names never match the mock taxonomy.
#   Placeholders such as "Dothideales.gen05" or "Xxx.fam.incertae.sedis" are
#   kept.
general_header_rules <- paste0(
  "/^>/{",
  "s/\\|k__/;k__/;",
  "s/__\\(([^];)]*)[])]/__\\1/g;",
  "s/__([^;(]+)\\([^;)]*\\)/__\\1/g;",
  "s/__([^;]+)\\.(s\\.str|nom\\.prov(is)?)\\.?(;|$)/__\\1\\4/g",
  "}"
)

# Returns a copy of `input` with `general_header_rules` applied, written in
# `tmpdir`; the caller must delete it.
general_headers_fixed <- function(input, tmpdir = tempdir()) {
  fixed <- tempfile(fileext = ".fasta", tmpdir = tmpdir)
  run_bash(sprintf(
    "sed -E %s %s > %s",
    shQuote(general_header_rules),
    shQuote(input),
    shQuote(fixed)
  ))
  fixed
}

# General FASTA -> sintax headers (">ID;tax=k:Fungi,p:Ascomycota,...").
# The rewritten copy is written next to the output: the EUKARYOME ITS release
# (1.5 GB) may not fit a RAM-backed tempdir().
derive_sintax <- function(input, output, input_format = "auto", force = FALSE) {
  if (skip_if_exists(output, "sintax", force = force)) {
    return(normalizePath(output))
  }
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  source_file <- general_headers_fixed(input, tmpdir = dirname(output))
  on.exit(unlink(source_file), add = TRUE)
  dbpq::format2sintax(
    fasta_db = source_file,
    input_format = input_format,
    output_path = output
  )
  normalizePath(output)
}

# General FASTA -> dada2 pure-rank headers (">Fungi;Ascomycota;...;") expected
# by dada2::assignTaxonomy().
derive_dada2 <- function(input, output, input_format = "auto", force = FALSE) {
  if (skip_if_exists(output, "dada2", force = force)) {
    return(normalizePath(output))
  }
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  source_file <- general_headers_fixed(input, tmpdir = dirname(output))
  on.exit(unlink(source_file), add = TRUE)
  dbpq::format2dada2(
    fasta_db = source_file,
    input_format = input_format,
    output_path = output
  )
  normalizePath(output)
}

# ---- filtering -------------------------------------------------------------

# Header pattern that matches exactly the kingdom `kingdom` in a given header
# format. Anchoring on the kingdom field keeps "f:Fungiidae" (corals) or
# "k:cf.Fungi" out of a Fungi-only database: a plain "Fungi" pattern let 965
# non-Fungi records into the former EUK_ITS_v2_Fungi (ROADMAP S7.2).
kingdom_pattern <- function(
  format = c("sintax", "dada2", "general"),
  kingdom = "Fungi"
) {
  format <- match.arg(format)
  switch(
    format,
    sintax = paste0("tax=k:", kingdom, ","),
    dada2 = paste0("^>", kingdom, ";"),
    general = paste0("k__", kingdom, ";")
  )
}

# Keep the records whose kingdom is `kingdom`. Delegates to dbpq::filter_db
# (grep on headers), with two-line normalisation of multi-line records.
derive_kingdom_only <- function(
  input,
  output,
  format = c("sintax", "dada2", "general"),
  kingdom = "Fungi",
  force = FALSE
) {
  format <- match.arg(format)
  if (skip_if_exists(output, paste0(kingdom, "_only"), force = force)) {
    return(normalizePath(output))
  }
  dbpq::filter_db(
    ref_fasta = input,
    pattern = kingdom_pattern(format, kingdom),
    output = output,
    force_two_lines_per_seq = TRUE
  )
  normalizePath(output)
}

# Inverse pattern (remove records whose header matches `pattern`, an awk
# regular expression). No dbpq equivalent yet — proposal in
# proposals_for_dbpq.md.
derive_no_pattern <- function(input, output, pattern = "Fungi", force = FALSE) {
  if (skip_if_exists(output, paste0("no_", pattern), force = force)) {
    return(normalizePath(output))
  }
  tmp_2line <- tempfile(fileext = ".fasta")
  on.exit(unlink(tmp_2line), add = TRUE)
  # Normalize to two-line format so awk sees exactly one sequence line per header.
  run_bash(sprintf(
    "cat %s | sed ':a;N;/>/!s/\\n//;ta;P;D' > %s",
    shQuote(input),
    shQuote(tmp_2line)
  ))
  # Print only records whose header does NOT match pattern.
  run_bash(sprintf(
    "awk '/^>/{keep=!/%s/} keep{print}' %s > %s",
    pattern,
    shQuote(tmp_2line),
    shQuote(output)
  ))
  normalizePath(output)
}

# ---- subsetting ------------------------------------------------------------

# First `n` records as a quick smoke-test subset.
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md.)
derive_mini <- function(input, output, n = 10000, force = FALSE) {
  if (skip_if_exists(output, "mini", force = force)) {
    return(normalizePath(output))
  }
  run_bash(sprintf(
    "awk '/^>/{count++} count > %d {exit} {print}' %s > %s",
    n,
    shQuote(input),
    shQuote(output)
  ))
  normalizePath(output)
}

# Cutadapt-trim every record to the region between the primers (produces the
# `_cut` variants). Delegates to dbpq::cutadapt_rm_primers_db(), which
# reverse-complements primer_rev; `...` goes to it (e.g. min_overlap, nproc,
# args_before_cutadapt for the conda env). discard_untrimmed = TRUE (the dbpq
# default) drops the records without primer sites; derive_all_variants() keeps
# them (config.R::cut_discard_untrimmed, ROADMAP S6.9).
derive_cutadapted <- function(
  input,
  output,
  primer_fw,
  primer_rev,
  discard_untrimmed = TRUE,
  force = FALSE,
  ...
) {
  if (skip_if_exists(output, "cutadapted", force = force)) {
    return(normalizePath(output))
  }
  dbpq::cutadapt_rm_primers_db(
    ref_fasta = input,
    output = output,
    primer_fw = primer_fw,
    primer_rev = primer_rev,
    discard_untrimmed = discard_untrimmed,
    return_file_path = TRUE,
    ...
  )
  normalizePath(output)
}

# ---- fake reference set ----------------------------------------------------

# Number of records per kingdom of a sintax-format fasta, read from the
# headers only. Names go through `aliases` (e.g. EUKARYOME's Straminipila
# counted as Stramenopila).
kingdom_counts <- function(fasta, aliases = kingdom_aliases) {
  out <- system2(
    "grep",
    c("-o", shQuote("tax=k:[^,;]*"), shQuote(fasta)),
    stdout = TRUE
  )
  kingdoms <- sub("^tax=k:", "", out)
  aliased <- kingdoms %in% names(aliases)
  kingdoms[aliased] <- aliases[kingdoms[aliased]]
  c(table(kingdoms))
}

# Kingdoms kept as groups of representative sequences: at least `min_records`
# records, and not matching `exclude` (placeholder labels, Fungi).
retained_kingdoms <- function(
  counts,
  min_records,
  exclude = rep_kingdom_exclude
) {
  keep <- counts >= min_records & !grepl(exclude, names(counts))
  sort(names(counts)[keep])
}

# Kingdoms retained in every release of `min_records` (named by source): the
# pool of the external negative controls.
external_control_kingdoms <- function(min_records = rep_kingdom_min_records) {
  per_source <- lapply(names(min_records), function(s) {
    retained_kingdoms(
      kingdom_counts(refseq_path(s, "sintax")),
      min_records[[s]]
    )
  })
  Reduce(intersect, per_source)
}

# Build a balanced reference from `input` (non-Fungi records in general
# UNITE-style headers): one record per phylum, then random fill to `n`
# records. `kingdoms`, when given, restricts the pool to those kingdoms first.
# Used by add_external_seq_pq() in pipelines/assign_taxo.R to inject
# negative-control taxa. (Proposal `subset_balanced_db()` in
# proposals_for_dbpq.md.)
derive_fake_ref <- function(
  input,
  output,
  n = 100,
  seed = targets_seed,
  kingdoms = NULL,
  force = FALSE
) {
  if (skip_if_exists(output, "fake_ref", force = force)) {
    return(normalizePath(output))
  }
  set.seed(seed)
  seqs <- Biostrings::readDNAStringSet(input)
  if (!is.null(kingdoms)) {
    record_kingdom <- stringr::str_match(names(seqs), "k__([^;]+);")[, 2]
    seqs <- seqs[record_kingdom %in% kingdoms]
    if (length(seqs) == 0) {
      stop("derive_fake_ref(): no record of `input` belongs to `kingdoms`.")
    }
  }
  phyla <- stringr::str_match(names(seqs), "p__\\s*(.*?)\\s*;c__")[, 2]

  one_per_phylum <- vapply(
    unique(phyla),
    function(p) sample(which(phyla == p), 1),
    integer(1)
  )

  remaining <- setdiff(seq_along(seqs), one_per_phylum)
  # sample.int() rather than sample(): a length-one `remaining` would otherwise
  # be read as 1:remaining, and n larger than the pool would error.
  n_topup <- min(length(remaining), max(0, n - length(one_per_phylum)))
  topup <- remaining[sample.int(length(remaining), n_topup)]
  picks <- unique(c(one_per_phylum, topup))

  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs[picks], output, width = 10000)
  normalizePath(output)
}

# ---- fungi + rep databases -------------------------------------------------

# Number of records drawn from each kingdom: `total` records spread as evenly
# as possible over `available` (records per kingdom, named), by water-filling.
# Every kingdom gets the same cap, a kingdom with fewer records than the cap
# gives them all, and the records still missing go one each to the largest
# kingdoms (docs/objectives_design.md §1.3).
water_fill <- function(available, total) {
  if (total >= sum(available)) {
    return(available)
  }
  lo <- 0L
  hi <- max(available)
  while (lo < hi) {
    mid <- (lo + hi + 1L) %/% 2L
    if (sum(pmin(available, mid)) <= total) {
      lo <- mid
    } else {
      hi <- mid - 1L
    }
  }
  alloc <- pmin(available, lo)
  above <- names(available)[available > lo]
  above <- above[order(-available[above], above)]
  extra <- utils::head(above, total - sum(alloc))
  alloc[extra] <- alloc[extra] + 1L
  alloc
}

# Rows of a fasta index whose record has an identity of at least
# `max_identity` with one of the records of `external_fasta` (the negative
# controls). The indexed records are the queries and the controls the database:
# searching 100 controls against a whole release the other way round takes
# tens of minutes, this way a few minutes at most. Returns positions in `index`.
records_near_external <- function(index, external_fasta, max_identity) {
  if (!dbpq::is_vsearch_installed()) {
    stop(
      "derive_fungi_rep() needs vsearch to keep the representatives away from ",
      "the external controls (config.R::rep_max_identity_to_external)."
    )
  }
  queries_file <- tempfile(fileext = ".fasta")
  hits_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(queries_file, hits_file)), add = TRUE)
  Biostrings::writeXStringSet(
    Biostrings::readBStringSet(index),
    queries_file,
    width = 20000
  )
  status <- system2(
    dbpq::find_vsearch(),
    c(
      "--usearch_global",
      shQuote(queries_file),
      "--db",
      shQuote(external_fasta),
      "--id",
      format(max_identity),
      "--maxaccepts",
      "1",
      "--userfields",
      "query",
      "--userout",
      shQuote(hits_file),
      "--threads",
      "4",
      "--quiet"
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  if (status != 0) {
    stop("vsearch failed (exit ", status, ") while checking the controls.")
  }
  if (file.size(hits_file) == 0) {
    return(integer(0))
  }
  which(index$desc %in% readLines(hits_file))
}

# `fungi + rep` database of `source` (objectives_design decisions 3, 13, 14,
# 21): the records of its `_Fungi` database, then `share` x (number of fungal
# records) non-fungal records, rounded up, drawn from the full release of the
# same source. The draw is limited to the kingdoms kept by retained_kingdoms()
# (aliases applied), spread by water_fill(), and never takes a record whose
# sequence is one of the external negative controls (`external_fasta`): the
# releases do not share identifiers, so the match is on the sequence. The
# dada2 and sintax files of a source hold the same records in the same order
# (checked on the record lengths and kingdoms), so the same draw is written in
# both formats. A drawn record whose identity with a control reaches
# `max_identity_to_external` is dropped and redrawn (vsearch), so no control has
# a close relative in the database. The mini_* files hold the first `mini_n` records of the
# `_Fungi` database and the first `share` x `mini_n` representatives (drawn in
# random order). The drawn records are listed in
# sources/<output_db>_reps.csv.
derive_fungi_rep <- function(
  source,
  output_db = paste0(source, db_suffix[["Fungi+rep"]]),
  share = rep_share,
  min_records = rep_kingdom_min_records[[source]],
  exclude = rep_kingdom_exclude,
  aliases = kingdom_aliases,
  external_fasta = here(fake_ref_fasta),
  max_identity_to_external = rep_max_identity_to_external,
  seed = targets_seed,
  mini_n = 10000L,
  force = FALSE
) {
  formats <- c("sintax", "dada2")
  outputs <- vapply(formats, \(fmt) refseq_path(output_db, fmt), character(1))
  minis <- vapply(
    formats,
    \(fmt) refseq_path(output_db, fmt, mini = TRUE),
    character(1)
  )
  reps_csv <- file.path(sources_dir, paste0(output_db, "_reps.csv"))
  all_files <- c(outputs, minis, reps_csv)
  if (!force && all(file.exists(all_files))) {
    message("- fungi_rep: ", output_db, " already exists; skipping")
    return(invisible(normalizePath(outputs)))
  }
  unlink(all_files)

  index <- lapply(
    stats::setNames(formats, formats),
    \(fmt) Biostrings::fasta.index(refseq_path(source, fmt))
  )
  kingdom <- stringr::str_match(index$sintax$desc, "tax=k:([^,;]*)")[, 2]
  kingdom_dada2 <- sub(";.*$", "", index$dada2$desc)
  if (
    nrow(index$sintax) != nrow(index$dada2) ||
      !identical(index$sintax$seqlength, index$dada2$seqlength) ||
      !identical(dplyr::coalesce(kingdom, ""), kingdom_dada2)
  ) {
    stop(
      "derive_fungi_rep(): the sintax and dada2 files of ",
      source,
      " do not hold the same records in the same order."
    )
  }

  total <- ceiling(share * sum(kingdom == "Fungi", na.rm = TRUE))
  aliased <- kingdom
  has_alias <- aliased %in% names(aliases)
  aliased[has_alias] <- aliases[aliased[has_alias]]
  kingdoms <- retained_kingdoms(c(table(aliased)), min_records, exclude)
  candidate <- which(aliased %in% kingdoms)

  external <- toupper(as.character(Biostrings::readBStringSet(external_fasta)))
  same_length <- candidate[
    index$sintax$seqlength[candidate] %in% nchar(external)
  ]
  if (length(same_length) > 0) {
    seqs <- Biostrings::readBStringSet(index$sintax[same_length, ])
    stopifnot(identical(names(seqs), index$sintax$desc[same_length]))
    candidate <- setdiff(
      candidate,
      same_length[toupper(as.character(seqs)) %in% external]
    )
  }

  near <- records_near_external(
    index$sintax[candidate, ],
    external_fasta,
    max_identity_to_external
  )
  n_near <- length(near)
  if (n_near > 0) {
    candidate <- candidate[-near]
  }

  available <- c(table(factor(aliased[candidate], levels = kingdoms)))
  alloc <- water_fill(available, total)
  set.seed(seed)
  drawn <- unlist(lapply(kingdoms, function(k) {
    pool <- candidate[aliased[candidate] == k]
    pool[sample.int(length(pool), alloc[[k]])]
  }))
  drawn <- drawn[sample.int(length(drawn))]
  message(
    "- fungi_rep: ",
    output_db,
    ": ",
    length(drawn),
    " representatives from ",
    length(kingdoms),
    " kingdoms (",
    length(candidate) - length(drawn),
    " candidates left, ",
    n_near,
    " candidates too close to a control)"
  )

  n_mini_reps <- ceiling(share * mini_n)
  for (fmt in formats) {
    reps <- Biostrings::readBStringSet(index[[fmt]][drawn, ])
    stopifnot(identical(names(reps), index[[fmt]]$desc[drawn]))
    fungi_file <- refseq_path(paste0(source, "_Fungi"), fmt)
    if (!file.copy(fungi_file, outputs[[fmt]])) {
      stop("Could not copy ", fungi_file)
    }
    Biostrings::writeXStringSet(
      reps,
      outputs[[fmt]],
      append = TRUE,
      width = 20000
    )
    derive_mini(fungi_file, minis[[fmt]], n = mini_n, force = TRUE)
    Biostrings::writeXStringSet(
      reps[seq_len(min(length(reps), n_mini_reps))],
      minis[[fmt]],
      append = TRUE,
      width = 20000
    )
  }
  utils::write.csv(
    data.frame(
      record = sub(";tax=.*$", "", index$sintax$desc[drawn]),
      kingdom = aliased[drawn],
      header = index$sintax$desc[drawn]
    ),
    reps_csv,
    row.names = FALSE
  )
  invisible(normalizePath(outputs))
}

# ---- orchestration ---------------------------------------------------------

# Step 2: build every file the pipelines read, from the general FASTA files of
# step 1. Idempotent; rerunning after a failure resumes where it stopped.
# force = TRUE deletes and rebuilds each derived file.
derive_all_variants <- function(
  force = FALSE,
  sources = reference_sources,
  dbs = dplyr::bind_rows(benchmark_dbs, zoom_dbs)
) {
  missing <- sources$source[
    !file.exists(refseq_path(sources$source, "general"))
  ]
  if (length(missing) > 0) {
    stop(
      "General FASTA missing for: ",
      paste(missing, collapse = ", "),
      ". Run download_reference_sources() first."
    )
  }

  # 2a. General FASTA -> dada2 and sintax formats, one pair per source.
  for (s in sources$source) {
    derive_dada2(
      refseq_path(s, "general"),
      refseq_path(s, "dada2"),
      force = force
    )
    derive_sintax(
      refseq_path(s, "general"),
      refseq_path(s, "sintax"),
      force = force
    )
  }

  # 2b. Kingdom-Fungi and primer-trimmed variants, only where `dbs` uses them.
  fungi_sources <- unique(dbs$source[
    dbs$simplification %in% c("Fungi", "Fungi+cut", "Fungi+rep")
  ])
  cut_sources <- unique(dbs$source[dbs$simplification == "Fungi+cut"])
  for (s in fungi_sources) {
    for (fmt in c("dada2", "sintax")) {
      derive_kingdom_only(
        refseq_path(s, fmt),
        refseq_path(paste0(s, "_Fungi"), fmt),
        format = fmt,
        kingdom = "Fungi",
        force = force
      )
    }
  }
  for (s in cut_sources) {
    for (fmt in c("dada2", "sintax")) {
      derive_cutadapted(
        refseq_path(paste0(s, "_Fungi"), fmt),
        refseq_path(paste0(s, "_Fungi_cut"), fmt),
        primer_fw = fw_primer_sequences,
        primer_rev = rev_primer_sequences,
        discard_untrimmed = cut_discard_untrimmed,
        force = force,
        min_overlap = primer_min_overlap,
        min_length = cut_min_length,
        nproc = n_threads,
        args_before_cutadapt = cutadapt_conda_prelude
      )
    }
  }

  # 2b-bis. Amplicon-specific trimmed variants, one per biological dataset
  # (config.R::bio_datasets, developer decision 2026-09-16). Each `_cut`
  # reference is trimmed with EXACTLY the primers of the dataset that will be
  # assigned against it, so queries and reference cover the same region: ITS1
  # for the Pauvert mock above, the full ITS (ITS9mun / ITS4ngsUni) for
  # Tedersoo. min_overlap follows that dataset's own reverse primer, not the
  # ITS1 one of `primer_min_overlap`.
  for (i in seq_len(nrow(bio_datasets))) {
    ds <- bio_datasets[i, ]
    if (is.na(ds$cut_suffix)) {
      next
    }
    for (s in cut_sources) {
      for (fmt in c("dada2", "sintax")) {
        derive_cutadapted(
          refseq_path(paste0(s, "_Fungi"), fmt),
          refseq_path(paste0(s, ds$cut_suffix), fmt),
          primer_fw = ds$fw_primer,
          primer_rev = ds$rev_primer,
          discard_untrimmed = cut_discard_untrimmed,
          force = force,
          min_overlap = nchar(ds$rev_primer),
          min_length = cut_min_length,
          nproc = n_threads,
          args_before_cutadapt = cutadapt_conda_prelude
        )
        derive_mini(
          refseq_path(paste0(s, ds$cut_suffix), fmt),
          refseq_path(paste0(s, ds$cut_suffix), fmt, mini = TRUE),
          force = force
        )
      }
    }
  }

  # 2c. Fake reference: non-Fungi records of `fake_ref_source`.
  pool <- derive_no_pattern(
    refseq_path(fake_ref_source, "general"),
    file.path(sources_dir, paste0(fake_ref_source, "_wo_Fungi.fasta")),
    pattern = kingdom_pattern("general", "Fungi"),
    force = force
  )
  # Controls come only from the kingdoms retained in every release
  # (config.R::rep_kingdom_min_records); counting them reads every sintax
  # header, so it only runs when the file is rebuilt.
  if (force || !file.exists(here(fake_ref_fasta))) {
    derive_fake_ref(
      pool,
      here(fake_ref_fasta),
      seed = targets_seed,
      kingdoms = external_control_kingdoms(),
      force = force
    )
  }

  # 2c-bis. Fungi + non-fungal representatives, after the fake reference,
  # whose records must never be drawn as representatives. derive_fungi_rep()
  # writes their mini_* files too.
  for (s in unique(dbs$source[dbs$simplification == "Fungi+rep"])) {
    derive_fungi_rep(s, force = force)
  }

  # 2d. mini_* subsets (first 10 000 records) of every other database.
  for (db in dbs$db[dbs$simplification != "Fungi+rep"]) {
    for (fmt in c("dada2", "sintax")) {
      derive_mini(
        refseq_path(db, fmt),
        refseq_path(db, fmt, mini = TRUE),
        force = force
      )
    }
  }

  invisible(NULL)
}
