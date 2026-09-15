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
#      in `benchmark_dbs` (kingdom Fungi only, then primer-trimmed), the fake
#      reference set and the mini_* smoke-test subsets.
#
#   source("make_databases.R")
#   download_reference_sources()
#   derive_all_variants()
#
# To benchmark a new release: edit `reference_sources` and `benchmark_dbs` in
# config.R (new versioned names), run both steps, then the pipelines.
#
# External tools: bash, sed, awk, grep; cutadapt in the conda env of
# config.R::cutadapt_conda_prelude (primer-trimmed variants); vsearch only for
# derive_clustered(), which the orchestration does not call.

library("here")
here::i_am("make_databases.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
load_pqverse("dbpq")

refseq_dir  <- here("data/data_raw/refseq")
sources_dir <- file.path(refseq_dir, "sources")
dada2_dir   <- file.path(refseq_dir, "dada2_format")
sintax_dir  <- file.path(refseq_dir, "sintax_format")

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
refseq_path <- function(name, format = c("dada2", "sintax", "general"), mini = FALSE) {
  format <- match.arg(format)
  dir <- switch(format, dada2 = dada2_dir, sintax = sintax_dir, general = sources_dir)
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
    source        = src$source,
    provider      = src$provider,
    release       = src$release,
    doi           = src$doi,
    url           = src$url,
    original_file = basename(fasta),
    bytes         = file.size(output),
    md5           = unname(tools::md5sum(output)),
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
      stop("Extracting ", basename(nested), " needs the 7z command-line tool (p7zip) on PATH.")
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
      "Expected one FASTA file in ", basename(archive), ", found: ",
      paste(basename(files), collapse = ", ")
    )
  }
  fasta
}

# Step 1 for every row of `sources`; (re)writes <dir>/manifest.csv from the
# per-source provenance files. Returns the general FASTA paths invisibly.
download_reference_sources <- function(sources = reference_sources,
                                       dir = sources_dir,
                                       force = FALSE) {
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
    utils::write.csv(manifest, file.path(dir, "manifest.csv"), row.names = FALSE)
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
    shQuote(general_header_rules), shQuote(input), shQuote(fixed)
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

# Strip parenthesized synonyms found in some EUKARYOME headers.
# Superseded by general_headers_fixed(), which keeps "g__(Candida)" as
# "g__Candida" instead of emptying it. Not called by the orchestration.
derive_no_parens <- function(input, output, force = FALSE) {
  if (skip_if_exists(output, "no_parens", force = force)) {
    return(normalizePath(output))
  }
  run_bash(sprintf('sed "s/([^)]*)//g" %s > %s',
                   shQuote(input), shQuote(output)))
  normalizePath(output)
}

# ---- filtering -------------------------------------------------------------

# Header pattern that matches exactly the kingdom `kingdom` in a given header
# format. Anchoring on the kingdom field keeps "f:Fungiidae" (corals) or
# "k:cf.Fungi" out of a Fungi-only database: a plain "Fungi" pattern let 965
# non-Fungi records into the former EUK_ITS_v2_Fungi (ROADMAP S7.2).
kingdom_pattern <- function(format = c("sintax", "dada2", "general"), kingdom = "Fungi") {
  format <- match.arg(format)
  switch(
    format,
    sintax  = paste0("tax=k:", kingdom, ","),
    dada2   = paste0("^>", kingdom, ";"),
    general = paste0("k__", kingdom, ";")
  )
}

# Keep the records whose kingdom is `kingdom`. Delegates to dbpq::filter_db
# (grep on headers), with two-line normalisation of multi-line records.
derive_kingdom_only <- function(input, output,
                                format = c("sintax", "dada2", "general"),
                                kingdom = "Fungi",
                                force = FALSE) {
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
    shQuote(input), shQuote(tmp_2line)
  ))
  # Print only records whose header does NOT match pattern.
  run_bash(sprintf(
    "awk '/^>/{keep=!/%s/} keep{print}' %s > %s",
    pattern, shQuote(tmp_2line), shQuote(output)
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
    n, shQuote(input), shQuote(output)
  ))
  normalizePath(output)
}

# Cluster at the given identity threshold via vsearch.
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md. Not called by
# the orchestration.)
derive_clustered <- function(input, output, identity = 0.99, force = FALSE) {
  if (skip_if_exists(output, "clustered", force = force)) {
    return(normalizePath(output))
  }
  if (!dbpq::is_vsearch_installed()) {
    stop("vsearch is required for derive_clustered() but was not found on PATH.")
  }
  status <- system2(
    dbpq::find_vsearch(),
    args = c("--cluster_fast", input,
             "--id", format(identity, nsmall = 2),
             "--centroids", output)
  )
  if (status != 0) {
    stop(sprintf("vsearch failed (exit %d) on %s", status, input))
  }
  normalizePath(output)
}

# Cutadapt-trim every record to the region between the primers (produces the
# `_cut` variants). Delegates to dbpq::cutadapt_rm_primers_db(), which
# reverse-complements primer_rev; `...` goes to it (e.g. min_overlap, nproc,
# args_before_cutadapt for the conda env). discard_untrimmed = TRUE (the dbpq
# default) drops the records without primer sites; derive_all_variants() keeps
# them (config.R::cut_discard_untrimmed, ROADMAP S6.9).
derive_cutadapted <- function(input, output, primer_fw, primer_rev,
                              discard_untrimmed = TRUE, force = FALSE, ...) {
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

# Build a balanced reference from `input` (non-Fungi records in general
# UNITE-style headers): one record per phylum, then random fill to `n`
# records. Used by add_external_seq_pq() in pipelines/assign_taxo.R to inject
# negative-control taxa. (Proposal `subset_balanced_db()` in
# proposals_for_dbpq.md.)
derive_fake_ref <- function(input, output, n = 100, seed = targets_seed, force = FALSE) {
  if (skip_if_exists(output, "fake_ref", force = force)) {
    return(normalizePath(output))
  }
  set.seed(seed)
  seqs  <- Biostrings::readDNAStringSet(input)
  phyla <- stringr::str_match(names(seqs), "p__\\s*(.*?)\\s*;c__")[, 2]

  one_per_phylum <- vapply(
    unique(phyla),
    function(p) sample(which(phyla == p), 1),
    integer(1)
  )

  remaining <- setdiff(seq_along(seqs), one_per_phylum)
  # sample.int() rather than sample(): a length-one `remaining` would otherwise
  # be read as 1:remaining, and n larger than the pool would error.
  n_topup   <- min(length(remaining), max(0, n - length(one_per_phylum)))
  topup     <- remaining[sample.int(length(remaining), n_topup)]
  picks     <- unique(c(one_per_phylum, topup))

  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs[picks], output, width = 10000)
  normalizePath(output)
}

# ---- orchestration ---------------------------------------------------------

# Step 2: build every file the pipelines read, from the general FASTA files of
# step 1. Idempotent; rerunning after a failure resumes where it stopped.
# force = TRUE deletes and rebuilds each derived file.
derive_all_variants <- function(force = FALSE,
                                sources = reference_sources,
                                dbs = benchmark_dbs) {
  missing <- sources$source[!file.exists(refseq_path(sources$source, "general"))]
  if (length(missing) > 0) {
    stop(
      "General FASTA missing for: ", paste(missing, collapse = ", "),
      ". Run download_reference_sources() first."
    )
  }

  # 2a. General FASTA -> dada2 and sintax formats, one pair per source.
  for (s in sources$source) {
    derive_dada2(refseq_path(s, "general"), refseq_path(s, "dada2"), force = force)
    derive_sintax(refseq_path(s, "general"), refseq_path(s, "sintax"), force = force)
  }

  # 2b. Kingdom-Fungi and primer-trimmed variants, only where `dbs` uses them.
  fungi_sources <- unique(dbs$source[dbs$simplification %in% c("Fungi", "Fungi+cut")])
  cut_sources   <- unique(dbs$source[dbs$simplification == "Fungi+cut"])
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

  # 2c. Fake reference: non-Fungi records of `fake_ref_source`.
  pool <- derive_no_pattern(
    refseq_path(fake_ref_source, "general"),
    file.path(sources_dir, paste0(fake_ref_source, "_wo_Fungi.fasta")),
    pattern = kingdom_pattern("general", "Fungi"),
    force = force
  )
  derive_fake_ref(pool, here(fake_ref_fasta), seed = targets_seed, force = force)

  # 2d. mini_* subsets (first 10 000 records) of every benchmarked database.
  for (db in dbs$db) {
    for (fmt in c("dada2", "sintax")) {
      derive_mini(refseq_path(db, fmt), refseq_path(db, fmt, mini = TRUE), force = force)
    }
  }

  invisible(NULL)
}
