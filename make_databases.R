# Reproducible derivation of every reference-database variant the benchmark
# consumes. Replaces the historical `some_bash_script` (now in Archives/).
#
# Most operations are delegated to `dbpq::` — only the steps for which dbpq
# does not (yet) have an equivalent stay local. Those are flagged in
# proposals_for_dbpq.md and should migrate to dbpq when added.
#
# Every function is idempotent: if the output file already exists, the
# function returns its path without rerunning. Delete the output first to
# force a rebuild.
#
# Sourcing this file has no side effects. Build everything with:
#   source("make_databases.R")
#   derive_all_variants()
#
# External tools required: bash, sed, grep, head, vsearch, plus cutadapt in
# the `cutadaptenv` conda env (see config.R::cutadapt_conda_prelude).

library("here")
library("dbpq")
here::i_am("make_databases.R")

# ---- helpers ----------------------------------------------------------------

skip_if_exists <- function(output, label = NULL) {
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
  status <- system2("bash", args = c("-c", cmd))
  if (status != 0) {
    stop(sprintf("bash command failed (exit %d): %s", status, cmd))
  }
  invisible(status)
}

# ---- format conversion (delegated to dbpq) ---------------------------------

# dbpq::format2sintax() auto-detects the input format (UNITE/dada2,
# Greengenes2, plain) and writes the sintax-tagged fasta. Subsumes both the
# legacy `derive_sintax_from_dada2` (sed pipeline) and `derive_sintax_from_plain`
# (>;tax= prefix) helpers.
derive_sintax <- function(input, output, input_format = "auto") {
  if (skip_if_exists(output, "sintax")) {
    return(normalizePath(output))
  }
  dbpq::format2sintax(
    fasta_db = input,
    input_format = input_format,
    output_path = output
  )
  normalizePath(output)
}

# Strip parenthesized synonyms found in EUKARYOME v1.9.3 headers.
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md.)
derive_no_parens <- function(input, output) {
  if (skip_if_exists(output, "no_parens")) {
    return(normalizePath(output))
  }
  run_bash(sprintf('sed "s/([^)]*)//g" %s > %s',
                   shQuote(input), shQuote(output)))
  normalizePath(output)
}

# ---- filtering -------------------------------------------------------------

# Keep only records whose header matches `pattern`. Delegates to
# dbpq::filter_db, which handles both two-line and multi-line fastas via
# `force_two_lines_per_seq`.
derive_pattern_only <- function(input, output,
                                pattern = "Fungi",
                                multi_line_records = FALSE) {
  if (skip_if_exists(output, paste0(pattern, "_only"))) {
    return(normalizePath(output))
  }
  dbpq::filter_db(
    ref_fasta = input,
    pattern = pattern,
    output = output,
    force_two_lines_per_seq = multi_line_records
  )
  normalizePath(output)
}

# Inverse pattern (remove records whose header matches). No dbpq equivalent
# yet — proposal in proposals_for_dbpq.md. `grep -vFf -` reads the matching
# block on stdin and excludes each line as a fixed string.
derive_no_pattern <- function(input, output, pattern = "Fungi") {
  if (skip_if_exists(output, paste0("no_", pattern))) {
    return(normalizePath(output))
  }
  run_bash(sprintf(
    "grep -vFf - %s < <(grep -A1 '%s' %s) > %s",
    shQuote(input), pattern, shQuote(input), shQuote(output)
  ))
  normalizePath(output)
}

# ---- subsetting ------------------------------------------------------------

# First `n` lines as a quick smoke-test subset.
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md.)
derive_mini <- function(input, output, n = 10000) {
  if (skip_if_exists(output, "mini")) {
    return(normalizePath(output))
  }
  run_bash(sprintf("head -n %d %s > %s",
                   n, shQuote(input), shQuote(output)))
  normalizePath(output)
}

# Cluster at the given identity threshold via vsearch.
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md.)
derive_clustered <- function(input, output, identity = 0.99) {
  if (skip_if_exists(output, "clustered")) {
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

# Cutadapt-trim the primer-flanked region from every record (produces the
# `_cut` DB variants). Delegates to dbpq::cutadapt_rm_primers_db.
derive_cutadapted <- function(input, output, primer_fw, primer_rev, ...) {
  if (skip_if_exists(output, "cutadapted")) {
    return(normalizePath(output))
  }
  dbpq::cutadapt_rm_primers_db(
    ref_fasta = input,
    output = output,
    primer_fw = primer_fw,
    primer_rev = primer_rev,
    return_file_path = TRUE,
    ...
  )
}

# ---- fake reference set ----------------------------------------------------

# Build a balanced non-Fungi reference: one record per phylum present in
# `input`, then random fill to `n` records. Used by `add_external_seq_pq` in
# the pipeline to inject negative-control taxa into `d_asv`.
# (No dbpq equivalent — proposal `subset_balanced_db()` in
# proposals_for_dbpq.md.)
derive_fake_ref <- function(input  = "data/data_raw/refseq/dada2_format/Unite_wo_fungi.fasta",
                            output = "data/data_raw/fake_ref/fake_ref_asv_100.fasta",
                            n = 100,
                            seed = 22) {
  if (skip_if_exists(output, "fake_ref")) {
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
  topup     <- sample(remaining, max(0, n - length(one_per_phylum)))
  picks     <- unique(c(one_per_phylum, topup))

  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs[picks], output, width = 10000)
  normalizePath(output)
}

# ---- orchestration ---------------------------------------------------------

# Build every variant the benchmark consumes. Each step is idempotent;
# rerunning after a failure resumes from where it stopped.
derive_all_variants <- function() {
  base       <- here("data/data_raw/refseq")
  dada2_dir  <- file.path(base, "dada2_format")
  sintax_dir <- file.path(base, "sintax_format")
  dir.create(sintax_dir, showWarnings = FALSE, recursive = TRUE)

  # 1. EUKARYOME v1.9.3: strip parens to get the canonical sintax-format DBs.
  for (root in c("EUK_SSU", "EUK_ITS")) {
    src <- file.path(base, sprintf("SINTAX_%s_v1.9.3.fasta", root))
    if (file.exists(src)) {
      derive_no_parens(src,
                       file.path(sintax_dir, sprintf("%s_v1_9_3.fasta", root)))
    }
  }

  # 2. Fungi-only EUK variants (two-line records).
  for (db in c("EUK_SSU_v1_9_3", "EUK_ITS_v1_9_3")) {
    src <- file.path(sintax_dir, sprintf("%s.fasta", db))
    if (file.exists(src)) {
      derive_pattern_only(
        src,
        file.path(sintax_dir, sprintf("%s_Fungi.fasta", db)),
        pattern = "Fungi",
        multi_line_records = FALSE
      )
    }
  }

  # 3. UNITE: extract non-Fungi records for the fake-reference seed pool.
  unite_src <- file.path(dada2_dir, "Unite.fasta")
  unite_wo  <- file.path(dada2_dir, "Unite_wo_fungi.fasta")
  if (file.exists(unite_src)) {
    derive_no_pattern(unite_src, unite_wo, pattern = "Fungi")
  }

  # 4. Fake reference: 100 non-Fungi taxa, one per phylum + random fill.
  if (file.exists(unite_wo)) {
    derive_fake_ref(
      input  = unite_wo,
      output = here("data/data_raw/fake_ref/fake_ref_asv_100.fasta")
    )
  }

  # 5. Cutadapt-trimmed `_cut` variants. Primers come from config.R.
  source(here("config.R"))
  for (db in c("EUK_SSU_v1_9_3", "EUK_SSU_v1_9_3_Fungi")) {
    src <- file.path(sintax_dir, sprintf("%s.fasta", db))
    if (file.exists(src)) {
      derive_cutadapted(
        input      = src,
        output     = file.path(sintax_dir, sprintf("%s_cut.fasta", db)),
        primer_fw  = fw_primer_sequences,
        primer_rev = rev_primer_sequences
      )
    }
  }

  # 6. Mini variants for smoke testing (`mini_<original>.fasta`).
  for (fmt_dir in c(dada2_dir, sintax_dir)) {
    fastas <- list.files(fmt_dir, pattern = "\\.fasta$", full.names = TRUE)
    for (src in fastas) {
      bn <- basename(src)
      if (startsWith(bn, "mini_") || bn == "Unite_wo_fungi.fasta") {
        next
      }
      derive_mini(src, file.path(fmt_dir, paste0("mini_", bn)))
    }
  }

  invisible(NULL)
}
