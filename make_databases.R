# Reproducible derivation of every reference-database variant the benchmark
# consumes. Replaces the historical `some_bash_script` (now in Archives/).
#
# Most operations are delegated to `dbpq::` — only the steps for which dbpq
# does not (yet) have an equivalent stay local. Those are flagged in
# proposals_for_dbpq.md and should migrate to dbpq when added.
#
# Every function is idempotent: if the output file already exists the function
# returns its path without rerunning (unless force = TRUE is passed, which
# deletes the existing file and rebuilds it).
#
# Sourcing this file has no side effects. Build everything with:
#   source("make_databases.R")
#   derive_all_variants()            # idempotent, resumes from last failure
#   derive_all_variants(force = TRUE) # rebuild all from scratch
#
# Three source files in data/data_raw/refseq/ drive all derivations:
#   Unite.fasta      — UNITE dada2/k__ format (input for dada2 and sintax)
#   Euk_ITS_v2.fasta — EUKARYOME ITS v2, k__ format
#   Euk_SSU_v2.fasta — EUKARYOME SSU v2, k__ format
#
# External tools required: bash, sed, grep, head, vsearch, plus cutadapt in
# the `cutadaptenv` conda env (see config.R::cutadapt_conda_prelude).

library("here")
library("dbpq")
here::i_am("make_databases.R")

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

# ---- format conversion (delegated to dbpq) ---------------------------------

# Plain file copy — for source files already in the right format (e.g. the
# UNITE dada2/k__ file copied verbatim into dada2_format/).
derive_copy <- function(input, output, force = FALSE) {
  if (skip_if_exists(output, "copy", force = force)) {
    return(normalizePath(output))
  }
  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  file.copy(input, output, overwrite = TRUE)
  normalizePath(output)
}

# dbpq::format2sintax() auto-detects the input format (UNITE/dada2 k__ style,
# Greengenes2, plain) and writes the sintax-tagged fasta (tax= style).
derive_sintax <- function(input, output, input_format = "auto", force = FALSE) {
  if (skip_if_exists(output, "sintax", force = force)) {
    return(normalizePath(output))
  }
  dbpq::format2sintax(
    fasta_db = input,
    input_format = input_format,
    output_path = output
  )
  normalizePath(output)
}

# dbpq::format2dada2() auto-detects the input format and writes the dada2
# pure-rank format (">Kingdom;Phylum;..." headers) expected by assignTaxonomy.
derive_dada2 <- function(input, output, input_format = "auto", force = FALSE) {
  if (skip_if_exists(output, "dada2", force = force)) {
    return(normalizePath(output))
  }
  dbpq::format2dada2(
    fasta_db = input,
    input_format = input_format,
    output_path = output
  )
  normalizePath(output)
}

# Strip parenthesized synonyms found in some EUKARYOME headers.
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md.)
derive_no_parens <- function(input, output, force = FALSE) {
  if (skip_if_exists(output, "no_parens", force = force)) {
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
                                multi_line_records = TRUE,
                                force = FALSE) {
  if (skip_if_exists(output, paste0(pattern, "_only"), force = force)) {
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
# yet — proposal in proposals_for_dbpq.md.
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

# First `n` lines as a quick smoke-test subset.
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
# (No dbpq equivalent yet — proposal in proposals_for_dbpq.md.)
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

# Cutadapt-trim the primer-flanked region from every record (produces the
# `_cut` DB variants). Delegates to dbpq::cutadapt_rm_primers_db.
derive_cutadapted <- function(input, output, primer_fw, primer_rev,
                              force = FALSE, ...) {
  if (skip_if_exists(output, "cutadapted", force = force)) {
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
                            seed = 22,
                            force = FALSE) {
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
  topup     <- sample(remaining, max(0, n - length(one_per_phylum)))
  picks     <- unique(c(one_per_phylum, topup))

  dir.create(dirname(output), showWarnings = FALSE, recursive = TRUE)
  Biostrings::writeXStringSet(seqs[picks], output, width = 10000)
  normalizePath(output)
}

# ---- orchestration ---------------------------------------------------------

# Build every variant the benchmark consumes from the three source files:
#   data/data_raw/refseq/Unite.fasta
#   data/data_raw/refseq/Euk_ITS_v2.fasta
#   data/data_raw/refseq/Euk_SSU_v2.fasta
#
# Each step is idempotent; rerunning after a failure resumes from where it
# stopped. Pass force = TRUE to delete and rebuild all derived files.
derive_all_variants <- function(force = FALSE) {
  source(here("config.R"))

  base       <- here("data/data_raw/refseq")
  dada2_dir  <- file.path(base, "dada2_format")
  sintax_dir <- file.path(base, "sintax_format")
  dir.create(dada2_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(sintax_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- 1. UNITE ----------------------------------------------------------
  unite_src <- file.path(base, "Unite.fasta")

  if (file.exists(unite_src)) {
    # 1a. Copy to dada2_format/ (UNITE k__ format, assignTaxonomy reads it
    #     directly without conversion).
    unite_dada2 <- derive_copy(
      unite_src,
      file.path(dada2_dir, "Unite.fasta"),
      force = force
    )

    # 1b. Non-Fungi records: seed pool for the fake-reference set.
    unite_wo <- derive_no_pattern(
      unite_dada2,
      file.path(dada2_dir, "Unite_wo_fungi.fasta"),
      pattern = "Fungi",
      force = force
    )

    # 1c. Convert to sintax format (k__ → tax= style) for vsearch/LCA/blastn.
    unite_sintax <- derive_sintax(
      unite_dada2,
      file.path(sintax_dir, "Unite.fasta"),
      force = force
    )

    # 1d. Fungi-only variants in both formats.
    derive_pattern_only(
      unite_dada2,
      file.path(dada2_dir, "Unite_Fungi.fasta"),
      pattern = "Fungi",
      force = force
    )
    derive_pattern_only(
      unite_sintax,
      file.path(sintax_dir, "Unite_Fungi.fasta"),
      pattern = "Fungi",
      force = force
    )

    # 1e. Fake reference: 100 non-Fungi taxa, one per phylum + random fill.
    derive_fake_ref(
      input  = unite_wo,
      output = here("data/data_raw/fake_ref/fake_ref_asv_100.fasta"),
      force  = force
    )
  } else {
    warning("Unite.fasta not found in ", base, " — skipping UNITE derivations.")
  }

  # ---- 2. EUKARYOME ITS v2 -----------------------------------------------
  euk_its_src <- file.path(base, "Euk_ITS_v2.fasta")

  if (file.exists(euk_its_src)) {
    # 2a. Convert to dada2 pure-rank format (">Kingdom;Phylum;..." headers).
    euk_its_dada2 <- derive_dada2(
      euk_its_src,
      file.path(dada2_dir, "EUK_ITS_v2.fasta"),
      force = force
    )

    # 2b. Convert to sintax format (tax= style).
    euk_its_sintax <- derive_sintax(
      euk_its_src,
      file.path(sintax_dir, "EUK_ITS_v2.fasta"),
      force = force
    )

    # 2c. Fungi-only variants.
    euk_its_dada2_fungi <- derive_pattern_only(
      euk_its_dada2,
      file.path(dada2_dir, "EUK_ITS_v2_Fungi.fasta"),
      pattern = "Fungi",
      force = force
    )
    euk_its_sintax_fungi <- derive_pattern_only(
      euk_its_sintax,
      file.path(sintax_dir, "EUK_ITS_v2_Fungi.fasta"),
      pattern = "Fungi",
      force = force
    )

    # 2d. Cutadapt-trimmed Fungi variants (ITS primers from config.R).
    derive_cutadapted(
      input      = euk_its_dada2_fungi,
      output     = file.path(dada2_dir, "EUK_ITS_v2_Fungi_cut.fasta"),
      primer_fw  = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      force      = force
    )
    derive_cutadapted(
      input      = euk_its_sintax_fungi,
      output     = file.path(sintax_dir, "EUK_ITS_v2_Fungi_cut.fasta"),
      primer_fw  = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      force      = force
    )
  } else {
    warning("Euk_ITS_v2.fasta not found in ", base, " — skipping EUK ITS derivations.")
  }

  # ---- 3. EUKARYOME SSU v2 -----------------------------------------------
  euk_ssu_src <- file.path(base, "Euk_SSU_v2.fasta")

  if (file.exists(euk_ssu_src)) {
    # 3a. Convert to dada2 pure-rank format.
    euk_ssu_dada2 <- derive_dada2(
      euk_ssu_src,
      file.path(dada2_dir, "EUK_SSU_v2.fasta"),
      force = force
    )

    # 3b. Convert to sintax format.
    euk_ssu_sintax <- derive_sintax(
      euk_ssu_src,
      file.path(sintax_dir, "EUK_SSU_v2.fasta"),
      force = force
    )

    # 3c. Fungi-only variants.
    euk_ssu_dada2_fungi <- derive_pattern_only(
      euk_ssu_dada2,
      file.path(dada2_dir, "EUK_SSU_v2_Fungi.fasta"),
      pattern = "Fungi",
      force = force
    )
    euk_ssu_sintax_fungi <- derive_pattern_only(
      euk_ssu_sintax,
      file.path(sintax_dir, "EUK_SSU_v2_Fungi.fasta"),
      pattern = "Fungi",
      force = force
    )

    # 3d. Cutadapt-trimmed Fungi variants (SSU primers from config.R).
    derive_cutadapted(
      input      = euk_ssu_dada2_fungi,
      output     = file.path(dada2_dir, "EUK_SSU_v2_Fungi_cut.fasta"),
      primer_fw  = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      force      = force
    )
    derive_cutadapted(
      input      = euk_ssu_sintax_fungi,
      output     = file.path(sintax_dir, "EUK_SSU_v2_Fungi_cut.fasta"),
      primer_fw  = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      force      = force
    )
  } else {
    warning("Euk_SSU_v2.fasta not found in ", base, " — skipping EUK SSU derivations.")
  }

  # ---- 4. Mini variants for smoke testing --------------------------------
  # Produces mini_<name>.fasta for every derived file (first 10 000 lines).
  # Skips Unite_wo_fungi.fasta (internal intermediate, not used by the
  # pipeline directly) and any file already prefixed with mini_.
  for (fmt_dir in c(dada2_dir, sintax_dir)) {
    fastas <- list.files(fmt_dir, pattern = "\\.fasta$", full.names = TRUE)
    for (src in fastas) {
      bn <- basename(src)
      if (startsWith(bn, "mini_") || bn == "Unite_wo_fungi.fasta") {
        next
      }
      derive_mini(src, file.path(fmt_dir, paste0("mini_", bn)), force = force)
    }
  }

  invisible(NULL)
}
