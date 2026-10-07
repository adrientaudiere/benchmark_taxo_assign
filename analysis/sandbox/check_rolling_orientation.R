# Read orientation of the Rolling et al. 2022 mock (ROADMAP §1, third mock;
# SRA PRJNA779405). The authors write that forward and reverse reads are mixed
# in R1 and R2 (their library preparation), and they split them by primer
# before DADA2. MiscMetabar::cutadapt_remove_primers(), as
# pipelines/dada2_bio.R calls it, looks for the forward primer in R1 and the
# reverse primer in R2 with --discard-untrimmed, so it may keep one
# orientation only. This script measures it on one run:
#
#   1. on the first 50 000 pairs, which primer starts R1 and R2;
#   2. cutadapt_remove_primers() as the pipeline calls it (fw, rev);
#   3. the same call with the primers swapped (rev, fw), which should keep the
#      reversed pairs, and whether the two passes share any read;
#   4. one call with `cutadapt_args = "--revcomp"` (cutadapt >= 4.6), checked
#      against the two passes.
#
# Reads data/data_raw/rawseq_rolling_mock/ (<run>_R1 / _R2.fastq.gz, from
# ENA); writes data/data_intermediate/rolling_orientation_test/. Run from the
# project root:
#   Rscript analysis/sandbox/check_rolling_orientation.R SRR16904671

suppressMessages({
  library(here)
  library(ShortRead)
  source(here("config.R"))
  source(here("R/load_pqverse.R"))
  load_pqverse("MiscMetabar")
})

run <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(run)) {
  run <- "SRR16904671"
}
raw_dir <- here("data/data_raw/rawseq_rolling_mock")
out_dir <- here("data/data_intermediate/rolling_orientation_test", run)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ITS-1-F and 5.8S-1R (Rolling et al. 2022, Methods); the wobble variant of
# ITS-1-F (S at position 19) is matched by the same pattern.
fw <- "CTTGGTCATTTAGAGGAASTAA"
rev <- "GTTCAAAGAYTCGATGATTCAC"

# ---- 1. which primer starts each read ----------------------------------------
first_reads <- function(path, n = 50000) {
  stream <- FastqStreamer(path, n = n)
  on.exit(close(stream))
  sread(yield(stream))
}
r1 <- first_reads(file.path(raw_dir, paste0(run, "_R1.fastq.gz")))
r2 <- first_reads(file.path(raw_dir, paste0(run, "_R2.fastq.gz")))
starts_with <- function(reads, primer, window = 60) {
  head_seq <- subseq(reads, 1, pmin(width(reads), window))
  vcountPattern(primer, head_seq, max.mismatch = 2, fixed = FALSE) > 0
}
classify <- function(reads) {
  is_fw <- starts_with(reads, fw)
  is_rev <- starts_with(reads, rev)
  ifelse(
    is_fw & !is_rev,
    "fw",
    ifelse(is_rev & !is_fw, "rev", ifelse(is_fw, "both", "none"))
  )
}
cat(
  "\nPrimer within the first 60 bp, first",
  length(r1),
  "pairs (rows R1, columns R2)\n"
)
print(table(R1 = classify(r1), R2 = classify(r2)))

# ---- 2 and 3. the pipeline's cutadapt call, then with the primers swapped -----
run_cutadapt <- function(label, primer_fw, primer_rev, cutadapt_args = "") {
  dir <- file.path(out_dir, label)
  unlink(dir, recursive = TRUE)
  cutadapt_remove_primers(
    path_to_fastq = raw_dir,
    pattern = run,
    primer_fw = primer_fw,
    primer_rev = primer_rev,
    folder_output = dir,
    nproc = n_threads,
    cutadapt_args = cutadapt_args,
    args_before_cutadapt = cutadapt_conda_prelude,
    verbose = FALSE
  )
  json <- jsonlite::read_json(Sys.glob(file.path(dir, "*.cutadapt.json"))[1])
  r1_out <- readFastq(Sys.glob(file.path(dir, "*_R1.fastq.gz"))[1])
  r2_out <- readFastq(Sys.glob(file.path(dir, "*_R2.fastq.gz"))[1])
  names_out <- as.character(id(r1_out))
  list(
    label = label,
    pairs_in = json$read_counts$input,
    pairs_out = json$read_counts$output,
    ids = sub(" .*$", "", names_out),
    rc = grepl(" rc$", names_out),
    r1 = as.character(sread(r1_out)),
    r2 = as.character(sread(r2_out))
  )
}
standard <- run_cutadapt("standard", fw, rev)
swapped <- run_cutadapt("swapped", rev, fw)
# 4. Cutadapt >= 4.6: --revcomp also tries each pair with R1 and R2 swapped and
# keeps the orientation that matches (suffix " rc" on the swapped pairs).
revcomp <- run_cutadapt("revcomp", fw, rev, cutadapt_args = "--revcomp")

cat("\ncutadapt_remove_primers(), pairs in and written\n")
for (x in list(standard, swapped, revcomp)) {
  cat(sprintf(
    "  %-9s %9d -> %9d (%.1f %%)\n",
    x$label,
    x$pairs_in,
    x$pairs_out,
    100 * x$pairs_out / x$pairs_in
  ))
}
shared <- length(intersect(standard$ids, swapped$ids))
both <- length(union(standard$ids, swapped$ids))
cat(sprintf(
  "  pairs kept by both passes: %d; by either: %d (%.1f %% of the input)\n",
  shared,
  both,
  100 * both / standard$pairs_in
))

# --revcomp against the two passes: same pairs kept, and the same trimmed
# reads. An unswapped pair must equal the standard pass. A swapped (rc) pair
# has its reads exchanged: its R1 is the R2 of the swapped pass (the read that
# starts with ITS-1-F) and its R2 the R1 of that pass, so every R1 of the
# --revcomp output starts after ITS-1-F.
cat("\n--revcomp against the two passes\n")
cat(sprintf(
  "  pairs written: %d, of which %d swapped (rc); in either pass: %d; only in --revcomp: %d; only in the passes: %d\n",
  length(revcomp$ids),
  sum(revcomp$rc),
  sum(revcomp$ids %in% union(standard$ids, swapped$ids)),
  length(setdiff(revcomp$ids, union(standard$ids, swapped$ids))),
  length(setdiff(union(standard$ids, swapped$ids), revcomp$ids))
))
same_reads <- function(pass, keep, exchanged) {
  ids <- intersect(revcomp$ids[keep], pass$ids)
  i <- match(ids, revcomp$ids)
  j <- match(ids, pass$ids)
  if (exchanged) {
    mean(revcomp$r1[i] == pass$r2[j] & revcomp$r2[i] == pass$r1[j])
  } else {
    mean(revcomp$r1[i] == pass$r1[j] & revcomp$r2[i] == pass$r2[j])
  }
}
cat(sprintf(
  "  identical pairs: %.4f of the unswapped pairs vs the standard pass, %.4f of the rc pairs vs the swapped pass (reads exchanged)\n",
  same_reads(standard, !revcomp$rc, exchanged = FALSE),
  same_reads(swapped, revcomp$rc, exchanged = TRUE)
))
