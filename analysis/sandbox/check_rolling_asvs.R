# ASVs of the Rolling et al. 2022 mock (OD 42) against the 12 expected ITS1
# amplicons of the authors, before any assignment: reads kept at each step,
# which expected sequences come out as an ASV (exact match), and the truth
# that R/mock_truth.R gives each ASV. Reads store_dada2_rolling_mock only.
#   Rscript analysis/sandbox/check_rolling_asvs.R

suppressMessages({
  library(here)
  library(targets)
  library(phyloseq)
  library(dplyr)
  library(Biostrings)
  source(here("config.R"))
  source(here("R/load_pqverse.R"))
  load_pqverse("dbpq")
  source(here("R/mock_truth.R"))
})
store <- here("store_dada2_rolling_mock")
cfg <- mock_sanger_tables[mock_sanger_tables$dataset == "rolling_mock", ]

cat("\nReads per step\n")
print(tar_read(track_df, store = store))

d_asv <- tar_read(d_asv, store = store)
asv_seq <- stats::setNames(as.character(d_asv@refseq), taxa_names(d_asv))
asv_reads <- taxa_sums(d_asv)
expected <- read_mock_sanger(cfg)
cat(
  "\n",
  ntaxa(d_asv),
  " ASVs, ",
  sum(asv_reads),
  " reads, ",
  nsamples(d_asv),
  " samples\n",
  sep = ""
)

# The pipeline trims one base at the start of R1 and of R2 (filter_trim(),
# trimLeft = 1), so an ASV is its amplicon less its first and last base.
expected_asv <- substr(expected$sequence, 2, nchar(expected$sequence) - 1)
cat("\nExpected amplicons recovered as an ASV (identical, ends trimmed)\n")
print(
  tibble(
    strain = expected$strain,
    length = nchar(expected_asv),
    asv = names(asv_seq)[match(expected_asv, asv_seq)],
    reads = unname(asv_reads[match(expected_asv, asv_seq)]),
    samples = vapply(
      match(expected_asv, asv_seq),
      \(k) if (is.na(k)) 0L else sum(otu_table(d_asv)[, k] > 0),
      integer(1)
    )
  ),
  n = Inf
)

truth <- build_mock_truth(d_asv, cfg)
cat("\nTruth of the ASVs (R/mock_truth.R)\n")
print(table(depth = truth$truth_depth, useNA = "ifany"))
cat("foreign:", sum(truth$foreign), "\n")
cat(
  "reads in ASVs identical to an expected amplicon:",
  round(100 * sum(asv_reads[asv_seq %in% expected_asv]) / sum(asv_reads), 1),
  "%\n"
)
truth_class <- ifelse(
  truth$foreign,
  "foreign",
  ifelse(is.na(truth$truth_depth), "no truth", truth$truth_depth)
)
cat("\nASVs and reads per truth class\n")
print(
  tibble(class = truth_class, n_reads = asv_reads[truth$taxon]) |>
    group_by(class) |>
    summarise(
      asv = n(),
      reads = sum(n_reads),
      share = round(100 * sum(n_reads) / sum(asv_reads), 2),
      max_reads = max(n_reads)
    )
)
