# Does ITSx "full" (ITS1 + 5.8S + ITS2) still speed blastn up on a full-ITS
# amplicon (developer's question, 2026-09-27)? The in silico set gained 1.1x
# only, but its records carry almost no flank; the Tedersoo Illumina mock
# (ITS9mun / ITS4ngsUni) carries 154 bp of 18S and 39 bp of 28S, which "full"
# removes, while the 5.8S stays. This extracts ITS1, ITS2 and full from the 238
# Tedersoo ASVs and times blastn on each against Unite_all_20250219_Fungi, one
# computation at a time, through the pipeline's compute_assignment().
#   Rscript analysis/sandbox/timing_blastn_itsx_regions.R

here::i_am("analysis/sandbox/timing_blastn_itsx_regions.R")
source(here::here("config.R"))
source(here::here("R/load_pqverse.R"))
source(here::here("R/assign_compute.R"))
suppressMessages(load_pqverse(c("MiscMetabar", "comparpq")))
suppressMessages(library(phyloseq))

asv <- targets::tar_read(
  d_asv,
  store = here::here("store_dada2_tedersoo_illumina")
)
ref <- here::here(
  "data/data_raw/refseq/sintax_format/Unite_all_20250219_Fungi.fasta"
)
inputs <- list(raw = asv)
# itsx_pq() replaced the former R/itsx.R prototype on 2026-09-27; units
# identical after trimming are removed (keeps the most abundant), as the
# prototype did when these timings were measured.
for (region in c("full", "ITS1", "ITS2")) {
  inputs[[region]] <- MiscMetabar::itsx_pq(
    asv,
    region = region,
    organism_groups = itsx_organism_groups,
    duplicated_seqs = "remove",
    nproc = n_threads,
    args_before_itsx = itsx_conda_prelude,
    verbose = FALSE
  )
}
for (input in names(inputs)) {
  pq <- inputs[[input]]
  elapsed <- system.time(
    hits <- compute_assignment(pq, "blastn", ref, nproc = 1)
  )[["elapsed"]]
  n_hits <- table(hits[[1]])
  message(sprintf(
    "%-5s median length %4.0f bp | blastn %6.1f s (%.2f s per ASV) | hits per ASV: median %4.0f",
    input,
    median(Biostrings::width(pq@refseq)),
    elapsed,
    elapsed / ntaxa(pq),
    median(n_hits)
  ))
}
