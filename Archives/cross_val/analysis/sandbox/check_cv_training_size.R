suppressMessages({
  library(dplyr)
  library(phyloseq)
})
source("config.R")
source("R/load_pqverse.R")
suppressMessages(load_pqverse(c("MiscMetabar", "comparpq", "tidypq", "dbpq")))
source("R/cross_val.R")
source("R/cv_queries.R")
source("R/create_fake_pq_from_refseq.R")

ref <- "data/data_raw/refseq/sintax_format/Unite_all_20250219.fasta"
cat("records in the reference file:", length(Biostrings::fasta.index(ref)$desc), "\n")

trace_write <- Biostrings::writeXStringSet
assign("writeXStringSet", function(x, filepath, ...) {
  if (grepl("cv_refseq_", filepath)) {
    cat("TRAIN fasta written:", length(x), "records\n")
  }
  trace_write(x, filepath, ...)
}, envir = globalenv())

trace_fake <- create_fake_pq_from_refseq
create_fake_pq_from_refseq <- function(...) {
  out <- trace_fake(...)
  cat("TRAIN queries in the fold:", ntaxa(out), "\n")
  out
}

res <- cross_val(
  ref,
  fold_number = 10,
  fold_tested = 1,
  method = "blastn",
  seed = targets_seed,
  max_seq = 500,
  nproc = 2,
  remove_tested_sequences = TRUE,
  primer_fw = fw_primer_sequences,
  primer_rev = rev_primer_sequences,
  primer_min_overlap = primer_min_overlap,
  oversample = cv_oversample,
  cutadapt_prelude = cutadapt_conda_prelude,
  vote_algorithm = "rel_majority",
  nb_voting = 100L,
  min_cover = blastn_min_cover
)
g <- res$metrics |> filter(name == "Genus")
cat("TRAIN Genus:", paste(g$metric, round(g$mean, 3), sep = "=", collapse = " "), "\n")
