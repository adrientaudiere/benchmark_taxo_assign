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

trace_blastn <- assign_blastn
assign_blastn <- function(...) {
  args <- list(...)
  cat("REPRO assign_blastn args:", paste(names(args), collapse = ", "), "\n")
  cat(
    "REPRO min_id:", if (is.null(args$min_id)) "default 95" else args$min_id,
    " min_cover:", if (is.null(args$min_cover)) "default 95" else args$min_cover,
    " vote:", if (is.null(args$vote_algorithm)) "default" else args$vote_algorithm,
    " nb_voting:", if (is.null(args$nb_voting)) "NULL" else args$nb_voting, "\n"
  )
  out <- do.call(trace_blastn, args)
  if (inherits(out, "phyloseq")) {
    tt <- as.data.frame(as(tax_table(out), "matrix"))
    cols <- grep("_blastn$", names(tt), value = TRUE)
    if (length(cols) > 0) {
      cat("REPRO blastn columns:", length(cols), " genus assigned:",
          round(100 * mean(!is.na(tt[[grep("^Genus", cols, value = TRUE)[1]]])), 1), "%\n")
    } else {
      cat("REPRO blastn added no column\n")
    }
  }
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
cat("REPRO Genus:", paste(g$metric, round(g$mean, 3), sep = "=", collapse = " "), "\n")
