suppressMessages({
  library(dplyr)
  library(phyloseq)
})
source("config.R")
source("R/load_pqverse.R")
suppressMessages(load_pqverse(c("MiscMetabar", "dbpq")))
source("R/cv_queries.R")
source("R/create_fake_pq_from_refseq.R")

ref <- "data/data_raw/refseq/sintax_format/Unite_all_20250219.fasta"
n_query <- 50

dna <- Biostrings::readDNAStringSet(ref)
set.seed(22)
pool <- names(dna)[sample(length(dna), 4 * n_query)]
queries <- trim_cv_queries(
  dna[match(pool, names(dna))],
  primer_fw = fw_primer_sequences,
  primer_rev = rev_primer_sequences,
  min_overlap = primer_min_overlap,
  nproc = 2,
  prelude = cutadapt_conda_prelude
)
queries <- head(queries[Biostrings::width(queries) >= 50], n_query)
cat("queries:", length(queries), " median query length:", median(Biostrings::width(queries)), "\n")

fake_pq <- create_fake_pq_from_refseq(queries)
training <- dna[!names(dna) %in% names(queries)]
train_file <- tempfile(fileext = ".fasta")
Biostrings::writeXStringSet(training, train_file)
cat("training records:", length(training), " median reference length:", median(Biostrings::width(training)), "\n")

tab <- blast_pq(
  physeq = fake_pq,
  fasta_for_db = train_file,
  unique_per_seq = FALSE,
  score_filter = FALSE,
  nproc = 2
)
cat("blast table rows:", nrow(tab), " columns:", paste(names(tab), collapse = " | "), "\n")
num <- function(x) suppressWarnings(as.numeric(x))
tab$cover <- num(tab$`Query cover`)
tab$id <- num(tab$`% id. match`)
tab$qlen <- num(tab$`Query seq. length`)
tab$slen <- num(tab$`Taxa seq. length`)
tab$alen <- num(tab$`Alignment length`)
cat("Query cover: median", median(tab$cover, na.rm = TRUE), " max", max(tab$cover, na.rm = TRUE), "\n")
cat("alignment / query length: median", round(median(tab$alen / tab$qlen, na.rm = TRUE), 3), "\n")
cat("alignment / subject length: median", round(median(tab$alen / tab$slen, na.rm = TRUE), 3), "\n")
best <- tab |> group_by(`Query name`) |> slice_max(num(`bit score`), n = 1, with_ties = FALSE) |> ungroup()
cat("queries with any hit:", nrow(best), "of", length(queries), "\n")
cat("best hit: median id", round(median(best$id), 1), " median cover", round(median(best$cover), 1), "\n")
kept <- tab |> filter(id >= 95, cover >= blastn_min_cover, num(`bit score`) >= 50, num(`e-value`) <= 1e-30)
cat("rows passing id 95 and cover", blastn_min_cover, ":", nrow(kept), " queries:", length(unique(kept$`Query name`)), "\n")
kept90 <- tab |> filter(id >= 90, cover >= blastn_min_cover)
cat("queries passing id 90 and cover", blastn_min_cover, ":", length(unique(kept90$`Query name`)), "\n")
unlink(train_file)
