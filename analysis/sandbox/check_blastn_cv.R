suppressMessages({
  library(dplyr)
  library(phyloseq)
})
source("config.R")
source("R/load_pqverse.R")
suppressMessages(load_pqverse(c("MiscMetabar", "dbpq")))
source("R/cv_queries.R")

ref <- "data/data_raw/refseq/sintax_format/Unite_all_20250219.fasta"
n_query <- 200

dna <- Biostrings::readDNAStringSet(ref)
set.seed(22)
pool <- names(dna)[sample(length(dna), 2 * n_query)]
queries <- trim_cv_queries(
  dna[match(pool, names(dna))],
  primer_fw = fw_primer_sequences,
  primer_rev = rev_primer_sequences,
  min_overlap = primer_min_overlap,
  nproc = 2,
  prelude = cutadapt_conda_prelude
)
queries <- head(queries[Biostrings::width(queries) >= 50], n_query)
cat("queries:", length(queries), " median length:", median(Biostrings::width(queries)), "\n")

genus_of <- function(x) sub("^.*,g:([^,]*).*$", "\\1", x)
query_genus <- setNames(genus_of(names(queries)), names(queries))

tmp_dir <- tempfile("blastcv_")
dir.create(tmp_dir)
query_file <- file.path(tmp_dir, "queries.fasta")
db_file <- file.path(tmp_dir, "training.fasta")
out_file <- file.path(tmp_dir, "hits.tsv")

qnames <- names(queries)
names(queries) <- paste0("q", seq_along(queries))
Biostrings::writeXStringSet(queries, query_file)

for (variant in c("standard", "leaked")) {
  training <- if (variant == "standard") dna[!names(dna) %in% qnames] else dna
  Biostrings::writeXStringSet(training, db_file)
  system2("makeblastdb", c("-in", db_file, "-dbtype", "nucl", "-out", file.path(tmp_dir, "db")), stdout = FALSE)
  system2("blastn", c(
    "-query", query_file, "-db", file.path(tmp_dir, "db"),
    "-outfmt", shQuote("6 qseqid sseqid pident qcovs evalue bitscore"),
    "-max_target_seqs", "10", "-num_threads", "2", "-out", out_file
  ))
  hits <- read.delim(out_file, header = FALSE, col.names = c("q", "s", "pident", "qcovs", "evalue", "bitscore"))
  hits$query_genus <- query_genus[qnames[as.integer(sub("^q", "", hits$q))]]
  hits$hit_genus <- genus_of(hits$s)
  best <- hits |> group_by(q) |> slice_max(bitscore, n = 1, with_ties = FALSE) |> ungroup()
  kept <- hits |> filter(pident >= 90, qcovs >= blastn_min_cover, bitscore >= 50, evalue <= 1e-30)
  kept_best <- kept |> group_by(q) |> slice_max(bitscore, n = 1, with_ties = FALSE) |> ungroup()
  cat("\nVARIANT", variant, "\n")
  cat("  queries with any hit:", length(unique(hits$q)), "of", length(queries), "\n")
  cat("  best hit identity: median", round(median(best$pident), 1), " 10th pct", round(quantile(best$pident, 0.1), 1), " share >= 90:", round(100 * mean(best$pident >= 90), 1), "%\n")
  cat("  best hit query cover: median", round(median(best$qcovs), 1), " share >=", blastn_min_cover, ":", round(100 * mean(best$qcovs >= blastn_min_cover), 1), "%\n")
  cat("  queries passing every filter (id 90, cover", blastn_min_cover, ", bitscore 50, evalue 1e-30):", nrow(kept_best), "\n")
  cat("  best hit shares the query genus:", round(100 * mean(best$hit_genus == best$query_genus, na.rm = TRUE), 1), "%\n")
  cat("  same genus among filtered best hits:", round(100 * mean(kept_best$hit_genus == kept_best$query_genus, na.rm = TRUE), 1), "%\n")
  cat("  queries whose genus exists elsewhere in the reference:", round(100 * mean(query_genus %in% genus_of(names(training))), 1), "%\n")
}
unlink(tmp_dir, recursive = TRUE)
