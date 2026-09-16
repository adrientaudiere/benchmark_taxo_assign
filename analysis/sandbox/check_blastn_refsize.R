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
pool_size <- 5413

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
genus_of <- function(x) sub("^.*,g:([^,]*).*$", "\\1", x)
qnames <- names(queries)
query_genus <- setNames(genus_of(qnames), paste0("q", seq_along(qnames)))
names(queries) <- paste0("q", seq_along(queries))

tmp_dir <- tempfile("refsize_")
dir.create(tmp_dir)
query_file <- file.path(tmp_dir, "queries.fasta")
db_file <- file.path(tmp_dir, "training.fasta")
out_file <- file.path(tmp_dir, "hits.tsv")
Biostrings::writeXStringSet(queries, query_file)

others <- dna[!names(dna) %in% qnames]

arms <- list(
  "full database" = others,
  "CV pool size (5413 records)" = others[sample(length(others), pool_size - n_query)]
)

for (arm in names(arms)) {
  training <- arms[[arm]]
  Biostrings::writeXStringSet(training, db_file)
  system2("makeblastdb", c("-in", db_file, "-dbtype", "nucl", "-out", file.path(tmp_dir, "db")), stdout = FALSE)
  t0 <- Sys.time()
  system2("blastn", c(
    "-query", query_file, "-db", file.path(tmp_dir, "db"),
    "-outfmt", shQuote("6 qseqid sseqid pident qcovs evalue bitscore"),
    "-max_target_seqs", "500", "-num_threads", "2", "-out", out_file
  ))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  hits <- read.delim(out_file, header = FALSE, col.names = c("q", "s", "pident", "qcovs", "evalue", "bitscore"))
  hits$hit_genus <- genus_of(hits$s)
  hits$query_genus <- query_genus[hits$q]
  cat("\nARM", arm, " references:", length(training), sprintf(" blast %.1f s\n", secs))
  for (id_cut in c(95, 90)) {
    voted <- hits |>
      filter(pident >= id_cut, qcovs >= blastn_min_cover, bitscore >= 50, evalue <= 1e-30) |>
      arrange(q, desc(bitscore)) |>
      group_by(q) |>
      slice_head(n = 100) |>
      summarise(
        query_genus = first(query_genus),
        top = {
          tb <- sort(table(hit_genus), decreasing = TRUE)
          if (length(tb) > 1 && tb[1] == tb[2]) NA_character_ else names(tb)[1]
        },
        .groups = "drop"
      )
    cat(sprintf(
      "  id >= %d : %5.1f %% of the 200 queries get a genus, good genus %5.1f %%\n",
      id_cut, 100 * sum(!is.na(voted$top)) / n_query,
      100 * mean(voted$top == voted$query_genus, na.rm = TRUE)
    ))
  }
}
unlink(tmp_dir, recursive = TRUE)
