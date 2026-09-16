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
genus_of <- function(x) sub("^.*,g:([^,]*).*$", "\\1", x)
qnames <- names(queries)
query_genus <- setNames(genus_of(qnames), paste0("q", seq_along(qnames)))

tmp_dir <- tempfile("blastvote_")
dir.create(tmp_dir)
query_file <- file.path(tmp_dir, "queries.fasta")
db_file <- file.path(tmp_dir, "training.fasta")
out_file <- file.path(tmp_dir, "hits.tsv")
names(queries) <- paste0("q", seq_along(queries))
Biostrings::writeXStringSet(queries, query_file)

training <- dna[!names(dna) %in% qnames]
Biostrings::writeXStringSet(training, db_file)
system2("makeblastdb", c("-in", db_file, "-dbtype", "nucl", "-out", file.path(tmp_dir, "db")), stdout = FALSE)
system2("blastn", c(
  "-query", query_file, "-db", file.path(tmp_dir, "db"),
  "-outfmt", shQuote("6 qseqid sseqid pident qcovs evalue bitscore"),
  "-max_target_seqs", "500", "-num_threads", "2", "-out", out_file
))
hits <- read.delim(out_file, header = FALSE, col.names = c("q", "s", "pident", "qcovs", "evalue", "bitscore"))
hits$hit_genus <- genus_of(hits$s)
hits$query_genus <- query_genus[hits$q]

vote_share <- function(h, id_cut, nb_voting = 100) {
  kept <- h |>
    filter(pident >= id_cut, qcovs >= blastn_min_cover, bitscore >= 50, evalue <= 1e-30) |>
    arrange(q, desc(bitscore)) |>
    group_by(q) |>
    slice_head(n = nb_voting)
  voted <- kept |>
    summarise(
      query_genus = first(query_genus),
      n_hits = n(),
      top = {
        tb <- sort(table(hit_genus), decreasing = TRUE)
        if (length(tb) > 1 && tb[1] == tb[2]) NA_character_ else names(tb)[1]
      },
      .groups = "drop"
    )
  cat(sprintf(
    "  id >= %2d : %3d queries with hits, %3d assigned by rel_majority (%4.1f %% of %d), good genus %4.1f %%, median hits kept %s\n",
    id_cut, nrow(voted), sum(!is.na(voted$top)), 100 * sum(!is.na(voted$top)) / n_query, n_query,
    100 * mean(voted$top == voted$query_genus, na.rm = TRUE), median(voted$n_hits)
  ))
}

cat("Standard variant, 200 trimmed queries, cover >=", blastn_min_cover, "\n")
for (id_cut in c(95, 90, 85, 80)) {
  vote_share(hits, id_cut)
}
unlink(tmp_dir, recursive = TRUE)
