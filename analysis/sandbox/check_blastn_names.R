suppressMessages({
  library(dplyr)
  library(phyloseq)
})
source("config.R")
source("R/load_pqverse.R")
suppressMessages(load_pqverse(c("MiscMetabar", "tidypq", "dbpq")))
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

fake_pq <- create_fake_pq_from_refseq(queries)
training <- dna[!names(dna) %in% names(queries)]
train_file <- tempfile(fileext = ".fasta")
Biostrings::writeXStringSet(training, train_file)

cat("names of queries (first):", head(names(queries), 1), "\n")
cat("taxa_names(fake_pq) (first):", head(taxa_names(fake_pq), 1), "\n")
cat("identical query names and taxa_names:", identical(unname(names(queries)), taxa_names(fake_pq)), "\n")

voted <- assign_blastn(
  fake_pq,
  ref_fasta = train_file,
  behavior = "return_matrix",
  nproc = 2,
  vote_algorithm = "rel_majority",
  nb_voting = 100L,
  min_cover = blastn_min_cover
)

if (is.list(voted) && !is.data.frame(voted)) {
  cat("return_matrix gave a list:", paste(names(voted), collapse = " | "), "\n")
  cat("raw_blast_table rows:", nrow(voted$raw_blast_table), "\n")
  voted <- voted$blast_table_per_query
}

if (is.null(voted)) {
  cat("assign_blastn returned NULL\n")
} else {
  cat("rows per query: median", median(table(voted$taxa_names)), " max", max(table(voted$taxa_names)), "\n")
  cat("voted rows:", nrow(voted), " columns:", paste(head(names(voted), 4), collapse = " | "), "\n")
  cat("voted taxa_names (first):", head(voted$taxa_names, 1), "\n")
  cat("voted names found in taxa_names(fake_pq):", sum(voted$taxa_names %in% taxa_names(fake_pq)), "of", nrow(voted), "\n")
  gcol <- grep("^Genus", names(voted), value = TRUE)[1]
  cat("genus assigned among voted rows:", sum(!is.na(voted[[gcol]])), "\n")
}

added <- assign_blastn(
  fake_pq,
  ref_fasta = train_file,
  behavior = "add_to_phyloseq",
  nproc = 2,
  vote_algorithm = "rel_majority",
  nb_voting = 100L,
  min_cover = blastn_min_cover
)
tt <- as.data.frame(as(tax_table(added), "matrix"))
gcol <- grep("^Genus.*_blastn$", names(tt), value = TRUE)[1]
cat("after add_to_phyloseq, genus assigned:", sum(!is.na(tt[[gcol]])), "of", ntaxa(added), "\n")
unlink(train_file)
