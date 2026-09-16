suppressMessages({
  library(dplyr)
  library(phyloseq)
})
source("config.R")
source("R/load_pqverse.R")
suppressMessages(load_pqverse(c("MiscMetabar", "tidypq")))
source("R/create_fake_pq_from_refseq.R")

ref <- "data/data_raw/refseq/sintax_format/mini_Unite_all_20250219.fasta"
dna <- Biostrings::readDNAStringSet(ref)
set.seed(22)
q <- dna[sample(length(dna), 2000)]
q <- q[!duplicated(as.character(q))]
fake_pq <- create_fake_pq_from_refseq(q)
truth <- tidypq::tax_table_to_df(fake_pq, convert = FALSE)

for (np in c(1, 4)) {
  res <- assign_sintax(fake_pq, ref_fasta = ref, nproc = np, behavior = "return_matrix", min_bootstrap = 0.5)
  tv <- res$taxo_value
  same_order <- identical(tv$taxa_names, taxa_names(fake_pq))
  n_match <- sum(tv$taxa_names %in% taxa_names(fake_pq))
  first_diff <- which(tv$taxa_names != taxa_names(fake_pq))[1]
  positional <- mean(tv$Genus == truth$Genus, na.rm = FALSE)
  joined <- tibble(taxa_names = taxa_names(fake_pq), truth = truth$Genus) |>
    left_join(select(tv, taxa_names, Genus), by = "taxa_names")
  cat(sprintf(
    "SINTAX nproc=%d rows=%d names matched=%d same order=%s first diff=%s | genus good positional=%.3f joined=%.3f\n",
    np, nrow(tv), n_match, same_order, first_diff,
    sum(tv$Genus == truth$Genus, na.rm = TRUE) / nrow(truth),
    sum(joined$Genus == joined$truth, na.rm = TRUE) / nrow(truth)
  ))
}
