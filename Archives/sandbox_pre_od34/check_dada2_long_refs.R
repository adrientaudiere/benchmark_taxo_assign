suppressMessages({
  library(targets)
  library(phyloseq)
})

Sys.setenv(TAR_PROJECT = "assign_taxo")
asv <- tar_read(d_asv_common)
seqs <- as.character(refseq(asv))

mini <- Biostrings::readDNAStringSet("data/data_raw/refseq/dada2_format/mini_EUK_ITS_v2.1_Fungi.fasta")
full_path <- "data/data_raw/refseq/dada2_format/EUK_ITS_v2.1.fasta"
idx <- Biostrings::fasta.index(full_path)
kingdom <- sub(";.*$", "", idx$desc)

long_rows <- which(idx$seqlength > 5000)
set.seed(22)
short_nonfungi_rows <- sample(which(idx$seqlength <= 5000 & kingdom != "Fungi"), length(long_rows))
long <- Biostrings::readDNAStringSet(idx[long_rows, ])
short_nonfungi <- Biostrings::readDNAStringSet(idx[sort(short_nonfungi_rows), ])

arms <- list(
  "A mini Fungi" = mini,
  "B mini Fungi + records > 5 kb" = c(mini, long),
  "C mini Fungi + non-Fungi records <= 5 kb" = c(mini, short_nonfungi)
)

for (arm in names(arms)) {
  path <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(arms[[arm]], path)
  t0 <- Sys.time()
  res <- dada2::assignTaxonomy(seqs, refFasta = path, minBoot = 50, multithread = 2)
  kingdoms <- sort(table(res[, 1], useNA = "ifany"), decreasing = TRUE)
  cat(sprintf(
    "DADA2 %-42s refs=%6d genus assigned=%5.1f %% kingdom assigned=%5.1f %% top kingdoms: %s (%.0f s)\n",
    arm, length(arms[[arm]]),
    100 * mean(!is.na(res[, "Genus"])), 100 * mean(!is.na(res[, 1])),
    paste(names(kingdoms)[1:min(3, length(kingdoms))], kingdoms[1:min(3, length(kingdoms))], collapse = ", "),
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))
}
