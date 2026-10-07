suppressMessages({
  library(targets)
  library(phyloseq)
})

# Reconciles the 2026-09-24 ladder (check_dada2_length_vs_breadth.R) with the
# overflow bug of dada2 1.40.0 (check_dada2_alias.R): its arms held about
# 18 000 taxonomies, far below 65 536, and still lost the genus bootstrap once
# non-fungal records were added (arm C). Those arms were queried with the
# untrimmed ASVs (d_asv_common), which carry the conserved SSU end and 5.8S
# start that every eukaryote shares. Same arms (same seed), queried with the
# untrimmed and with the ITSx-trimmed real ASVs of Pauvert (ITS1, OD 34).
#   Rscript analysis/sandbox/check_dada2_ladder_itsx.R

queries <- list(
  untrimmed = as.character(refseq(tar_read(
    d_asv_common,
    store = "store_assign_taxo"
  ))),
  itsx_ITS1 = as.character(refseq(tar_read(
    d_ITS1_asv_real,
    store = "store_assign_taxo"
  )))
)

mini <- Biostrings::readDNAStringSet(
  "data/data_raw/refseq/dada2_format/mini_EUK_ITS_v2.1_Fungi.fasta"
)
full_path <- "data/data_raw/refseq/dada2_format/EUK_ITS_v2.1.fasta"
idx <- Biostrings::fasta.index(full_path)
kingdom <- sub(";.*$", "", idx$desc)
long_rows <- which(idx$seqlength > 5000)
set.seed(22)
short_nonfungi_rows <- sample(
  which(idx$seqlength <= 5000 & kingdom != "Fungi"),
  length(long_rows)
)
short_nonfungi <- Biostrings::readDNAStringSet(idx[sort(short_nonfungi_rows), ])

arms <- list(
  "A mini Fungi" = mini,
  "C + non-Fungi records <= 5 kb" = c(mini, short_nonfungi)
)

for (arm in names(arms)) {
  path <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(arms[[arm]], path)
  n_tax <- length(unique(names(arms[[arm]])))
  for (q in names(queries)) {
    res <- dada2::assignTaxonomy(
      queries[[q]],
      refFasta = path,
      minBoot = 0,
      outputBootstraps = TRUE,
      multithread = 3
    )
    boot <- res$boot
    cat(sprintf(
      "%-32s %-10s taxonomies %6d | genus boot >= 50: %5.1f %% | median boot Kingdom %3.0f Genus %3.0f | kingdom Fungi %5.1f %%\n",
      arm,
      q,
      n_tax,
      100 * mean(boot[, "Genus"] >= 50),
      median(boot[, "Kingdom"]),
      median(boot[, "Genus"]),
      100 * mean(res$tax[, "Kingdom"] == "Fungi", na.rm = TRUE)
    ))
  }
  unlink(path)
}
