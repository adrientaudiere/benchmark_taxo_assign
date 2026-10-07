# Why does dada2 abstain on the full EUK_ITS_v2.1 (ROADMAP C22-P5)? The
# 2026-09-16 test (check_dada2_long_refs.R) added either the records over 5 kb
# (arm B) or as many non-Fungi records <= 5 kb (arm C) to mini_EUK_ITS_v2.1_Fungi,
# but left length and taxonomic breadth confounded. Arm D separates them: the
# same records as B, each cut to one random window with a length drawn from the
# records of C, so B and D hold the same taxa and differ only in length.
#
#   B vs D: effect of length at constant taxa
#   D vs C: effect of which taxa are added, at typical length
#
# Run from the project root (a few minutes; dada2 on ~18 000 references):
#   Rscript analysis/sandbox/check_dada2_length_vs_breadth.R

suppressMessages({
  library(targets)
  library(phyloseq)
})

queries <- list(
  pauvert = as.character(refseq(tar_read(
    d_asv_common,
    store = "store_assign_taxo"
  ))),
  insilico = as.character(refseq(tar_read(
    d_asv_common,
    store = "store_assign_taxo_insilico"
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
long <- Biostrings::readDNAStringSet(idx[long_rows, ])
short_nonfungi <- Biostrings::readDNAStringSet(idx[sort(short_nonfungi_rows), ])

# One random window per long record, its length drawn from arm C.
window_len <- sample(Biostrings::width(short_nonfungi), length(long))
start <- vapply(
  seq_along(long),
  \(i) sample.int(Biostrings::width(long)[i] - window_len[i] + 1, 1),
  integer(1)
)
long_cut <- Biostrings::subseq(long, start = start, width = window_len)

arms <- list(
  "A mini Fungi" = mini,
  "B + records > 5 kb" = c(mini, long),
  "D + same records cut to typical length" = c(mini, long_cut),
  "C + other non-Fungi records <= 5 kb" = c(mini, short_nonfungi)
)
message(
  length(long),
  " records over 5 kb (median ",
  median(Biostrings::width(long)),
  " bp); cut to a median of ",
  median(Biostrings::width(long_cut)),
  " bp"
)

for (arm in names(arms)) {
  path <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(arms[[arm]], path)
  for (q in names(queries)) {
    res <- dada2::assignTaxonomy(
      queries[[q]],
      refFasta = path,
      minBoot = 0,
      outputBootstraps = TRUE,
      multithread = 2
    )
    boot <- res$boot
    cat(sprintf(
      "%-40s %-9s refs %6d | genus boot >= 50: %5.1f %% | median boot Kingdom %3.0f Genus %3.0f | raw kingdom Fungi %5.1f %%\n",
      arm,
      q,
      length(arms[[arm]]),
      100 * mean(boot[, "Genus"] >= 50),
      median(boot[, "Kingdom"]),
      median(boot[, "Genus"]),
      100 * mean(res$tax[, "Kingdom"] == "Fungi", na.rm = TRUE)
    ))
  }
  unlink(path)
}
