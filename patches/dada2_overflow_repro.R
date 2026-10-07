# Minimal reproducible example of the dada2 assignTaxonomy() overflow
# (patches/dada2_1.40.0_taxonomy_overflow.patch), for the upstream report.
# A reference of 65 536 random "Dummy" sequences, each with its own taxonomy,
# followed by one "Real;Target" sequence (the 65 537th taxonomy). The query is
# that sequence itself. About 17 GB of RAM (65 537 x 4^8 floats).
#   Rscript --no-init-file patches/dada2_overflow_repro.R   (Bioconductor dada2)
#   Rscript patches/dada2_overflow_repro.R                  (patched, via .Rprofile)

library(dada2)
message("dada2 ", packageVersion("dada2"), " from ", find.package("dada2"))

set.seed(1)
random_seq <- function(n, len = 300) {
  vapply(
    seq_len(n),
    \(i) {
      paste(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = "")
    },
    character(1)
  )
}
n_dummy <- 65536
target <- random_seq(1)
ref <- Biostrings::DNAStringSet(c(random_seq(n_dummy), target))
names(ref) <- c(
  sprintf("Dummy;D%05d;", seq_len(n_dummy)),
  "Real;Target;"
)
ref_file <- tempfile(fileext = ".fasta")
Biostrings::writeXStringSet(ref, ref_file)

res <- assignTaxonomy(
  target,
  refFasta = ref_file,
  taxLevels = c("Kingdom", "Genus"),
  minBoot = 0,
  outputBootstraps = TRUE
)
print(data.frame(taxonomy = res$tax[1, ], bootstrap = res$boot[1, ]))
# Expected: Real / Target with bootstrap 100 (the query is identical to it).
# Bioconductor dada2 1.40.0: a Dummy taxonomy with a low bootstrap, because
# taxonomy 65 537 is scored with the k-mer profile of taxonomy 1.
