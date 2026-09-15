# OTU taxonomy from the member ASVs (R/otu_taxonomy.R, ROADMAP Q5). The
# resolution is tested on a hand-made membership; the vsearch clustering is
# skipped when vsearch is not on the PATH. Run from the project root:
#   Rscript tests/test_otu_taxonomy.R

library("testthat")
library("here")
library("phyloseq")

here::i_am("tests/test_otu_taxonomy.R")
source(here("R/load_pqverse.R"))
load_pqverse(c("MiscMetabar", "tidypq"))
source(here("R/otu_taxonomy.R"))

seq_a <- paste(rep("ACGTTGACCATGGTACCAGTTAGCCATAGG", 8), collapse = "")
seq_b <- sub("^A", "G", seq_a)
seq_c <- paste(rep("TTGCAAGGCTTACCGATTAGCATCGGTAAC", 8), collapse = "")
seq_fake <- paste(rep("GGATCCTTAAGCTAGCAATTGCCGGTTAAC", 8), collapse = "")

tiny_pq <- function() {
  taxa <- c("Taxa_1", "Taxa_2", "Taxa_3", "fake_1")
  otu <- otu_table(
    matrix(
      c(10, 2, 5, 0, 4, 1, 3, 0),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("s1", "s2"), taxa)
    ),
    taxa_are_rows = FALSE
  )
  tax <- tax_table(matrix(
    c(
      "G1",
      "S1",
      "G1",
      "G1",
      "S2",
      NA,
      "G2",
      "S3",
      "G2",
      NA,
      NA,
      NA
    ),
    nrow = 4,
    byrow = TRUE,
    dimnames = list(taxa, c("Genus_x", "Species_x", "Genus_y"))
  ))
  refs <- Biostrings::DNAStringSet(c(
    Taxa_1 = seq_a,
    Taxa_2 = seq_b,
    Taxa_3 = seq_c,
    fake_1 = seq_fake
  ))
  phyloseq(otu, tax, refseq(refs))
}

test_that("resolve_otu_taxonomy keeps the values all member ASVs agree on", {
  membership <- tibble::tibble(
    taxon = c("Taxa_1", "Taxa_2", "Taxa_3"),
    cluster = c(0L, 0L, 1L),
    archetype = c("Taxa_1", "Taxa_1", "Taxa_3")
  )
  out <- resolve_otu_taxonomy(tiny_pq(), membership)
  expect_setequal(taxa_names(out), c("Taxa_1", "Taxa_3", "fake_1"))
  tt <- tidypq::tax_table_to_df(out, convert = FALSE)
  otu_1 <- tt[tt$taxon == "Taxa_1", ]
  expect_equal(otu_1$Genus_x, "G1")
  expect_true(is.na(otu_1$Species_x))
  expect_equal(otu_1$Genus_y, "G1")
  expect_equal(tt$Species_x[tt$taxon == "Taxa_3"], "S3")
  expect_true(all(is.na(unlist(tt[tt$taxon == "fake_1", -1]))))
  expect_equal(
    unname(taxa_sums(out)[c("Taxa_1", "Taxa_3", "fake_1")]),
    c(17, 8, 0)
  )
  expect_equal(as.character(out@refseq[["Taxa_1"]]), seq_a)
})

test_that("otu_membership clusters near-identical sequences with vsearch", {
  skip_if(Sys.which("vsearch") == "", "vsearch not on the PATH")
  pq <- prune_taxa(c("Taxa_1", "Taxa_2", "Taxa_3"), tiny_pq())
  membership <- otu_membership(pq)
  expect_setequal(membership$taxon, taxa_names(pq))
  cl <- stats::setNames(membership$cluster, membership$taxon)
  expect_equal(cl[["Taxa_1"]], cl[["Taxa_2"]])
  expect_false(cl[["Taxa_1"]] == cl[["Taxa_3"]])
  expect_equal(
    unique(membership$archetype[membership$cluster == cl[["Taxa_1"]]]),
    "Taxa_1"
  )
})
