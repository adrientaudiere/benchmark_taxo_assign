# ITSx helpers of R/itsx.R. The replacement logic is tested without ITSx; the
# ITSx run is skipped when the itsxenv conda env is missing. Run from the
# project root:
#   Rscript tests/test_itsx.R

library("testthat")
library("here")
library("phyloseq")

here::i_am("tests/test_itsx.R")
source(here("config.R"))
source(here("R/itsx.R"))

# Taxa_1 of the mock-community d_asv: 45 bp of 18S, ITS1, start of 5.8S.
asv_1 <- paste0(
  "AGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTACAGAGTTCATGCCCGAAAGGGTAGACCTCCCACCC",
  "TTGTGTATTATTACTTTGTTGCTTTGGCGAGCTGCCTTCGGGCCTTGTATGCTCGCCAGAGAATACCAAAACTCTTTTTA",
  "TTAATGTCGTCTGAGTACTATATAATAGTTAAAACTTTCAACAACGGATCTCTTGGTTCT"
)
its1_1 <- paste0(
  "CAGAGTTCATGCCCGAAAGGGTAGACCTCCCACCCTTGTGTATTATTACTTTGTTGCTTTGGCGAGCTGCCTTCGGGCCTTG",
  "TATGCTCGCCAGAGAATACCAAAACTCTTTTTATTAATGTCGTCTGAGTACTATATAATAGTTA"
)
not_its <- paste(rep("ACGTTGCA", 30), collapse = "")
# Taxa_1 of the Tedersoo Illumina mock d_asv (full-ITS amplicon): ITSx 1.1.3
# finds SSU 1-154, ITS1 155-316, 5.8S 317-474, ITS2 475-656, LSU 657-695.
asv_full_its <- paste0(
  "TACTACCGATTGAATGGCTTAGTGAGGTCTCCGGATTAGCTTTGGCGCACCGGCAACGGAATGCTATTGCTGAGAAGTTG",
  "ATCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTATTGAAA",
  "TAAACCTGATGGGTTGTTGCTGGTTCTCTAGGGAGCATGTGCACACCTTGTCATCTTTATATCTCCACCTGTGCACCTTT",
  "TGTAGACCTGAAAGGTCTATGTTGCTTCATTTACCCCAATGTATGTCAATAGAATGTTGTGCCTATATAATATATACAAC",
  "TTTCAGCAACGGATCTCTTGGCTCTCGCATCGATGAAGAACGCAGCGAAATGCGATAAGTAATGTGAATTGCAGAATTCA",
  "GTGAATCATCGAATCTTTGAACGCACCTTGCGCTCCTTGGTATTCCGAGGAGCATGCCTGTTTGAGTGTCATTAATATAT",
  "CAACCTTCTCTTTTTGAGTGGTTTGGATGTGGGGGTTTGCTGGCCTCTTAAAAGGTCTTGGCTCTCCTGAAATACATTAG",
  "CAGAACAACCCTGTTCATTGGTGTGATAACTATCTACGCTATTGAATGTGAAGGGCAGTTTTGCTTTCTAACAGTCCTTG",
  "GACAAGCTCATCATTAATGTGACCTCAAATCAGGTAGGACTACCCGCTGAACTTA"
)

tiny_pq <- function() {
  otu <- otu_table(
    matrix(c(5, 3), nrow = 1, dimnames = list("s1", c("Taxa_1", "Taxa_2"))),
    taxa_are_rows = FALSE
  )
  phyloseq(
    otu,
    refseq(Biostrings::DNAStringSet(c(Taxa_1 = asv_1, Taxa_2 = not_its)))
  )
}

test_that("itsx_replace_refseq swaps detected sequences and keeps the taxa set", {
  itsx <- list(
    sequences = Biostrings::DNAStringSet(c(Taxa_1 = its1_1)),
    region = "ITS1"
  )
  expect_message(out <- itsx_replace_refseq(tiny_pq(), itsx), "1 of 2 taxa")
  expect_equal(taxa_names(out), c("Taxa_1", "Taxa_2"))
  expect_equal(as.character(out@refseq[["Taxa_1"]]), its1_1)
  expect_equal(as.character(out@refseq[["Taxa_2"]]), not_its)
})

test_that("itsx_replace_refseq can drop undetected taxa", {
  itsx <- list(
    sequences = Biostrings::DNAStringSet(c(Taxa_1 = its1_1)),
    region = "ITS1"
  )
  out <- suppressMessages(itsx_replace_refseq(
    tiny_pq(),
    itsx,
    keep_undetected = FALSE
  ))
  expect_equal(taxa_names(out), "Taxa_1")
  expect_equal(as.character(out@refseq[["Taxa_1"]]), its1_1)
})

test_that("taxa identical after trimming: the less abundant is dropped from both inputs", {
  # Taxa_1 and Taxa_2 differ only in the 18S flank; both trim to its1_1.
  otu <- otu_table(
    matrix(
      c(3, 5, 1),
      nrow = 1,
      dimnames = list("s1", c("Taxa_1", "Taxa_2", "Taxa_3"))
    ),
    taxa_are_rows = FALSE
  )
  pq <- phyloseq(
    otu,
    refseq(Biostrings::DNAStringSet(
      c(Taxa_1 = asv_1, Taxa_2 = paste0("GGGG", asv_1), Taxa_3 = not_its)
    ))
  )
  itsx <- list(
    sequences = Biostrings::DNAStringSet(c(Taxa_1 = its1_1, Taxa_2 = its1_1)),
    region = "ITS1"
  )
  expect_error(itsx_replace_refseq(pq, itsx), "itsx_duplicated_taxa")
  expect_message(dropped <- itsx_duplicated_taxa(pq, itsx), "1 taxa removed")
  expect_equal(dropped, "Taxa_1")
  kept <- prune_taxa(setdiff(taxa_names(pq), dropped), pq)
  out <- itsx_replace_refseq(kept, itsx)
  expect_equal(taxa_names(out), c("Taxa_2", "Taxa_3"))
  expect_equal(as.character(out@refseq[["Taxa_2"]]), its1_1)
  itsx_one <- list(
    sequences = Biostrings::DNAStringSet(c(Taxa_1 = its1_1)),
    region = "ITS1"
  )
  expect_length(itsx_duplicated_taxa(tiny_pq(), itsx_one), 0)
})

test_that("run_itsx extracts ITS1 from a real ASV and ignores a non-ITS sequence", {
  itsx_ok <- system2(
    "bash",
    c("-c", shQuote(paste0(itsx_conda_prelude, "ITSx --help"))),
    stdout = FALSE,
    stderr = FALSE
  ) ==
    0
  skip_if_not(itsx_ok, "ITSx (itsxenv conda env) not available")
  res <- run_itsx(
    Biostrings::DNAStringSet(c(Taxa_1 = asv_1, Taxa_2 = not_its)),
    region = itsx_region,
    organism_groups = itsx_organism_groups,
    prelude = itsx_conda_prelude
  )
  expect_equal(names(res$sequences), "Taxa_1")
  expect_equal(as.character(res$sequences[["Taxa_1"]]), its1_1)
  expect_true("Taxa_1" %in% res$positions$taxon)
  expect_equal(
    res$positions$ITS1[res$positions$taxon == "Taxa_1"],
    "ITS1: 46-191"
  )
})

test_that("run_itsx region = 'full' returns ITS1 + 5.8S + ITS2 of a full-ITS ASV", {
  itsx_ok <- system2(
    "bash",
    c("-c", shQuote(paste0(itsx_conda_prelude, "ITSx --help"))),
    stdout = FALSE,
    stderr = FALSE
  ) ==
    0
  skip_if_not(itsx_ok, "ITSx (itsxenv conda env) not available")
  res <- run_itsx(
    Biostrings::DNAStringSet(c(Taxa_1 = asv_full_its)),
    region = "full",
    organism_groups = itsx_organism_groups,
    prelude = itsx_conda_prelude
  )
  expect_equal(res$region, "full")
  # ITSx ends the full span one base after the ITS2 end of its positions table.
  expect_equal(
    as.character(res$sequences[["Taxa_1"]]),
    substr(asv_full_its, 155, 657)
  )
})
