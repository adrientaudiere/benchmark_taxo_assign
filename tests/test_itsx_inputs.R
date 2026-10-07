# The ITSx inputs of the grid (R/itsx_inputs.R, OD 34) and the extraction by
# MiscMetabar::itsx_pq() that replaced the former R/itsx.R prototype. The
# expected sequences were checked against ITSx 1.1.3 when the prototype was
# written. Skipped when the itsxenv conda env is missing. Run from the project
# root:
#   Rscript tests/test_itsx_inputs.R

library("testthat")
library("here")
library("phyloseq")

here::i_am("tests/test_itsx_inputs.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
load_pqverse(c("MiscMetabar", "comparpq"))
source(here("R/itsx_inputs.R"))

itsx_ok <- MiscMetabar::is_itsx_installed(args_before_itsx = itsx_conda_prelude)

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

tiny_pq <- function(seqs, counts) {
  otu <- otu_table(
    matrix(counts, nrow = 1, dimnames = list("s1", names(seqs))),
    taxa_are_rows = FALSE
  )
  tax <- tax_table(matrix(
    "Fungi",
    nrow = length(seqs),
    dimnames = list(names(seqs), "Kingdom")
  ))
  # comparpq's control builders need a sample_data slot with 2 variables.
  sam <- sample_data(data.frame(site = "a", depth = 1, row.names = "s1"))
  phyloseq(otu, tax, sam, refseq(Biostrings::DNAStringSet(seqs)))
}

test_that("itsx_pq() cuts the ITS1 and the full span exactly as the former prototype", {
  skip_if_not(itsx_ok, "ITSx (itsxenv conda env) not available")
  its1 <- MiscMetabar::itsx_pq(
    tiny_pq(c(Taxa_1 = asv_1, Taxa_2 = not_its), c(5, 3)),
    region = "ITS1",
    args_before_itsx = itsx_conda_prelude,
    verbose = FALSE
  )
  expect_equal(as.character(its1@refseq[["Taxa_1"]]), its1_1)
  expect_equal(as.character(its1@refseq[["Taxa_2"]]), not_its)
  expect_true(is.na(its1@tax_table["Taxa_2", "ITSx_origin"]))
  full <- MiscMetabar::itsx_pq(
    tiny_pq(c(Taxa_1 = asv_full_its), 5),
    region = "full",
    args_before_itsx = itsx_conda_prelude,
    verbose = FALSE
  )
  # ITSx ends the full span one base after the ITS2 end of its positions table.
  expect_equal(
    as.character(full@refseq[["Taxa_1"]]),
    substr(asv_full_its, 155, 657)
  )
  expect_equal(as.character(full@tax_table["Taxa_1", "ITSx_origin"]), "Fungi")
})

test_that("with_controls() adds fakes from the trimmed units and trimmed externals with their origin", {
  skip_if_not(itsx_ok, "ITSx (itsxenv conda env) not available")
  real <- MiscMetabar::itsx_pq(
    tiny_pq(c(Taxa_1 = asv_1, Taxa_2 = asv_full_its), c(5, 3)),
    region = "ITS1",
    add_detected = TRUE,
    args_before_itsx = itsx_conda_prelude,
    verbose = FALSE
  )
  ext <- itsx_controls(
    Biostrings::DNAStringSet(c(ctrl_its = asv_full_its, ctrl_random = not_its)),
    region = "ITS1",
    args_before_itsx = itsx_conda_prelude
  )
  expect_setequal(names(ext$seqs), c("ctrl_its", "ctrl_random"))
  expect_equal(unname(ext$origin[["ctrl_its"]]), "Fungi")
  expect_true(is.na(ext$origin[["ctrl_random"]]))
  expect_equal(
    unname(ext$detected[c("ctrl_its", "ctrl_random")]),
    c("TRUE", "FALSE")
  )
  expect_equal(
    as.character(ext$seqs[["ctrl_its"]]),
    substr(asv_full_its, 155, 316)
  )

  out <- with_controls(real, ext, prop_fake = 0.5, seed = 22)
  names_out <- taxa_names(out)
  expect_equal(sum(grepl("^fake_", names_out)), 1)
  expect_setequal(
    grep("^external_", names_out, value = TRUE),
    c("external_ctrl_its", "external_ctrl_random")
  )
  expect_equal(
    as.character(out@tax_table["external_ctrl_its", "ITSx_origin"]),
    "Fungi"
  )
  expect_equal(
    as.character(out@tax_table[
      c("external_ctrl_its", "external_ctrl_random"),
      "ITSx_detected"
    ]),
    c("TRUE", "FALSE")
  )
  fake <- grep("^fake_", names_out, value = TRUE)
  expect_true(all(is.na(out@tax_table[fake, ])))
  expect_true(
    Biostrings::width(out@refseq[fake]) %in%
      Biostrings::width(real@refseq)
  )
  expect_equal(
    out@tax_table[c("Taxa_1", "Taxa_2"), "ITSx_origin"],
    real@tax_table[c("Taxa_1", "Taxa_2"), "ITSx_origin"]
  )
  again <- with_controls(real, ext, prop_fake = 0.5, seed = 22)
  expect_identical(
    as.character(again@refseq[fake]),
    as.character(out@refseq[fake])
  )
})
