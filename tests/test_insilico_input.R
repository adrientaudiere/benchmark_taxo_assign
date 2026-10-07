# Input records of the in silico set (R/insilico_input.R, D1b, OD 29). Run
# from the project root:
#   Rscript tests/test_insilico_input.R

library("testthat")
library("here")

here::i_am("tests/test_insilico_input.R")
source(here("R/mock_truth.R"))
source(here("R/insilico_input.R"))

set.seed(3)
rand_seq <- function(n) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}
header <- function(name, phylum, genus, species, family = "Fam") {
  paste0(
    name,
    "|ACC|SH|reps|k__Fungi;p__",
    phylum,
    ";c__Cls;o__Ord;f__",
    family,
    ";g__",
    genus,
    ";s__",
    species
  )
}

records <- c(
  header("a1", "Ascomycota", "Alpha", "Alpha_one"),
  header("a2", "Ascomycota", "Alpha", "Alpha_two"),
  header("a3", "Ascomycota", "Alpha", "Alpha_three"),
  header("a4", "Ascomycota", "Alpha", "Alpha_one"),
  header("b1", "Basidiomycota", "Beta", "Beta_one"),
  header("g1", "Glomeromycota", "Gamma", "Gamma_one"),
  header("sp", "Ascomycota", "Alpha", "Alpha_sp"),
  header(
    "is",
    "Ascomycota",
    "Fam_gen_Incertae_sedis",
    "Fam_gen_Incertae_sedis_one"
  ),
  header(
    "fam",
    "Ascomycota",
    "Delta",
    "Delta_one",
    family = "Ord_fam_Incertae_sedis"
  ),
  header("long", "Basidiomycota", "Beta", "Beta_long"),
  header("short", "Basidiomycota", "Beta", "Beta_short"),
  sub("k__Fungi", "k__Metazoa", header("met", "Nematoda", "Eps", "Eps_one"))
)
lengths <- c(400, 400, 400, 400, 400, 400, 400, 400, 400, 900, 150, 400)
seqs <- Biostrings::DNAStringSet(vapply(lengths, rand_seq, character(1)))
names(seqs) <- records

test_that("the three filters and the kingdom filter are applied", {
  out <- names(select_insilico_records(seqs, n = 100))
  expect_false(any(grepl("^(sp|is|fam|long|short|met)\\|", out)))
  expect_true(all(grepl("^(a[1-4]|b1|g1)\\|", out)))
})

test_that("one record per species", {
  out <- names(select_insilico_records(seqs, n = 100))
  species <- sub("^.*s__", "", out)
  expect_equal(anyDuplicated(species), 0L)
  expect_length(out, 5)
})

test_that("phyla are drawn round-robin before a second record of any", {
  out <- names(select_insilico_records(seqs, n = 3))
  phyla <- sub("^.*p__([^;]*);.*$", "\\1", out)
  expect_setequal(phyla, c("Ascomycota", "Basidiomycota", "Glomeromycota"))
})

test_that("the draw depends on the seed only", {
  expect_identical(
    names(select_insilico_records(seqs, n = 4, seed = 1)),
    names(select_insilico_records(seqs, n = 4, seed = 1))
  )
})

test_that("an empty pool stops", {
  expect_error(select_insilico_records(seqs[7:12]), "no record passes")
})
