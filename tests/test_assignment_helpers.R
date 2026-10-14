# Preparation of the assignment columns before scoring (R/assignment_helpers.R,
# moved out of analysis/01_load_and_clean.qmd on 2026-09-23). Run from the
# project root:
#   Rscript tests/test_assignment_helpers.R

library("testthat")
library("here")
library("phyloseq")

here::i_am("tests/test_assignment_helpers.R")
source(here("R/load_pqverse.R"))
load_pqverse(c("MiscMetabar", "comparpq", "tidypq"))
source(here("R/assignment_helpers.R"))

ranks <- c("Kingdom", "Genus", "Gen_sp")
taxa <- c("t1", "t2", "t3")
tt <- matrix(
  c(
    "Fungi",
    "Fungi",
    "Fungi",
    "Alpha",
    "Fungi_gen_Incertae_sedis",
    "Beta",
    "Alpha_one",
    "Alpha_sp",
    "two",
    "Fungi",
    "Fungi",
    NA,
    "Alpha",
    "unclassified_Beta",
    NA,
    "Alpha_one",
    NA,
    NA
  ),
  nrow = 3,
  dimnames = list(
    taxa,
    c(
      "Kingdom_m__db",
      "Genus_m__db",
      "Species_m__db",
      "Kingdom_n__db",
      "Genus_n__db",
      "Species_n__db"
    )
  )
)
otu <- matrix(1:6, nrow = 3, dimnames = list(taxa, c("s1", "s2")))
physeq <- phyloseq(otu_table(otu, taxa_are_rows = TRUE), tax_table(tt))

test_that("clean_taxo() blanks placeholders, incertae sedis and unclassified", {
  out <- as.data.frame(tax_table(clean_taxo(physeq)), stringsAsFactors = FALSE)
  expect_equal(out$Species_m__db, c("Alpha_one", NA, "two"))
  expect_equal(out$Genus_m__db, c("Alpha", NA, "Beta"))
  expect_equal(out$Genus_n__db, c("Alpha", NA, NA))
})

test_that("gen_sp_for() keeps a binomial and prefixes a bare epithet", {
  binomials <- data.frame(
    Genus_x = c("Alpha", "Beta"),
    Species_x = c("Alpha_one", "Beta_two")
  )
  expect_equal(gen_sp_for(binomials, "x"), c("Alpha_one", "Beta_two"))
  epithets <- data.frame(Genus_x = c("Alpha", NA), Species_x = c("one", "two"))
  expect_equal(gen_sp_for(epithets, "x"), c("Alpha_one", NA))
})

test_that("add_gen_sp() adds one Gen_sp column per assignment", {
  out <- colnames(tax_table(add_gen_sp(clean_taxo(physeq))))
  expect_true(all(c("Gen_sp_m__db", "Gen_sp_n__db") %in% out))
})

test_that("build_ranks_df() maps ranks to columns and guards missing ones", {
  prepared <- add_gen_sp(clean_taxo(physeq))
  rdf <- build_ranks_df(prepared, c("m__db", "n__db"), ranks = ranks)
  expect_equal(dim(rdf), c(3L, 2L))
  expect_equal(rdf$m__db, c("Kingdom_m__db", "Genus_m__db", "Gen_sp_m__db"))
  expect_error(
    build_ranks_df(prepared, c("m__db", "absent__db"), ranks = ranks),
    "missing"
  )
  expect_message(
    rdf <- build_ranks_df(
      prepared,
      c("m__db", "absent__db"),
      ranks = ranks,
      allow_empty = TRUE
    ),
    "1 empty assignments dropped"
  )
  expect_equal(names(rdf), "m__db")
})

test_that("assignment_labels() ignores the ITSx columns", {
  df <- data.frame(
    taxon = "t1",
    Kingdom = "Fungi",
    Genus_m__db = "Alpha",
    Species_m__db = "Alpha_one",
    ITSx_origin = "Fungi",
    ITSx_detected = "TRUE"
  )
  expect_equal(assignment_labels(df), "m__db")
})

test_that("add_consensus_columns() votes over the voters only", {
  # Two of four voters say Alpha: no absolute majority. The bare Genus column
  # and a strategy run before abs_majority must not vote.
  vote_tt <- matrix(
    c("Alpha", "Alpha", "Beta", "Gamma", "Alpha", "Alpha"),
    nrow = 1,
    dimnames = list(
      "t1",
      c(paste0("Genus_in_v", 1:4), "Genus_in_other", "Genus")
    )
  )
  vote_pq <- phyloseq(
    otu_table(
      matrix(1, 1, 1, dimnames = list("t1", "s1")),
      taxa_are_rows = TRUE
    ),
    tax_table(vote_tt)
  )
  configs <- list(
    list(method = "rel_majority", strict = FALSE, nb = 1L, suffix = "rel"),
    list(method = "abs_majority", strict = FALSE, nb = 1L, suffix = "abs")
  )
  out <- suppressMessages(add_consensus_columns(
    vote_pq,
    voter_names = paste0("in_v", 1:4),
    input = "in",
    configs = configs,
    preference_pattern = "v1",
    ranks = "Genus"
  ))
  tt_out <- as(tax_table(out), "matrix")
  expect_equal(unname(tt_out[, "Genus_in_rel"]), "Alpha")
  expect_true(is.na(tt_out[, "Genus_in_abs"]))
  expect_equal(unname(tt_out[, "Genus_in_other"]), "Alpha")
  expect_error(
    suppressMessages(add_consensus_columns(
      vote_pq,
      voter_names = c("in_v1", "in_absent"),
      input = "in",
      configs = configs,
      preference_pattern = "v1",
      ranks = "Genus"
    )),
    "missing"
  )
})
