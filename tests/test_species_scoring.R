# Species calls from exact hits (R/species_scoring.R, OD 43). Run from the
# project root:
#   Rscript tests/test_species_scoring.R

library("testthat")
library("here")

here::i_am("tests/test_species_scoring.R")
source(here("R/species_scoring.R"))

hits <- tibble::tibble(
  db = "D",
  query = c("q1", "q1", "q2", "q2", "q3", "q4", "q4"),
  binomial = c("A_b", "A_b", "A_b", "A_c", NA, "A_b", NA)
)

test_that("calls_from_hits applies the unanimity rule, unnamed hits aside", {
  out <- calls_from_hits(hits, "vsearch")
  out <- out[order(out$query), ]
  expect_equal(
    out$status,
    c("one species", "several species", "only unnamed hits", "one species")
  )
  expect_equal(out$call, c("A_b", NA, NA, "A_b"))
  expect_equal(out$n_binomials, c(1L, 2L, 0L, 1L))
})

test_that("calls_from_assign_species reads Genus and the joined epithets", {
  assigned <- tibble::tibble(
    db = "D",
    query = c("q1", "q2", "q3", "q4"),
    genus = c("A", "A", NA, NA),
    epithets = c("b", "b/c", "b/d", NA)
  )
  out <- calls_from_assign_species(assigned)
  expect_equal(out$query, c("q1", "q2", "q3"))
  expect_equal(out$call, c("A_b", NA, NA))
  expect_equal(
    out$status,
    c("one species", "several species", "several species")
  )
})

test_that("species_calls gives every query of every method a status", {
  assigned <- tibble::tibble(
    db = "D",
    query = "q1",
    genus = "A",
    epithets = "b"
  )
  out <- species_calls(
    assigned,
    hits,
    hits[0, ],
    queries = c("q1", "q2", "q3", "q4", "q5")
  )
  expect_equal(nrow(out), 15)
  blastn <- out[out$method == "blastn", ]
  expect_true(all(blastn$status == "no hit"))
  expect_equal(
    as.character(out$status[out$method == "vsearch" & out$query == "q5"]),
    "no hit"
  )
})

test_that("near_calls applies the rule at each match level", {
  near <- tibble::tibble(
    db = "D",
    query = c("q1", "q1", "q2"),
    binomial = c("A_b", "A_c", "A_b"),
    diffs = c(0L, 1L, 3L),
    exact_cover = c(TRUE, FALSE, TRUE)
  )
  out <- near_calls(near, near[0, ], queries = c("q1", "q2"))
  expect_equal(nrow(out), 12)
  v <- out[out$method == "vsearch", ]
  v <- v[order(v$query, v$match), ]
  expect_equal(v$call, c("A_b", NA, "A_b", NA, NA, "A_b"))
  expect_equal(
    as.character(v$status),
    c(
      "one species",
      "several species",
      "one species",
      "no hit",
      "no hit",
      "one species"
    )
  )
  expect_true(all(out$status[out$method == "blastn"] == "no hit"))
})

test_that("add_exact_columns writes the calls on every taxon of a query", {
  skip_if_not_installed("phyloseq")
  otu <- phyloseq::otu_table(
    matrix(1, 3, 1, dimnames = list(c("t1", "t2", "t3"), "s")),
    taxa_are_rows = TRUE
  )
  tax <- phyloseq::tax_table(matrix(
    "x",
    3,
    1,
    dimnames = list(c("t1", "t2", "t3"), "Kingdom")
  ))
  physeq <- phyloseq::phyloseq(otu, tax)
  calls <- tibble::tibble(
    db = "D",
    query = c("q1", "q2"),
    method = "vsearch",
    call = c("A_b", NA)
  )
  map <- tibble::tibble(
    taxon = c("t1", "t3", "t2"),
    query = c("q1", "q1", "q2")
  )
  out <- add_exact_columns(physeq, calls, map, "ITS1_asv")
  near <- add_exact_columns(
    physeq,
    dplyr::mutate(calls, match = "90 % covered"),
    map,
    "ITS1_asv"
  )
  expect_true(
    "Gen_sp_ITS1_asv_vsearch__D___cover90" %in%
      colnames(phyloseq::tax_table(near$physeq))
  )
  col <- "Gen_sp_ITS1_asv_vsearch__D___exact"
  tt <- methods::as(phyloseq::tax_table(out$physeq), "matrix")
  expect_equal(unname(tt[, col]), c("A_b", NA, "A_b"))
  expect_true(all(is.na(tt[, "no_rank"])))
  expect_equal(out$ranks_df[[1]], c(rep("no_rank", 6), col))
})
