# Species by exact match (R/species_exact.R, OD 43). Run from the project
# root:
#   Rscript tests/test_species_exact.R

library("testthat")
library("here")

here::i_am("tests/test_species_exact.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
load_pqverse("dbpq")
source(here("R/species_exact.R"))

set.seed(4)
rand_seq <- function(n) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}

test_that("binomial_of_sintax_header() reads UNITE and EUKARYOME species", {
  headers <- c(
    "Abrothallus_subhalei|MT1|SH1|refs;tax=k:Fungi,g:Abrothallus,s:Abrothallus_subhalei",
    "EUK1;tax=k:Fungi,p:Ascomycota,g:Aldinomyces,s:tarquinii",
    "EUK2;tax=k:Fungi,g:Candida,s:Candida_sp",
    "EUK3;tax=k:Fungi,g:unclassified,s:unclassified",
    "EUK4;tax=k:Fungi,g:Mortierella",
    "EUK5;tax=k:Fungi,g:Fusarium,s:Fusarium_incertae_sedis"
  )
  expect_equal(
    binomial_of_sintax_header(headers),
    c("Abrothallus_subhalei", "Aldinomyces_tarquinii", NA, NA, NA, NA)
  )
})

test_that("pool_queries() keeps distinct ACGT sequences and maps them back", {
  shared <- rand_seq(60)
  inputs <- list(
    "pauvert|ITS1_asv" = Biostrings::DNAStringSet(c(
      Taxa_1 = shared,
      Taxa_2 = rand_seq(50),
      Taxa_3 = paste0(rand_seq(20), "N", rand_seq(20))
    )),
    "tedersoo_illumina|ITS1_asv" = Biostrings::DNAStringSet(c(Taxa_9 = shared))
  )
  pooled <- pool_queries(inputs)
  expect_length(pooled$seqs, 2)
  expect_equal(nrow(pooled$map), 3)
  expect_equal(
    pooled$map$query[pooled$map$taxon == "Taxa_1"],
    pooled$map$query[pooled$map$taxon == "Taxa_9"]
  )
  expect_equal(attr(pooled$map, "dropped")$taxon, "Taxa_3")
  expect_true(all(pooled$map$query %in% names(pooled$seqs)))
})

test_that("query_diffs() counts differences over the whole query", {
  q <- rand_seq(100)
  mism <- q
  substr(mism, 50, 50) <- if (substr(q, 50, 50) == "A") "C" else "A"
  ins <- paste0(substr(q, 1, 60), "T", substr(q, 61, 100))
  flank <- rand_seq(30)
  records <- c(
    paste0(flank, q, flank),
    paste0(flank, mism, flank),
    paste0(flank, ins, flank),
    paste0(flank, substr(q, 1, 98))
  )
  out <- query_diffs(
    rep(q, 4),
    records,
    qlo = c(1, 1, 1, 1),
    qhi = c(100, 100, 100, 98),
    tlo = c(31, 31, 31, 31),
    thi = c(130, 130, 131, 128),
    approx = c(0, 1, 1, 2),
    max_diffs = 1
  )
  expect_equal(out, c(0L, 1L, 1L, 2L))
})

test_that("exact_cover() accepts a record that stops early, not one that differs", {
  q <- rand_seq(100)
  other_tail <- paste(
    ifelse(strsplit(substr(q, 96, 100), "")[[1]] == "A", "C", "A"),
    collapse = ""
  )
  records <- c(
    substr(q, 5, 100), # starts after the query: 96 % covered
    substr(q, 15, 100), # 86 % covered
    paste0("A", q, "A"), # holds the query
    paste0(substr(q, 1, 95), other_tail, rand_seq(20)), # last 5 bases differ
    paste0(rand_seq(20), substr(q, 1, 92)), # ends before the query: 92 %
    substr(q, 4, 97) # inside the query: 94 %
  )
  expect_equal(
    exact_cover(rep(q, 6), records, min_cover = 0.9),
    c(TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  )
})

# Six records around one 200 bp query: r1 holds it, r2 with one mismatch,
# r3 with three, r4 ends after its first 185 bases (92.5 % covered at 100 %),
# r5 with an N, r6 holds it under a header ending in "." (BLAST+ drops that
# dot from its ids) and names no species, r7 holds its first 185 bases then
# different ones (not a near match: the record does not stop early).
near_db <- function() {
  q <- rand_seq(200)
  swap <- function(x, at) {
    for (i in at) {
      substr(x, i, i) <- if (substr(x, i, i) == "A") "C" else "A"
    }
    x
  }
  with_n <- q
  substr(with_n, 100, 100) <- "N"
  db <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">r1;tax=k:Fungi,g:Genus,s:Genus_one",
      paste0(rand_seq(100), q, rand_seq(100)),
      ">r2;tax=k:Fungi,g:Genus,s:Genus_two",
      paste0(rand_seq(100), swap(q, 100), rand_seq(100)),
      ">r3;tax=k:Fungi,g:Genus,s:Genus_three",
      paste0(rand_seq(100), swap(q, c(50, 100, 150)), rand_seq(100)),
      ">r4 with space;tax=k:Fungi,g:Genus,s:Genus_four",
      paste0(rand_seq(100), substr(q, 1, 185)),
      ">r5;tax=k:Fungi,g:Genus,s:Genus_five",
      paste0(rand_seq(100), with_n, rand_seq(100)),
      ">r6;tax=k:Fungi,g:Genus,s:sp.nov.inval.",
      paste0(rand_seq(50), q, rand_seq(50)),
      ">r7;tax=k:Fungi,g:Genus,s:Genus_seven",
      paste0(rand_seq(100), swap(q, 186:200), rand_seq(100))
    ),
    db
  )
  list(seqs = Biostrings::DNAStringSet(c(q1 = q)), db = db)
}
near_expected <- tibble::tibble(
  binomial = c("Genus_one", "Genus_two", "Genus_four", "Genus_five", NA),
  diffs = c(0L, 1L, 15L, 1L, 0L),
  exact_cover = c(TRUE, FALSE, TRUE, FALSE, TRUE)
)

test_that("vsearch_near_hits() keeps one difference or 90 % covered at 100 %", {
  skip_if_not(dbpq::is_vsearch_installed())
  x <- near_db()
  hits <- vsearch_near_hits(x$seqs, x$db)
  hits <- hits[order(match(hits$binomial, near_expected$binomial)), ]
  expect_equal(hits[, names(near_expected)], near_expected)
})

test_that("blastn_near_hits() gives the same pairs from its own search", {
  skip_if(!nzchar(Sys.which("blastn")) || !nzchar(Sys.which("makeblastdb")))
  x <- near_db()
  hits <- blastn_near_hits(x$seqs, x$db)
  hits <- hits[order(match(hits$binomial, near_expected$binomial)), ]
  expect_equal(hits[, names(near_expected)], near_expected)
  expect_true(all(hits$n_targets >= 6L))
})

test_that("vsearch_exact_hits() finds contained queries, not ambiguous ones", {
  skip_if_not(dbpq::is_vsearch_installed())
  q <- rand_seq(200)
  with_n <- q
  substr(with_n, 100, 100) <- "N"
  db <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">r1;tax=k:Fungi,g:Genus,s:Genus_one",
      paste0(rand_seq(150), q, rand_seq(250)),
      ">r2;tax=k:Fungi,g:Genus,s:Genus_two",
      paste0(rand_seq(150), with_n, rand_seq(250)),
      ">r3 with space;tax=k:Fungi,g:Genus,s:Genus_three",
      q
    ),
    db
  )
  hits <- vsearch_exact_hits(Biostrings::DNAStringSet(c(q1 = q)), db)
  expect_setequal(hits$binomial, c("Genus_one", "Genus_three"))
  expect_true(all(hits$query == "q1"))
})
