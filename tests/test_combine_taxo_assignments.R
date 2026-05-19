# Verifies combine_taxo_assignments() reproduces the column set and values of
# the sequential `previous_target` chain used in the original
# script_assign_taxo.R. Run from the project root:
#   Rscript tests/test_combine_taxo_assignments.R

library("phyloseq")
library("testthat")
library("here")

here::i_am("tests/test_combine_taxo_assignments.R")
source(here("R/combine_taxo_assignments.R"))

# Build a tiny phyloseq with three ASVs and three base ranks.
make_base_pq <- function() {
  tt <- matrix(
    c("Fungi", "Ascomycota",   "Eurotiomycetes",
      "Fungi", "Basidiomycota", "Agaricomycetes",
      "Fungi", "Ascomycota",   "Sordariomycetes"),
    nrow = 3, byrow = TRUE,
    dimnames = list(
      c("ASV1", "ASV2", "ASV3"),
      c("Kingdom", "Phylum", "Class")
    )
  )
  otu <- matrix(
    c(10, 20, 30, 40, 50, 60),
    nrow = 3,
    dimnames = list(c("ASV1", "ASV2", "ASV3"), c("S1", "S2"))
  )
  phyloseq(tax_table(tt), otu_table(otu, taxa_are_rows = TRUE))
}

# Stand-in for add_new_taxonomy_pq(..., suffix = "_<full_name>"): appends three
# suffixed columns to whatever tax_table is already there.
simulate_add_taxo <- function(pq, full_name,
                              values = c("Fungi", "Ascomycota", "Eurotiomycetes")) {
  current <- as.matrix(unclass(pq@tax_table))
  new_cols <- matrix(
    rep(values, each = nrow(current)),
    nrow = nrow(current),
    dimnames = list(
      rownames(current),
      paste0(c("Kingdom", "Phylum", "Class"), "_", full_name)
    )
  )
  pq@tax_table <- tax_table(cbind(current, new_cols))
  pq
}

test_that("combine_taxo_assignments matches the chain-accumulator output", {
  base_pq <- make_base_pq()
  full_names <- c(
    "dada2__Unite___0.5",
    "sintax__EUK_ITS_v1_9_3___0.5",
    "blastn__Unite___0.5...rel_majority...100"
  )

  # Old: chain — each step adds to the prior.
  chain_pq <- Reduce(simulate_add_taxo, full_names, init = base_pq)

  # New: each assignment built from the same base, then combined.
  assignments <- lapply(full_names, simulate_add_taxo, pq = base_pq)
  combined <- do.call(combine_taxo_assignments, c(list(base_pq), assignments))

  expect_equal(colnames(combined@tax_table), colnames(chain_pq@tax_table))
  expect_equal(
    as.matrix(unclass(combined@tax_table)),
    as.matrix(unclass(chain_pq@tax_table))
  )
})

test_that("combine_taxo_assignments handles a single assignment", {
  base_pq <- make_base_pq()
  pq_a <- simulate_add_taxo(base_pq, "only_method")
  combined <- combine_taxo_assignments(base_pq, pq_a)
  expect_equal(
    as.matrix(unclass(combined@tax_table)),
    as.matrix(unclass(pq_a@tax_table))
  )
})

test_that("combine_taxo_assignments preserves taxa_names and otu_table", {
  base_pq <- make_base_pq()
  pq_a <- simulate_add_taxo(base_pq, "m1")
  pq_b <- simulate_add_taxo(base_pq, "m2")
  combined <- combine_taxo_assignments(base_pq, pq_a, pq_b)

  expect_equal(taxa_names(combined), taxa_names(base_pq))
  expect_equal(as.matrix(otu_table(combined)), as.matrix(otu_table(base_pq)))
})

test_that("combine_taxo_assignments is order-stable (column order = arg order)", {
  base_pq <- make_base_pq()
  pq_a <- simulate_add_taxo(base_pq, "m1")
  pq_b <- simulate_add_taxo(base_pq, "m2")
  pq_c <- simulate_add_taxo(base_pq, "m3")

  combined <- combine_taxo_assignments(base_pq, pq_a, pq_b, pq_c)
  expect_equal(
    colnames(combined@tax_table),
    c("Kingdom", "Phylum", "Class",
      "Kingdom_m1", "Phylum_m1", "Class_m1",
      "Kingdom_m2", "Phylum_m2", "Class_m2",
      "Kingdom_m3", "Phylum_m3", "Class_m3")
  )
})
