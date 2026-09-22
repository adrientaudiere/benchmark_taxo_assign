# One computation per method, several derived rows (R/assign_compute.R): a
# derived lca or blastn row must equal the direct assignment with the same
# parameters. Skipped without vsearch / BLAST+. Run from the project root:
#   Rscript tests/test_assign_compute.R

library("testthat")
library("here")

here::i_am("tests/test_assign_compute.R")
source(here("R/load_pqverse.R"))
load_pqverse(c("MiscMetabar", "tidypq"))
source(here("R/assign_compute.R"))

data("data_fungi_mini", package = "MiscMetabar")
ref <- system.file(
  "extdata",
  "mini_UNITE_fungi.fasta.gz",
  package = "MiscMetabar"
)
pq <- phyloseq::prune_taxa(
  phyloseq::taxa_names(data_fungi_mini)[1:30],
  data_fungi_mini
)

new_columns <- function(res, suffix) {
  df <- tidypq::tax_table_to_df(res, convert = FALSE)
  df[, grep(paste0(suffix, "$"), names(df), value = TRUE), drop = FALSE]
}

test_that("lca rows derived from one search equal add_new_taxonomy_pq() at each lca_cutoff", {
  skip_if_not(MiscMetabar::is_vsearch_installed(), "vsearch not available")
  computed <- compute_assignment(pq, "lca", ref_fasta = ref)
  expect_named(computed, c("query", "id", "target"))
  for (cutoff in c(0.8, 1)) {
    suffix <- paste0("_lca_", cutoff)
    derived <- derive_assignment(
      pq,
      computed,
      "lca",
      suffix,
      lca_cutoff = cutoff
    )
    direct <- MiscMetabar::add_new_taxonomy_pq(
      pq,
      ref_fasta = ref,
      method = "lca",
      suffix = suffix,
      clean_pq = FALSE,
      lca_cutoff = cutoff
    )
    expect_equal(new_columns(derived, suffix), new_columns(direct, suffix))
  }
})

test_that("blastn rows apply their min_id to the stored hits", {
  skip_if(Sys.which("blastn") == "", "BLAST+ not available")
  computed <- compute_assignment(pq, "blastn", ref_fasta = ref)
  suffix <- "_blastn_90"
  derived <- derive_assignment(
    pq,
    computed,
    "blastn",
    suffix,
    vote_algorithm = "rel_majority",
    nb_voting = 100,
    min_cover = 80,
    min_id = 90
  )
  direct <- MiscMetabar::assign_blastn(
    pq,
    ref_fasta = ref,
    behavior = "add_to_phyloseq",
    suffix = suffix,
    vote_algorithm = "rel_majority",
    nb_voting = 100,
    min_cover = 80,
    min_id = 90
  )
  expect_equal(new_columns(derived, suffix), new_columns(direct, suffix))
  strict <- derive_assignment(
    pq,
    computed,
    "blastn",
    "_blastn_97",
    vote_algorithm = "rel_majority",
    nb_voting = 100,
    min_cover = 80,
    min_id = 97
  )
  # A stricter min_id assigns fewer values (none at all adds no column).
  expect_lte(
    sum(!is.na(new_columns(strict, "_blastn_97"))),
    sum(!is.na(new_columns(derived, suffix)))
  )
})
