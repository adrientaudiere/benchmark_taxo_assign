# One computation per method, several derived rows (R/assign_compute.R): a
# derived lca or blastn row must equal the direct assignment with the same
# parameters. Skipped without vsearch / BLAST+. Run from the project root:
#   Rscript tests/test_assign_compute.R

library("testthat")
library("here")

here::i_am("tests/test_assign_compute.R")
# config.R: the lca search settings (lca_top_hits_only, lca_id,
# lca_maxaccepts, OD 33) and sintax_seed that compute_assignment() reads.
source(here("config.R"))
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

test_that("lca rows derived from one search equal a direct vsearch --lcaout search at each lca_cutoff", {
  skip_if_not(MiscMetabar::is_vsearch_installed(), "vsearch not available")
  computed <- compute_assignment(pq, "lca", ref_fasta = ref)
  expect_named(computed, c("query", "id", "target"))
  for (cutoff in c(0.51, 1)) {
    suffix <- paste0("_lca_", cutoff)
    derived <- derive_assignment(
      pq,
      computed,
      "lca",
      suffix,
      lca_cutoff = cutoff
    )
    direct <- MiscMetabar::assign_vsearch_lca(
      pq,
      ref_fasta = ref,
      behavior = "add_to_phyloseq",
      suffix = suffix,
      clean_pq = FALSE,
      top_hits_only = lca_top_hits_only,
      id = lca_id,
      maxaccepts = lca_maxaccepts,
      lca_cutoff = cutoff,
      verbose = FALSE
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
