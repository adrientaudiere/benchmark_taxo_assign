# Pins the target names produced by R/values_map.R. Renaming a target silently
# orphans its store object, so any change here must be deliberate.
# Run from the project root:
#   Rscript tests/test_values_map.R

library("testthat")
library("here")
library("dplyr")

here::i_am("tests/test_values_map.R")
source(here("config.R"))
source(here("R/values_map.R"))

# The formula used by the pipelines before S2.1 (kept verbatim as the oracle).
legacy_full_name <- function(method, db, min_bootstrap, vote_algorithm, nb_voting) {
  gsub(
    "...NA...NA", "",
    paste0(method, "__", db, "___", min_bootstrap, "...", vote_algorithm, "...", nb_voting)
  )
}

test_that("build_values_map reproduces the legacy target names", {
  vm <- build_values_map(mini_db = FALSE)
  expect_equal(nrow(vm), 12 * length(db_list))
  expect_false(any(duplicated(vm$full_name)))
  expect_equal(
    vm$full_name,
    legacy_full_name(vm$method, vm$db, vm$min_bootstrap, vm$vote_algorithm, vm$nb_voting)
  )
  expect_true(all(c(
    "dada2__Unite___0.5",
    "sintax__EUK_ITS_v2_Fungi_cut___0.4",
    "blastn__EUK_SSU_v2_Fungi_cut___0.5...rel_majority...100"
  ) %in% vm$full_name))
})

test_that("db_path follows the method and the mini_db flag", {
  vm <- build_values_map(dbs = "Unite", mini_db = FALSE)
  expect_equal(unique(vm$db_path[vm$method == "dada2"]),
               "data/data_raw/refseq/dada2_format/Unite.fasta")
  expect_equal(unique(vm$db_path[vm$method == "sintax"]),
               "data/data_raw/refseq/sintax_format/Unite.fasta")
  vm_mini <- build_values_map(dbs = "Unite", mini_db = TRUE)
  expect_true(all(grepl("/mini_Unite\\.fasta$", vm_mini$db_path)))
  expect_equal(vm_mini$full_name, vm$full_name)
})

test_that("build_cv_values_map reproduces the legacy CV target names", {
  cv <- build_cv_values_map(mini_db = FALSE)
  expect_equal(nrow(cv), 4 * length(db_list) * 2)
  expect_false(any(duplicated(cv$full_name)))
  expect_true(all(c(
    "cv__dada2__Unite__standard",
    "cv__blastn__EUK_SSU_v2_Fungi_cut__leaked"
  ) %in% cv$full_name))
  expect_equal(unique(cv$vote_algorithm[cv$method == "blastn"]), "rel_majority")
  expect_true(all(is.na(cv$vote_algorithm[cv$method != "blastn"])))
})

test_that("db_meta covers every database of db_list", {
  expect_setequal(db_meta$db, db_list)
})
