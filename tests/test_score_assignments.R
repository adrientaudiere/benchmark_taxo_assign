# Grouping of the assignment columns by the treatment of the external controls
# (R/score_assignments.R, ROADMAP 0.8a), and the scoring call it wraps. Run
# from the project root:
#   Rscript tests/test_score_assignments.R

library("testthat")
library("here")
library("phyloseq")

here::i_am("tests/test_score_assignments.R")
source(here("R/load_pqverse.R"))
load_pqverse("comparpq")
source(here("R/score_assignments.R"))

db_meta_test <- tibble::tibble(
  db = c(
    "Unite_all_20250219",
    "Unite_all_20250219_Fungi",
    "Unite_all_20250219_Fungi_rep"
  ),
  db_base = "Unite_all_20250219",
  simplification = c("full", "Fungi", "Fungi+rep")
)

test_that("assignment_db() reads the database out of a target name", {
  expect_equal(
    assignment_db("dada2__Unite_all_20250219_Fungi___0.5"),
    "Unite_all_20250219_Fungi"
  )
  expect_equal(
    assignment_db("blastn__Unite_all_20250219___0.5...rel_majority...100...95"),
    "Unite_all_20250219"
  )
  expect_equal(
    assignment_db("otu_lca__Unite_all_20250219_Fungi___1"),
    "Unite_all_20250219_Fungi"
  )
  # A `mini_` prefix is stripped as a safety net: the *_mini projects name
  # their targets after the production database (`full_name_for()` takes `db`,
  # not `db_name`), so it should never appear here.
  expect_equal(
    assignment_db("lca__mini_Unite_all_20250219_Fungi___1"),
    "Unite_all_20250219_Fungi"
  )
  expect_true(is.na(assignment_db("rel_majority_consensus")))
})

test_that("each database gets the scoring its simplification imposes", {
  groups <- scoring_groups(
    c(
      "dada2__Unite_all_20250219___0.5",
      "sintax__Unite_all_20250219_Fungi___0.5",
      "lca__Unite_all_20250219_Fungi_rep___1"
    ),
    db_meta_test
  )
  expect_equal(nrow(groups), 3)
  expect_equal(
    groups$external_scoring,
    c("aside", "matrix", "aside")
  )
  expect_equal(groups$kingdom_only, c(FALSE, FALSE, TRUE))
})

test_that("a consensus column is scored once per requested scoring", {
  groups <- scoring_groups(
    c("dada2__Unite_all_20250219_Fungi___0.5", "unanimity_consensus"),
    db_meta_test
  )
  expect_equal(sum(groups$assignment == "unanimity_consensus"), 2)
  consensus <- groups[groups$assignment == "unanimity_consensus", ]
  expect_setequal(consensus$external_scoring, c("matrix", "aside"))
  # ext_correct stays at the kingdom for the "aside" reading only.
  expect_equal(
    consensus$kingdom_only[consensus$external_scoring == "aside"],
    TRUE
  )
  one <- scoring_groups(
    "unanimity_consensus",
    db_meta_test,
    consensus_scorings = "matrix"
  )
  expect_equal(nrow(one), 1)
})

test_that("an unknown database is an error, not a silent NA", {
  expect_error(
    scoring_groups("dada2__Not_a_database___0.5", db_meta_test),
    "not in db_meta"
  )
})

# A hand-made object: two real units with a truth, one shuffled and one
# external control, and one assignment column per simplification level.
score_mock <- function() {
  taxa <- c("Taxa_1", "Taxa_2", "fake_1", "external_1")
  otu <- otu_table(
    matrix(1, nrow = 1, ncol = 4, dimnames = list("s1", taxa)),
    taxa_are_rows = FALSE
  )
  tax <- tax_table(as.matrix(data.frame(
    row.names = taxa,
    Kingdom_dada2__Unite_all_20250219___0.5 = c(
      "Fungi",
      "Fungi",
      NA,
      "Viridiplantae"
    ),
    Genus_dada2__Unite_all_20250219___0.5 = c("Alpha", "Beta", NA, "Pinus"),
    Kingdom_dada2__Unite_all_20250219_Fungi___0.5 = c(
      "Fungi",
      "Fungi",
      NA,
      "Fungi"
    ),
    Genus_dada2__Unite_all_20250219_Fungi___0.5 = c("Alpha", "Beta", NA, "Zeta")
  )))
  phyloseq(otu, tax)
}

score_truth <- data.frame(
  taxon = c("Taxa_1", "Taxa_2"),
  Kingdom = "Fungi",
  Genus = c("Alpha", "Alpha"),
  truth_depth = "Genus"
)

score_external <- data.frame(
  taxon = "external_1",
  Kingdom = "Viridiplantae",
  Genus = "Pinus"
)

test_that("score_assignments() scores each column under its own rule", {
  assignments <- c(
    "dada2__Unite_all_20250219___0.5",
    "dada2__Unite_all_20250219_Fungi___0.5"
  )
  ranks_df <- as.data.frame(lapply(
    stats::setNames(assignments, assignments),
    \(fn) paste0(c("Kingdom", "Genus"), "_", fn)
  ))
  res <- score_assignments(
    score_mock(),
    ranks_df = ranks_df,
    truth = score_truth,
    db_meta = db_meta_test,
    external_truth = score_external,
    truth_ranks = c("Kingdom", "Genus")
  )
  cell <- function(assignment, metric) {
    res$values[
      res$method_db == assignment &
        res$tax_level == "Genus" &
        res$metrics == metric
    ]
  }
  # The whole release can hold the control: it stays out of the matrix
  # (2 real + 1 shuffled = 3 units), the Fungi database puts it in (4).
  expect_equal(
    sum(unlist(lapply(c("TP", "FP", "FN", "TN"), \(m) {
      cell(assignments[[1]], m)
    }))),
    3
  )
  expect_equal(
    sum(unlist(lapply(c("TP", "FP", "FN", "TN"), \(m) {
      cell(assignments[[2]], m)
    }))),
    4
  )
  expect_equal(
    unique(res$external_scoring[res$method_db == assignments[[1]]]),
    "aside"
  )
  expect_equal(
    unique(res$external_scoring[res$method_db == assignments[[2]]]),
    "matrix"
  )
  # ext_correct is computed at every rank for the whole release, and the
  # control keeps its own lineage there.
  expect_equal(
    res$values[
      res$method_db == assignments[[1]] &
        res$metrics == "ext_correct" &
        res$tax_level == "Genus"
    ],
    1
  )
  # ext_fungi sees the Fungi database calling the control Fungi.
  expect_equal(
    res$values[
      res$method_db == assignments[[2]] & res$metrics == "ext_fungi"
    ],
    1
  )
})
