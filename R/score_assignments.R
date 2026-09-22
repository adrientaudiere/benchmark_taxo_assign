# Scoring of the assignment columns against the per-unit truth of a mock
# (ROADMAP 0.8a): groups the columns by the treatment their database imposes on
# the `external_` controls, then calls comparpq::tc_metrics_unit() once per
# group.
#
# Why the grouping (docs/objectives_design.md decisions 6, 15, 22): naming an
# external control is an error only when the reference database cannot hold it.
#
#   simplification | external_scoring | ext_correct
#   ---------------|------------------|-------------
#   Fungi          | "matrix"         | every rank
#   full           | "aside"          | every rank
#   Fungi+rep      | "aside"          | Kingdom only
#
# A consensus column has no database: it votes over the nine, so both readings
# are defensible and the developer asked (2026-09-22) for both — every
# consensus column is scored twice and the rows are told apart by the
# `external_scoring` column of the result.
#
# Requires comparpq (tc_metrics_unit) and dplyr.

# Database an assignment column was computed with, from its name
# ("dada2__Unite_all_20250219_Fungi___0.5" -> "Unite_all_20250219_Fungi").
# The `mini_` prefix of the smoke-test databases is dropped so the result joins
# `config.R::db_meta`. A consensus column carries no "__" and gets NA.
assignment_db <- function(assignment) {
  db <- sub(".*__", "", sub("___.*", "", assignment))
  db <- sub("^mini_", "", db)
  ifelse(grepl("__", assignment), db, NA_character_)
}

# One row per assignment column: the `external_scoring` it must be scored with
# and whether `ext_correct` is restricted to the kingdom. Consensus columns
# (database NA) get one row per element of `consensus_scorings`.
scoring_groups <- function(
  assignments,
  db_meta,
  consensus_scorings = c("matrix", "aside")
) {
  rules <- tibble::tibble(
    simplification = c("full", "Fungi", "Fungi+rep"),
    external_scoring = c("aside", "matrix", "aside"),
    kingdom_only = c(FALSE, FALSE, TRUE)
  )
  out <- tibble::tibble(
    assignment = assignments,
    db = assignment_db(assignments)
  ) |>
    dplyr::left_join(
      dplyr::select(db_meta, "db", "simplification"),
      by = "db"
    ) |>
    dplyr::left_join(rules, by = "simplification")

  unknown <- out$assignment[!is.na(out$db) & is.na(out$simplification)]
  if (length(unknown) > 0) {
    stop(
      "scoring_groups(): database not in db_meta for ",
      paste(utils::head(unknown, 3), collapse = ", ")
    )
  }

  consensus <- out[is.na(out$db), ]
  out <- out[!is.na(out$db), ]
  if (nrow(consensus) > 0) {
    consensus <- dplyr::bind_rows(lapply(consensus_scorings, function(scoring) {
      consensus$external_scoring <- scoring
      consensus$kingdom_only <- scoring == "aside"
      consensus
    }))
  }
  dplyr::bind_rows(out, consensus)
}

# Score every column of `ranks_df` against `truth`, one comparpq call per
# scoring group. `ranks_df` has one row per rank of `truth_ranks` and one
# column per assignment, as tc_metrics_unit() expects; a consensus column is
# scored once per element of `consensus_scorings`, so the result can hold two
# rows per (consensus column, rank, metric), told apart by `external_scoring`.
score_assignments <- function(
  physeq,
  ranks_df,
  truth,
  db_meta,
  external_truth = NULL,
  truth_ranks = c(
    "Kingdom",
    "Phylum",
    "Class",
    "Order",
    "Family",
    "Genus",
    "Species"
  ),
  consensus_scorings = c("matrix", "aside"),
  kingdom_rank = "Kingdom",
  ...
) {
  groups <- scoring_groups(
    names(ranks_df),
    db_meta,
    consensus_scorings = consensus_scorings
  )
  keys <- unique(groups[, c("external_scoring", "kingdom_only")])
  res <- lapply(seq_len(nrow(keys)), function(i) {
    cols <- groups$assignment[
      groups$external_scoring == keys$external_scoring[[i]] &
        groups$kingdom_only == keys$kingdom_only[[i]]
    ]
    out <- comparpq::tc_metrics_unit(
      physeq,
      ranks_df = ranks_df[, cols, drop = FALSE],
      truth = truth,
      truth_ranks = truth_ranks,
      external_truth = external_truth,
      external_correct_ranks = if (keys$kingdom_only[[i]]) {
        kingdom_rank
      } else {
        NULL
      },
      external_scoring = keys$external_scoring[[i]],
      kingdom_rank = kingdom_rank,
      ...
    )
    out$external_scoring <- keys$external_scoring[[i]]
    out
  })
  dplyr::bind_rows(res)
}
