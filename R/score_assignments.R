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
# Requires comparpq (tc_metrics_unit), dplyr and tidyr.

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
  # "Fungi+rep+cut" (final zoom, OD 40) holds the same representatives as
  # "Fungi+rep", so its controls are read the same way.
  rules <- tibble::tibble(
    simplification = c("full", "Fungi", "Fungi+rep", "Fungi+rep+cut"),
    external_scoring = c("aside", "matrix", "aside", "aside"),
    kingdom_only = c(FALSE, FALSE, TRUE, TRUE)
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
  dplyr::bind_rows(res) |>
    real_units_headline()
}

# The cells (TP, FP, FN) of the real units at one rank, for the paired
# bootstrap intervals of OD 50 (R/paired_bootstrap.R): one row per (column of
# `ranks_df`, scored real unit), under one truth rule (`no_truth`). The
# controls are dropped (OD 49: the headline leaves them out), and so are the
# units left out of the matrix at this rank.
real_unit_cells <- function(
  physeq,
  ranks_df,
  truth,
  truth_ranks,
  rank = "Genus",
  no_truth = "aside"
) {
  comparpq::tc_cells_unit(
    physeq,
    ranks_df = ranks_df,
    truth = truth,
    truth_ranks = truth_ranks,
    ranks = rank,
    no_truth = no_truth
  ) |>
    dplyr::filter(!.data$control, !is.na(.data$cell)) |>
    dplyr::select("method_db", "taxon", "cell")
}

# OD 49 (C22-P8): the headline F1, PPV and FDR are scored on the real units
# only, the controls reporting through their own rates (NA_fake, ext_fungi,
# ext_correct, ctrl_assigned_*), so the scores do not depend on how many
# controls were added. The controls enter the matrix as FP or TN only, so the
# real-units values follow from the cells tc_metrics_unit() returns:
# FP_real = misassign_seq x n_real, checked against FP minus the named
# controls of the matrix. The values with the controls in the matrix stay as
# F1_with_ctrl, PPV_with_ctrl and FDR_with_ctrl. TPR never involves the
# controls; MCC, ACC and TNR keep them, as they need the TN cell that only
# the controls fill.
real_units_headline <- function(res) {
  ids <- c("method_db", "tax_level", "external_scoring")
  headline <- c(
    F1_score = "F1_with_ctrl",
    PPV = "PPV_with_ctrl",
    FDR = "FDR_with_ctrl"
  )
  cells <- res |>
    dplyr::filter(
      .data$metrics %in%
        c(
          "TP",
          "FP",
          "FN",
          "misassign_seq",
          "n_real",
          "ctrl_assigned_fake",
          "ctrl_assigned_ext"
        )
    ) |>
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(ids),
      names_from = "metrics",
      values_from = "values"
    ) |>
    dplyr::mutate(
      FP_real = ifelse(
        .data$n_real > 0,
        round(.data$misassign_seq * .data$n_real),
        0
      ),
      FP_ctrl = .data$ctrl_assigned_fake +
        ifelse(.data$external_scoring == "matrix", .data$ctrl_assigned_ext, 0)
    )
  if (!isTRUE(all(cells$FP_real + cells$FP_ctrl == cells$FP))) {
    stop(
      "real_units_headline(): FP is not the real misassignments plus the ",
      "named controls; has comparpq::tc_metrics_unit() changed its cells?"
    )
  }
  real <- cells |>
    dplyr::transmute(
      dplyr::across(dplyr::all_of(ids)),
      F1_score = 2 * .data$TP / (2 * .data$TP + .data$FP_real + .data$FN),
      PPV = .data$TP / (.data$TP + .data$FP_real),
      FDR = .data$FP_real / (.data$FP_real + .data$TP)
    ) |>
    tidyr::pivot_longer(
      dplyr::all_of(names(headline)),
      names_to = "metrics",
      values_to = "real_value"
    )
  with_ctrl <- res |>
    dplyr::filter(.data$metrics %in% names(headline)) |>
    dplyr::mutate(metrics = unname(headline[.data$metrics]))
  res |>
    dplyr::left_join(real, by = c(ids, "metrics")) |>
    dplyr::mutate(
      values = ifelse(
        .data$metrics %in% names(headline),
        .data$real_value,
        .data$values
      )
    ) |>
    dplyr::select(-"real_value") |>
    dplyr::bind_rows(with_ctrl)
}
