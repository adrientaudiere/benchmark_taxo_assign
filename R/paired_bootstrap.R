# Paired bootstrap intervals on the differences of Genus F1 and PPV between
# assignment columns scored on the same units (OD 50, accepts C22-P11).
#
# Input: the per-unit cells of the real units (comparpq::tc_cells_unit(),
# through real_unit_cells() in R/score_assignments.R), one row per
# (column, unit) with the columns
#   column  the assignment column (any id, e.g. "ITS1:dada2__EUK...");
#   unit_id the scored unit, unique across the inputs that are resampled
#           together (e.g. "OTU:Taxa_1" and "ASV:Taxa_1" are two units);
#   group   the resampling unit: the unit itself, or the OTU cluster it
#           belongs to (Q5), so that the units of a cluster are drawn together;
#   cell    "TP", "FP" or "FN".
# A draw resamples the groups with replacement; every column is scored on the
# same draw, so a difference between columns is paired. The intervals are
# percentile intervals over `n_boot` draws.
#
# Requires stats and withr.

# Genus F1 and PPV of every column on each of `n_boot` draws, and on the
# observed units. Returns a list: `F1` and `PPV` (n_boot x columns matrices)
# and `point` (a list of two named vectors, one value per column).
boot_scores <- function(cells, n_boot = 1e4, seed = 1) {
  if (!all(c("column", "unit_id", "group", "cell") %in% names(cells))) {
    stop("boot_scores(): `cells` needs column, unit_id, group and cell.")
  }
  if (anyNA(cells$cell) || !all(cells$cell %in% c("TP", "FP", "FN"))) {
    stop("boot_scores(): `cell` must be TP, FP or FN (real units only).")
  }
  units <- unique(cells[, c("unit_id", "group")])
  if (anyDuplicated(units$unit_id)) {
    stop("boot_scores(): a unit belongs to two groups.")
  }
  columns <- unique(cells$column)
  groups <- unique(units$group)

  # One indicator matrix per cell: units x columns, 0 where a unit is not
  # scored by a column (an input it does not belong to).
  indicator <- function(which_cell) {
    m <- matrix(
      0,
      nrow = nrow(units),
      ncol = length(columns),
      dimnames = list(units$unit_id, columns)
    )
    hit <- cells[cells$cell == which_cell, ]
    m[cbind(match(hit$unit_id, units$unit_id), match(hit$column, columns))] <- 1
    m
  }
  tp <- indicator("TP")
  fp <- indicator("FP")
  fn <- indicator("FN")

  scores <- function(TP, FP, FN) {
    list(F1 = 2 * TP / (2 * TP + FP + FN), PPV = TP / (TP + FP))
  }
  point <- scores(colSums(tp), colSums(fp), colSums(fn))

  # Weight of each group on each draw (groups x n_boot), then of each unit.
  weights <- withr::with_seed(seed, {
    stats::rmultinom(
      n_boot,
      length(groups),
      rep(1 / length(groups), length(groups))
    )
  })
  unit_weights <- weights[match(units$group, groups), , drop = FALSE]
  draws <- scores(
    crossprod(unit_weights, tp),
    crossprod(unit_weights, fp),
    crossprod(unit_weights, fn)
  )
  list(F1 = draws$F1, PPV = draws$PPV, point = point)
}

# The paired difference mean(columns `a`) - mean(columns `b`) of `metric`
# ("F1" or "PPV"): its observed value and its 95 % percentile interval over
# the draws of boot_scores(). A draw where a PPV is undefined (no unit named)
# is dropped.
boot_contrast <- function(boot, a, b, metric = c("F1", "PPV"), level = 0.95) {
  metric <- match.arg(metric)
  if (!all(c(a, b) %in% colnames(boot[[metric]]))) {
    stop(
      "boot_contrast(): unknown column ",
      setdiff(c(a, b), colnames(boot[[metric]]))[1]
    )
  }
  draws <- rowMeans(boot[[metric]][, a, drop = FALSE]) -
    rowMeans(boot[[metric]][, b, drop = FALSE])
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  bounds <- stats::quantile(draws, probs, na.rm = TRUE, names = FALSE)
  data.frame(
    metric = metric,
    delta = mean(boot$point[[metric]][a]) - mean(boot$point[[metric]][b]),
    lo = bounds[[1]],
    hi = bounds[[2]]
  )
}

# The paired difference max(columns `a`) - max(columns `b`) of `metric`,
# with the best column of each set chosen again on every draw, so that the
# interval carries the selection of the two winners (e.g. the best consensus
# strategy against the best single column, Q3): an interval on two columns
# fixed after choosing them on the observed units would be too narrow, and
# the more so the more candidates `b` holds.
boot_contrast_best <- function(
  boot,
  a,
  b,
  metric = c("F1", "PPV"),
  level = 0.95
) {
  metric <- match.arg(metric)
  if (!all(c(a, b) %in% colnames(boot[[metric]]))) {
    stop(
      "boot_contrast_best(): unknown column ",
      setdiff(c(a, b), colnames(boot[[metric]]))[1]
    )
  }
  row_max <- function(m) {
    apply(m, 1, \(x) suppressWarnings(max(x, na.rm = TRUE)))
  }
  draws <- row_max(boot[[metric]][, a, drop = FALSE]) -
    row_max(boot[[metric]][, b, drop = FALSE])
  draws[!is.finite(draws)] <- NA
  probs <- c((1 - level) / 2, 1 - (1 - level) / 2)
  bounds <- stats::quantile(draws, probs, na.rm = TRUE, names = FALSE)
  data.frame(
    metric = metric,
    delta = max(boot$point[[metric]][a], na.rm = TRUE) -
      max(boot$point[[metric]][b], na.rm = TRUE),
    lo = bounds[[1]],
    hi = bounds[[2]]
  )
}
