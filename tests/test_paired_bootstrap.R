# Paired bootstrap intervals (R/paired_bootstrap.R, OD 50). Run from the
# project root:
#   Rscript tests/test_paired_bootstrap.R

library("testthat")
library("here")

here::i_am("tests/test_paired_bootstrap.R")
source(here("R/paired_bootstrap.R"))

# Two columns on four units: `a` right on three, `b` on two.
boot_cells <- function(group = c("u1", "u2", "u3", "u4")) {
  data.frame(
    column = rep(c("a", "b"), each = 4),
    unit_id = rep(c("u1", "u2", "u3", "u4"), 2),
    group = rep(group, 2),
    cell = c("TP", "TP", "TP", "FN", "TP", "TP", "FP", "FN")
  )
}

test_that("the observed scores are the F1 and PPV of the cells", {
  boot <- boot_scores(boot_cells(), n_boot = 50, seed = 1)
  expect_equal(unname(boot$point$F1[["a"]]), 2 * 3 / (2 * 3 + 0 + 1))
  expect_equal(unname(boot$point$PPV[["b"]]), 2 / 3)
  expect_equal(dim(boot$F1), c(50, 2))
  res <- boot_contrast(boot, "a", "b")
  expect_equal(res$delta, boot$point$F1[["a"]] - boot$point$F1[["b"]])
  expect_true(res$lo <= res$delta && res$delta <= res$hi)
})

test_that("a column against itself gives an interval of zero", {
  boot <- boot_scores(boot_cells(), n_boot = 50, seed = 1)
  res <- boot_contrast(boot, "a", "a", metric = "PPV")
  expect_equal(c(res$delta, res$lo, res$hi), c(0, 0, 0))
})

test_that("the units of a group are drawn together", {
  cells <- data.frame(
    column = "a",
    unit_id = c("u1", "u2", "u3", "u4"),
    group = c("g1", "g1", "g2", "g2"),
    cell = c("TP", "FN", "TP", "FN")
  )
  # Each group brings one TP and one FN: every draw scores F1 = 2/3.
  boot <- boot_scores(cells, n_boot = 50, seed = 1)
  expect_equal(unique(round(as.vector(boot$F1), 10)), round(2 / 3, 10))
  # Drawn unit by unit, the draws vary.
  cells$group <- cells$unit_id
  expect_gt(length(unique(boot_scores(cells, n_boot = 50, seed = 1)$F1)), 1)
})

test_that("the same seed gives the same draws", {
  expect_identical(
    boot_scores(boot_cells(), n_boot = 20, seed = 3)$F1,
    boot_scores(boot_cells(), n_boot = 20, seed = 3)$F1
  )
})

test_that("cells other than TP, FP and FN are refused", {
  cells <- boot_cells()
  cells$cell[1] <- "TN"
  expect_error(boot_scores(cells, n_boot = 5), "real units only")
  expect_error(
    boot_contrast(boot_scores(boot_cells(), n_boot = 5), "a", "z"),
    "unknown column"
  )
})

test_that("boot_contrast_best chooses the best column again on each draw", {
  boot <- boot_scores(boot_cells(), n_boot = 200, seed = 1)
  res <- boot_contrast_best(boot, "a", c("a", "b"))
  # `a` against the best of {a, b}: never above 0, and 0 on the observed units.
  expect_equal(res$delta, 0)
  expect_lte(res$hi, 0)
  # Against a set of one column it is the fixed contrast.
  expect_equal(
    boot_contrast_best(boot, "a", "b"),
    boot_contrast(boot, "a", "b")
  )
})
