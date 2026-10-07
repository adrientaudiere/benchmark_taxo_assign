# ETAT.md helpers (R/etat_helpers.R). Run from the project root:
#   Rscript tests/test_etat_helpers.R

library("testthat")
library("here")

here::i_am("tests/test_etat_helpers.R")
source(here("R/etat_helpers.R"))

header <- c("# ETAT", "", "Intro.", "")

test_that("prepend_etat puts the new block after the header", {
  old <- c(header, "## 2026-10-01", "", "old", "")
  out <- prepend_etat(old, c("## 2026-10-07", "", "new", ""), "2026-10-07")
  expect_equal(out[1:4], header)
  expect_equal(out[5], "## 2026-10-07")
  expect_equal(
    grep("^## ", out, value = TRUE),
    c("## 2026-10-07", "## 2026-10-01")
  )
})

test_that("prepend_etat replaces a block of the same date", {
  old <- c(
    header,
    "## 2026-10-07",
    "",
    "first",
    "",
    "## 2026-10-01",
    "",
    "old",
    ""
  )
  out <- prepend_etat(old, c("## 2026-10-07", "", "second", ""), "2026-10-07")
  expect_false("first" %in% out)
  expect_true(all(c("second", "old") %in% out))
  expect_equal(sum(out == "## 2026-10-07"), 1)
})

test_that("prepend_etat starts a file without dated blocks", {
  out <- prepend_etat(header, c("## 2026-10-07", ""), "2026-10-07")
  expect_equal(out, c(header, "## 2026-10-07", ""))
})

test_that("test_counts counts files, blocks and expectations", {
  dir <- tempfile()
  dir.create(dir)
  writeLines(
    c("test_that('a', {", "  expect_equal(1, 1)", "  expect_true(TRUE)", "})"),
    file.path(dir, "test_x.R")
  )
  writeLines("x <- 1", file.path(dir, "helper.R"))
  expect_equal(
    test_counts(dir),
    c(files = 1, tests = 1, expectations = 2)
  )
})
