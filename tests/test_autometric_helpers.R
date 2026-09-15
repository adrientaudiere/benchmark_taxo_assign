# Per-target autometric logging helpers (R/autometric_helpers.R). Run from
# the project root:
#   Rscript tests/test_autometric_helpers.R

library("testthat")
library("here")

here::i_am("tests/test_autometric_helpers.R")
source(here("R/autometric_helpers.R"))

test_that("summarise_autometric_costs keeps only the newest log file per phase", {
  log_df <- tibble::tibble(
    phase    = c(rep("A", 3), rep("A", 2), rep("B", 2), rep("__DEFAULT__", 2)),
    time     = c(0, 10, 20,   100, 101,   0, 5,   0, 1),
    resident = c(100, 200, 300,  50, 60,   10, 20,  1, 1),
    cpu      = c(50, 50, 50,     80, 80,   10, 10,  0, 0),
    log_file = c(rep("A__20260101_000000.txt", 3), rep("A__20260102_000000.txt", 2),
                 rep("B__20260101_000000.txt", 2), rep("A__20260102_000000.txt", 2))
  )
  out <- summarise_autometric_costs(log_df)
  expect_equal(sort(out$phase), c("A", "B"))
  a <- out[out$phase == "A", ]
  expect_equal(a$log_file, "A__20260102_000000.txt")
  expect_equal(a$n_samples, 2)
  expect_equal(a$wall_time_s, 1)
  expect_equal(a$peak_resident_mb, 60)
  expect_equal(a$mean_cpu_pct, 80)
  b <- out[out$phase == "B", ]
  expect_equal(b$wall_time_s, 5)
})

test_that("read_autometric_dir fails loudly on an empty directory", {
  expect_error(read_autometric_dir(tempfile()), "No non-empty autometric log file")
})

test_that("read_autometric_dir skips empty log files (targets shorter than one sample)", {
  skip_if_not_installed("autometric")
  dir <- file.path(tempfile(), "am")
  with_autometric("long", { Sys.sleep(1.2); 1 }, dir = dir)
  file.create(file.path(dir, "short__20260101_000000.txt"))
  expect_length(list.files(dir), 2)
  logged <- read_autometric_dir(dir)
  expect_equal(unique(logged$phase), "long")
})

test_that("with_autometric returns the value and writes one file per call", {
  skip_if_not_installed("autometric")
  dir <- file.path(tempfile(), "am")
  value <- with_autometric("toy", { Sys.sleep(1.2); 42 }, dir = dir)
  expect_equal(value, 42)
  files <- list.files(dir)
  expect_length(files, 1)
  expect_match(files, "^toy__[0-9]{8}_[0-9]{6}\\.txt$")
  logged <- read_autometric_dir(dir)
  expect_true(all(logged$phase == "toy"))
  expect_equal(unique(logged$log_file), files)
})
