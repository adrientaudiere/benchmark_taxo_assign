# Per-target autometric logging helpers (R/autometric_helpers.R). Run from
# the project root:
#   Rscript tests/test_autometric_helpers.R

library("testthat")
library("here")

here::i_am("tests/test_autometric_helpers.R")
source(here("R/autometric_helpers.R"))

# The columns read_autometric_dir() adds: `time` is relative to the first
# sample of its file, so `abs_time` puts every file on one clock.
with_clock <- function(log_df) {
  log_df$run_started <- do.call(
    c,
    lapply(log_df$log_file, \(f) run_started_of(f, f))
  )
  log_df$abs_time <- as.numeric(log_df$run_started) + log_df$time
  log_df
}

test_that("summarise_autometric_costs keeps only the newest run per phase", {
  log_df <- with_clock(tibble::tibble(
    phase = c(rep("A", 3), rep("A", 2), rep("B", 2), rep("__DEFAULT__", 2)),
    time = c(0, 10, 20, 0, 1, 0, 5, 0, 1),
    resident = c(100, 200, 300, 50, 60, 10, 20, 1, 1),
    cpu = c(50, 50, 50, 80, 80, 10, 10, 0, 0),
    log_file = c(
      rep("A__20260101_000000.txt", 3),
      rep("A__20260102_000000.txt", 2),
      rep("B__20260101_000000.txt", 2),
      rep("A__20260102_000000.txt", 2)
    )
  ))
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

# S1.6: a logger leaked by target A keeps writing after B starts, so part of
# B's samples sit in A's file, labelled B. They must be reunited with B's own
# file, and an older run of B must stay out; the former
# `log_file == max(log_file)` picked a file by spelling instead.
test_that("summarise_autometric_costs reunites a leaked copy and drops an old run", {
  log_df <- with_clock(tibble::tibble(
    phase = c("A", "A", "B", "B", "B", "B", "B", "B"),
    time = c(0, 60, 120, 400, 0, 280, 0, 5),
    resident = c(1000, 1100, 9000, 3000, 2000, 2500, 99000, 99000),
    cpu = 50,
    log_file = c(
      rep("A__20260922_120000.txt", 4),
      rep("B__20260922_120200.txt", 2),
      rep("B__20260901_080000.txt", 2)
    )
  ))
  b <- summarise_autometric_costs(log_df)
  b <- b[b$phase == "B", ]
  expect_equal(b$n_log_files, 2)
  expect_equal(b$peak_resident_mb, 9000)
  expect_equal(b$wall_time_s, 280)
  expect_false(grepl("20260901", b$log_file))
})

test_that("read_autometric_dir fails loudly on an empty directory", {
  expect_error(
    read_autometric_dir(tempfile()),
    "No non-empty autometric log file"
  )
})

test_that("read_autometric_dir skips empty log files (targets shorter than one sample)", {
  skip_if_not_installed("autometric")
  dir <- file.path(tempfile(), "am")
  with_autometric(
    "long",
    {
      Sys.sleep(1.2)
      1
    },
    dir = dir
  )
  file.create(file.path(dir, "short__20260101_000000.txt"))
  expect_length(list.files(dir, pattern = "\\.txt$"), 2)
  logged <- read_autometric_dir(dir)
  expect_equal(unique(logged$phase), "long")
})

test_that("with_autometric returns the value and writes one file per call", {
  skip_if_not_installed("autometric")
  dir <- file.path(tempfile(), "am")
  value <- with_autometric(
    "toy",
    {
      Sys.sleep(1.2)
      42
    },
    dir = dir
  )
  expect_equal(value, 42)
  files <- list.files(dir, pattern = "\\.txt$")
  expect_length(files, 1)
  expect_match(files, "^toy__[0-9]{8}_[0-9]{6}\\.txt$")
  stem <- sub("\\.txt$", "", files)
  expect_true(file.exists(file.path(dir, paste0(stem, "__cpu.csv"))))
  expect_true(file.exists(file.path(dir, paste0(stem, "__tree.csv"))))
  logged <- read_autometric_dir(dir)
  expect_true(all(logged$phase == "toy"))
  expect_equal(unique(logged$log_file), files)
})

test_that("with_autometric measures the memory and CPU of child processes", {
  skip_if_not_installed("autometric")
  dir <- file.path(tempfile(), "am")
  child <- c(
    "-e",
    shQuote(paste(
      "x <- numeric(5e7); x[] <- 1;",
      "t0 <- proc.time()[[3]];",
      "while (proc.time()[[3]] - t0 < 2) sum(sqrt(seq_len(1e5)))"
    ))
  )
  with_autometric(
    "child",
    system2(file.path(R.home("bin"), "Rscript"), child),
    dir = dir
  )
  costs <- add_tree_costs(tibble::tibble(phase = "child"), dir)
  # the child's 400 MB vector and ~2 s of CPU, invisible to autometric
  expect_gt(costs$peak_tree_mb, 300)
  expect_gt(costs$cpu_child_s, 1)
  expect_gte(costs$wall_s, 2)
})

test_that("add_tree_costs keeps the newest run of each phase", {
  dir <- tempfile()
  dir.create(dir)
  write_run <- function(phase, stamp, wall, rss) {
    stem <- file.path(dir, paste0(phase, "__", stamp))
    write.csv(
      data.frame(phase = phase, cpu_self_s = 1, cpu_child_s = 2, wall_s = wall),
      paste0(stem, "__cpu.csv"),
      row.names = FALSE
    )
    write.csv(
      data.frame(epoch = c(1, 2), rss_mb = rss, n_procs = 2L),
      paste0(stem, "__tree.csv"),
      row.names = FALSE
    )
  }
  write_run("p", "20261001_100000", wall = 10, rss = c(100, 200))
  write_run("p", "20261007_100000", wall = 20, rss = c(300, 500))
  costs <- add_tree_costs(tibble::tibble(phase = c("p", "q")), dir)
  expect_equal(costs$wall_s, c(20, NA))
  expect_equal(costs$peak_tree_mb, c(500, NA))
  expect_equal(costs$n_tree_samples, c(2L, NA))
})

test_that("add_tree_costs keeps a phase too short for autometric", {
  dir <- tempfile()
  dir.create(dir)
  stem <- file.path(dir, "short__20261007_100000")
  write.csv(
    data.frame(
      phase = "short",
      cpu_self_s = 0.1,
      cpu_child_s = 0.6,
      wall_s = 0.8
    ),
    paste0(stem, "__cpu.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(epoch = 1, rss_mb = 950, n_procs = 3L),
    paste0(stem, "__tree.csv"),
    row.names = FALSE
  )
  costs <- add_tree_costs(
    tibble::tibble(phase = "long", wall_time_s = 60),
    dir
  )
  short <- costs[costs$phase == "short", ]
  expect_equal(nrow(costs), 2)
  expect_true(is.na(short$wall_time_s))
  expect_equal(short$wall_s, 0.8)
  expect_equal(short$peak_tree_mb, 950)
})
