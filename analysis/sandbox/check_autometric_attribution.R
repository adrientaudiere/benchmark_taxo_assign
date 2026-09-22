suppressMessages({
  library(dplyr)
  library(targets)
})
source("R/autometric_helpers.R")

# S1.6: does summarise_autometric_costs() attribute cost to the right phase?
#
# The check is a diff against the implementation that was in place until
# 2026-09-22 (kept below as `old_summarise()`), plus a validation against an
# independent measurement of the same quantity: `targets::tar_meta()$seconds`.
# A log file can hold samples of a phase other than the one it is named after,
# so the old `log_file == max(log_file)` selected a run by spelling; the phase
# column itself has always been right.

dir <- "data/data_final/autometric/assign_taxo"
store <- "store_assign_taxo"

old_summarise <- function(log_df) {
  log_df |>
    filter(!is.na(phase), phase != "", phase != "__DEFAULT__") |>
    group_by(phase) |>
    filter(log_file == max(log_file)) |>
    summarise(
      wall_time_s = as.numeric(max(time) - min(time)),
      peak_resident_mb = max(resident, na.rm = TRUE),
      .groups = "drop"
    )
}

d <- read_autometric_dir(dir)
message(
  "rows: ",
  nrow(d),
  " | files: ",
  dplyr::n_distinct(d$log_file),
  " | phases: ",
  dplyr::n_distinct(d$phase)
)

multi <- d |>
  group_by(log_file) |>
  summarise(n_phase = n_distinct(phase), .groups = "drop") |>
  filter(n_phase > 1)
message("log files holding more than one phase: ", nrow(multi))

o <- old_summarise(d)
n <- summarise_autometric_costs(d)
cmp <- full_join(o, n, by = "phase", suffix = c("_old", "_new")) |>
  mutate(
    d_peak = peak_resident_mb_new - peak_resident_mb_old,
    d_wall = wall_time_s_new - wall_time_s_old
  )
message(
  "phases with a corrected peak: ",
  sum(abs(cmp$d_peak) > 1e-6, na.rm = TRUE),
  " (understated by the old code: ",
  sum(cmp$d_peak > 0, na.rm = TRUE),
  ")"
)
message(
  "phases with a corrected wall_time_s: ",
  sum(abs(cmp$d_wall) > 1e-6, na.rm = TRUE)
)

message("\nworst understatements of the old code (peak, MB)")
cmp |>
  filter(d_peak > 1e-6) |>
  mutate(ratio = peak_resident_mb_new / peak_resident_mb_old) |>
  arrange(desc(ratio)) |>
  head(5) |>
  transmute(
    phase = substr(phase, 1, 46),
    peak_old = round(peak_resident_mb_old, 1),
    peak_new = round(peak_resident_mb_new, 1),
    ratio = round(ratio, 1)
  ) |>
  as.data.frame() |>
  print()

# Independent validation: targets records its own wall time per target.
m <- tar_meta(store = store) |>
  filter(!is.na(seconds), !is.na(time)) |>
  select(phase = name, tar_seconds = seconds)
j <- m |>
  inner_join(select(o, phase, wall_old = wall_time_s), by = "phase") |>
  inner_join(select(n, phase, wall_new = wall_time_s), by = "phase") |>
  mutate(err_old = wall_old - tar_seconds, err_new = wall_new - tar_seconds)

message("\nagainst targets::seconds, by how long the target actually ran")
j |>
  mutate(
    bucket = cut(
      tar_seconds,
      c(-Inf, 1, 10, 60, Inf),
      labels = c("< 1 s", "1-10 s", "10-60 s", "> 60 s")
    )
  ) |>
  group_by(bucket) |>
  summarise(
    n = n(),
    med_err_old = round(median(abs(err_old)), 1),
    med_err_new = round(median(abs(err_new)), 1),
    off10_old = sum(abs(err_old) > 10),
    off10_new = sum(abs(err_new) > 10),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()

# The fix must not lose ground on the targets that carry the figures: every
# computation runs for minutes, and those are the rows benchmark_costs joins.
long <- filter(j, tar_seconds > 60)
stopifnot(sum(abs(long$err_new) > 10) <= sum(abs(long$err_old) > 10))
message("\nlong targets no worse than before: OK")
