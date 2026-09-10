# Per-target resource logging with {autometric} that also works inside crew
# workers.
#
# Why: `autometric::log_start()` only logs the process it is called in. The
# pipelines run each assignment target in a separate crew worker process, so a
# logger started in the main `_targets` process never sees
# `log_phase_set(full_name)` calls made in the workers (ROADMAP S1.1). The
# helpers below start a logger *inside* the target, write one file per target
# and per run, and aggregate the newest file per phase.
#
# Layout: <dir>/<phase>__<YYYYmmdd_HHMMSS>.txt. The timestamp in the file name
# is what `summarise_autometric_costs()` uses to pick the latest run.

# Run `expr` while logging this process' CPU / memory every `seconds` seconds
# under the phase name `phase`. Returns the value of `expr`.
with_autometric <- function(phase, expr, dir, seconds = 1) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(
    dir,
    paste0(phase, "__", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  )
  autometric::log_start(path = path, seconds = seconds)
  on.exit(autometric::log_stop(), add = TRUE)
  autometric::log_phase_set(phase)
  force(expr)
}

# Read every log file of `dir` into one tibble (autometric defaults: time in
# seconds, memory in megabytes). Adds `log_file` so runs can be told apart.
read_autometric_dir <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No autometric log file found in ", dir)
  }
  rows <- lapply(files, function(f) {
    df <- autometric::log_read(f)
    df$log_file <- basename(f)
    df
  })
  dplyr::bind_rows(rows)
}

# One row per phase, computed on the newest log file of that phase only.
# `wall_time_s` is the span between the first and last sample of the run, so
# its resolution is the `seconds` interval of with_autometric() (1 s): fine for
# assignments that take minutes to hours, meaningless for sub-second targets.
summarise_autometric_costs <- function(log_df) {
  log_df |>
    dplyr::filter(!is.na(phase), phase != "", phase != "__DEFAULT__") |>
    dplyr::group_by(phase) |>
    dplyr::filter(log_file == max(log_file)) |>
    dplyr::summarise(
      wall_time_s      = as.numeric(max(time) - min(time)),
      peak_resident_mb = max(resident, na.rm = TRUE),
      mean_cpu_pct     = mean(cpu, na.rm = TRUE),
      n_samples        = dplyr::n(),
      log_file         = dplyr::first(log_file),
      .groups          = "drop"
    )
}
