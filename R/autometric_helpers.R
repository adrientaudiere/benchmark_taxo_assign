# Per-target resource logging with {autometric} that also works inside crew
# workers.
#
# Why: `autometric::log_start()` only logs the process it is called in. The
# pipelines run each assignment target in a separate crew worker process, so a
# logger started in the main `_targets` process never sees
# `log_phase_set(full_name)` calls made in the workers (ROADMAP S1.1). The
# helpers below start a logger *inside* the target, write one file per target
# and per run, and aggregate the newest *run* of each phase.
#
# Layout: <dir>/<phase>__<YYYYmmdd_HHMMSS>.txt. The timestamp in the file name
# is what `summarise_autometric_costs()` uses to pick the latest run; the phase
# part of the name is documentation only, because a file can hold samples of a
# phase other than the one it is named after (S1.6).
#
# autometric samples the R process only, so it misses the external programs
# that sintax, lca (vsearch) and blastn (BLAST+) run as children of the worker:
# their peak memory read as the idle R session (~0.9 GB) and their CPU as ~0
# (code review §1.3). Two files beside each log fill that gap, read by
# add_tree_costs():
#   <phase>__<stamp>__tree.csv  memory of the worker and all its descendants,
#                               sampled from another process (R/sample_tree.R);
#   <phase>__<stamp>__cpu.csv   CPU seconds of the worker and of the children
#                               it waited for (proc.time()), and wall seconds.

# Run `expr` while logging this process' CPU / memory every `seconds` seconds
# under the phase name `phase`, and the memory and CPU of its process tree.
# Returns the value of `expr`.
with_autometric <- function(phase, expr, dir, seconds = 1) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(
    dir,
    paste0(phase, "__", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt")
  )
  # S1.6a: a logger leaked by a previous target of this worker (interrupted
  # before the on.exit() below) would make log_start() a silent no-op; a stop
  # with nothing running does nothing.
  autometric::log_stop()
  autometric::log_start(path = path, seconds = seconds)
  on.exit(autometric::log_stop(), add = TRUE)
  autometric::log_phase_set(phase)
  sampler <- processx::process$new(
    file.path(R.home("bin"), "Rscript"),
    c(
      here::here("R/sample_tree.R"),
      Sys.getpid(),
      sub("\\.txt$", "__tree.csv", path),
      seconds
    )
  )
  on.exit(sampler$kill(), add = TRUE)
  cpu0 <- proc.time()
  value <- force(expr)
  cpu <- proc.time() - cpu0
  utils::write.csv(
    data.frame(
      phase = phase,
      cpu_self_s = cpu[["user.self"]] + cpu[["sys.self"]],
      cpu_child_s = cpu[["user.child"]] + cpu[["sys.child"]],
      wall_s = cpu[["elapsed"]]
    ),
    sub("\\.txt$", "__cpu.csv", path),
    row.names = FALSE
  )
  value
}

# Read every log file of `dir` into one tibble (autometric defaults: time in
# seconds, memory in megabytes). Adds `log_file` so runs can be told apart.
# Empty files are skipped: a target that ends before the first sampling
# interval (e.g. sintax on a mini_* database, ~1-3 s) writes no row, and
# autometric::log_read() errors on an empty file.
read_autometric_dir <- function(dir) {
  files <- list.files(dir, pattern = "\\.txt$", full.names = TRUE)
  files <- files[file.size(files) > 0]
  if (length(files) == 0) {
    stop("No non-empty autometric log file found in ", dir)
  }
  rows <- lapply(files, function(f) {
    df <- autometric::log_read(f)
    df$log_file <- basename(f)
    df$run_started <- run_started_of(basename(f), f)
    # `autometric::log_read()` returns `time` in seconds since the *first
    # sample of that file*, not an epoch, so times from two files are not
    # comparable as they stand. `abs_time` puts every sample on one clock
    # (epoch seconds, +/- 1 s: the file name carries whole seconds).
    df$abs_time <- as.numeric(df$run_started) + df$time
    df
  })
  dplyr::bind_rows(rows)
}

# Start time of the run a log file belongs to, parsed from the `__YYYYmmdd_HHMMSS`
# suffix written by with_autometric(). Falls back to the file mtime for a name
# that does not carry one, so a hand-renamed file never silently sorts first.
run_started_of <- function(base, path) {
  stamp <- sub("^.*__(\\d{8}_\\d{6})\\.txt$", "\\1", base)
  if (identical(stamp, base)) {
    return(file.mtime(path))
  }
  parsed <- as.POSIXct(stamp, format = "%Y%m%d_%H%M%S", tz = "")
  if (is.na(parsed)) {
    file.mtime(path)
  } else {
    parsed
  }
}

# One row per phase, computed on the newest run of that phase.
#
# Attribution is by the `phase` column, never by the file name (S1.6). A crew
# worker interrupted before `with_autometric()`'s on.exit() leaves its logger
# running, and `autometric::log_phase_set()` is process-global, so the *next*
# target's samples are written into the *previous* target's file, correctly
# labelled with the next target's phase. The file name is therefore the only
# thing that lies; the phase column is right in every sample.
#
# The former implementation kept `log_file == max(log_file)`, an alphabetical
# sort over names that begin with the phase, so a phase whose rows had leaked
# into another phase's file could select by spelling instead of by time — and
# could return a stale run as the current one. The run is now chosen by
# `run_started` (parsed from the file name), and the rows of that run are taken
# from every file whose samples overlap it in time, so a leaked copy is
# reunited with the phase's own file instead of replacing it.
#
# `wall_time_s` is the span between the first and last sample of the run, so
# its resolution is the `seconds` interval of with_autometric() (1 s): fine for
# assignments that take minutes to hours, meaningless for sub-second targets.
# `n_samples` counts samples, not instants: when two loggers of the same worker
# sampled the phase concurrently they interleave at different sub-second
# offsets, so the count is inflated and `n_log_files` says by how many files.
summarise_autometric_costs <- function(log_df) {
  rows <- dplyr::filter(
    log_df,
    !is.na(phase),
    phase != "",
    phase != "__DEFAULT__"
  )
  if (nrow(rows) == 0) {
    return(rows[0, ])
  }
  phases <- sort(unique(rows$phase))
  out <- lapply(phases, function(p) {
    this <- rows[rows$phase == p, ]
    newest <- this[this$run_started == max(this$run_started), ]
    window <- range(newest$abs_time)
    # Any other file holding this phase in the same execution overlaps the
    # newest file on the shared clock; a genuinely earlier run of the phase
    # does not, and is left out.
    overlaps <- vapply(
      split(this, this$log_file),
      \(df) max(df$abs_time) >= window[1] && min(df$abs_time) <= window[2],
      logical(1)
    )
    run <- this[this$log_file %in% names(overlaps)[overlaps], ]
    tibble::tibble(
      phase = p,
      wall_time_s = as.numeric(max(run$abs_time) - min(run$abs_time)),
      peak_resident_mb = max(run$resident, na.rm = TRUE),
      mean_cpu_pct = mean(run$cpu, na.rm = TRUE),
      n_samples = nrow(run),
      n_log_files = length(unique(run$log_file)),
      log_file = paste(sort(unique(run$log_file)), collapse = " + ")
    )
  })
  dplyr::bind_rows(out)
}

# Adds to the rows of summarise_autometric_costs() the cost of each phase's
# process tree, from the newest `__cpu.csv` of that phase in `dir` and the
# `__tree.csv` of the same run (written by with_autometric()):
#   wall_s        elapsed seconds around the computation (proc.time()), valid
#                 also below the 1 s sampling interval;
#   cpu_self_s    CPU seconds of the R worker;
#   cpu_child_s   CPU seconds of the children it waited for (vsearch, BLAST+);
#   peak_tree_mb  highest summed RSS of the worker and its descendants (an
#                 upper bound: shared pages count once per process);
#   n_tree_samples number of memory samples behind that peak.
# Phases logged before these files existed get NA. A phase shorter than the
# 1 s sampling interval has no autometric row but keeps its tree costs: the
# join is full, so its autometric columns are NA and `wall_s` is not.
add_tree_costs <- function(costs, dir) {
  cpu_files <- list.files(dir, pattern = "__cpu\\.csv$", full.names = TRUE)
  empty <- tibble::tibble(
    phase = character(),
    wall_s = numeric(),
    cpu_self_s = numeric(),
    cpu_child_s = numeric(),
    peak_tree_mb = numeric(),
    n_tree_samples = integer()
  )
  tree <- if (length(cpu_files) == 0) {
    empty
  } else {
    cpu <- dplyr::bind_rows(lapply(cpu_files, function(f) {
      df <- utils::read.csv(f)
      df$stamp <- sub("^.*__(\\d{8}_\\d{6})__cpu\\.csv$", "\\1", basename(f))
      df$tree_file <- sub("__cpu\\.csv$", "__tree.csv", f)
      df
    }))
    cpu <- cpu |>
      dplyr::group_by(phase) |>
      dplyr::slice_max(stamp, n = 1, with_ties = FALSE) |>
      dplyr::ungroup()
    peaks <- lapply(cpu$tree_file, function(f) {
      samples <- if (file.exists(f)) utils::read.csv(f) else NULL
      if (is.null(samples) || nrow(samples) == 0) {
        return(c(NA_real_, 0))
      }
      c(max(samples$rss_mb), nrow(samples))
    })
    cpu$peak_tree_mb <- vapply(peaks, `[`, numeric(1), 1)
    cpu$n_tree_samples <- as.integer(vapply(peaks, `[`, numeric(1), 2))
    dplyr::select(cpu, dplyr::all_of(names(empty)))
  }
  dplyr::full_join(costs, tree, by = "phase")
}
