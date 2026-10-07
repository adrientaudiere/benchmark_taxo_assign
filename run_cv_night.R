### Overnight cross-validation run, staged by memory (ROADMAP S8.6).
###
### Stage 1 builds every outdated cross_val target EXCEPT the six dada2 cells
### on the EUKARYOME databases. Those are the ones that OOM-killed the run of
### 2026-09-14: `cv_reduce_reference = FALSE` keeps the whole reference in the
### worker, and dada2 on EUK_ITS_v2.1 (1 877 003 records) was measured at
### >= 45 GB before the kill, then at 58.1 GB on the assignment side once the
### S1.6 attribution bug was fixed. With `cv_n_workers_fast = 2` the machine
### would be asked for 58.1 + 6.2 + 5.9 = 70.2 GB against 62 GB installed.
###
### Stage 1 peak, from docs/diagnostics_cv_and_databases.md section 7:
###   dada2 x Unite_s_all (23.5) + lca (6.2) + sintax (5.9) = 35.6 GB.
###
### Stage 2 (the six heavy cells) runs only when CV_HEAVY=TRUE, and only makes
### sense with the desktop closed: it needs ~58 GB to itself. Running it in the
### same tar_make() as stage 1 would let targets schedule a heavy dada2 cell
### beside two fast workers, which is the configuration that was killed.

### Deliberately does not use R/run_project.R: that helper calls tar_prune()
### after tar_make(), which has no business running between two stages of a
### partial build. Prune by hand once the CV is complete.
Sys.setenv(TAR_PROJECT = "cross_val")

### The script keeps its own log, so an overnight run is readable whatever it
### was launched from and needs no shell redirection.
log_path <- paste0(
  "run_cross_val_",
  format(Sys.time(), "%Y%m%d_%H%M%S"),
  ".log"
)
### Closed at the end of the script, not with on.exit(): on.exit() is a
### function-frame hook and does nothing at top level. R flushes and closes the
### connection when the process exits, so the log survives an error too.
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

heavy_pattern <- "^cv__dada2__EUK_ITS_v2\\.1"
log_line <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}
log_line("log: ", log_path)

outdated <- targets::tar_outdated(reporter = "silent")
heavy <- grep(heavy_pattern, outdated, value = TRUE)
### Aggregators depend on EVERY cv__ target, so naming one in stage 1 makes
### targets build the heavy cells as its dependencies — which is exactly what
### happened on the 2026-09-22 attempt: stage 1 dispatched
### `cv__dada2__EUK_ITS_v2.1__standard` seconds after starting. Stage 1 is
### therefore restricted to individual cv__ cells; the aggregators are built
### once the heavy cells exist.
aggregators <- setdiff(outdated, grep("^cv__", outdated, value = TRUE))
light <- setdiff(outdated, c(heavy, aggregators))

log_line(
  "outdated: ",
  length(outdated),
  " | stage 1: ",
  length(light),
  " | deferred (heavy): ",
  length(heavy),
  " | deferred (aggregators): ",
  length(aggregators)
)
if (length(heavy) > 0) {
  log_line("deferred cells: ", paste(heavy, collapse = ", "))
}

### tar_make() evaluates `names` inside its callr subprocess, which cannot see
### a variable of this session: `names = all_of(light)` fails with "objet
### 'light' introuvable". The selection is therefore passed through a file whose
### path travels in an environment variable, and read on the other side.
make_these <- function(selection, label) {
  if (length(selection) == 0) {
    log_line(label, ": nothing to build")
    return(invisible(NULL))
  }
  path <- file.path(tempdir(), paste0("cv_sel_", label, ".rds"))
  saveRDS(selection, path)
  Sys.setenv(CV_SEL_FILE = path)
  log_line(label, " starting: ", length(selection), " targets")
  targets::tar_make(
    names = tidyselect::all_of(readRDS(Sys.getenv("CV_SEL_FILE"))),
    reporter = "timestamp"
  )
  log_line(label, " done")
}

make_these(light, "stage 1")

run_heavy <- toupper(Sys.getenv("CV_HEAVY", "FALSE")) %in% c("TRUE", "1", "YES")
if (run_heavy && length(heavy) > 0) {
  log_line("stage 2: heavy dada2 x EUKARYOME cells, ~58 GB each, one at a time")
  make_these(heavy, "stage 2")
  ### Only now can the aggregators run: they need every cell to exist.
  make_these(aggregators, "stage 3 (aggregators)")
} else if (length(heavy) > 0) {
  log_line(
    "stage 2 skipped. To run it with the machine otherwise idle: CV_HEAVY=TRUE Rscript run_cv_night.R"
  )
}

still <- targets::tar_outdated(reporter = "silent")
log_line("finished. still outdated: ", length(still))
if (length(still) > 0) {
  log_line("remaining: ", paste(still, collapse = ", "))
}

sink(type = "message")
sink()
close(log_con)
