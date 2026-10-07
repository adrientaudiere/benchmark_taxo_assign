### Overnight runs, one dataset after the other, in one log.
###
### For each dataset given on the command line, in order, builds its DADA2 (or
### published-OTU) store, then its assignment grid (R/run_project.R,
### dataset_projects()): "pauvert" is `dada2` + `assign_taxo`, any other name
### `dada2_<dataset>` + `assign_taxo_<dataset>` (korhonen_wood, gssp_air,
### insilico, tedersoo_illumina, rolling_mock). The DADA2 step is quick when
### its store is up to date; when it fails or leaves targets outdated, the
### assignment grid of that dataset is skipped. No pruning: run_project(prune =
### TRUE) by hand once the run is checked.
###
### Since 2026-09-27 every project assigns its ITSx inputs only (OD 34). Since 2026-09-29 the computations run one at a time
### (config.R single_worker = TRUE): with the dada2 worker beside two fast
### workers, dada2 x full EUK_ITS_v2.1 (about 57 GB) got the 2026-09-28 run
### killed by the OOM killer. Run it with the graphical session closed, e.g.
### from a systemd timer or over ssh inside tmux (S8.6).
###
###   Rscript run_bio_night.R pauvert korhonen_wood
###
### "species" runs the species section (pipelines/species.R, OD 43), which
### reads every assignment store: put it last.
###
### A dataset whose run fails is logged and the next one still runs.

source("R/run_project.R")

datasets <- commandArgs(trailingOnly = TRUE)
if (length(datasets) == 0) {
  datasets <- c("pauvert", "korhonen_wood")
}

dir.create("logs", showWarnings = FALSE)
log_path <- file.path(
  "logs",
  paste0("run_bio_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)
log_con <- file(log_path, open = "wt")
sink(log_con, split = TRUE)
sink(log_con, type = "message")

log_line <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ...)
}
log_line("log: ", log_path, " | datasets: ", paste(datasets, collapse = ", "))

for (dataset in datasets) {
  for (project in dataset_projects(dataset)) {
    Sys.setenv(TAR_PROJECT = project)
    log_line(
      "== ",
      project,
      " | outdated targets: ",
      length(targets::tar_outdated(reporter = "silent"))
    )
    run <- run_project(project, prune = FALSE, reporter = "timestamp")
    log_line(
      "== ",
      project,
      " finished (",
      run$status,
      "). still outdated: ",
      length(run$still_outdated)
    )
    if (length(run$still_outdated) > 0) {
      log_line(
        "remaining: ",
        paste(utils::head(run$still_outdated, 20), collapse = ", ")
      )
    }
    # An assignment grid on a DADA2 store that did not finish would read a
    # stale d_asv and be invalidated by the next DADA2 run: skip the rest of
    # this dataset.
    if (run$status != "ok" || length(run$still_outdated) > 0) {
      log_line("skipping the rest of ", dataset)
      break
    }
  }
}

log_line("all datasets done")

# Today's figures in front of ETAT.md (tools/etat.R, a few minutes).
etat_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  "tools/etat.R",
  stdout = FALSE,
  stderr = FALSE
)
log_line("ETAT.md ", if (etat_status == 0) "updated" else "NOT updated")
sink(type = "message")
sink()
close(log_con)
