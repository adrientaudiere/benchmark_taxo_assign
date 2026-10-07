# Build {targets} projects of _targets.yaml. Sourcing this file runs nothing,
# unlike make.R and run_bio_night.R, which both use it:
#   source("R/run_project.R")
#   run_project("assign_taxo_mini")

# Build one project. With `prune = TRUE`, then remove the store objects that no
# longer belong to it (targets of databases or parameters removed from
# config.R), but only when the build left nothing outdated, so a failed run
# never prunes. Returns "ok" or the error message, and the number of targets
# still outdated.
run_project <- function(name, prune = TRUE, reporter = "verbose") {
  withr::with_envvar(c(TAR_PROJECT = name), {
    message(
      "== targets project: ",
      name,
      " (store: ",
      targets::tar_config_get("store"),
      ")"
    )
    status <- tryCatch(
      {
        targets::tar_make(reporter = reporter)
        "ok"
      },
      error = function(e) conditionMessage(e)
    )
    still <- targets::tar_outdated(reporter = "silent")
    if (prune && status == "ok" && length(still) == 0) {
      targets::tar_prune()
    }
    list(status = status, still_outdated = still)
  })
}

# The projects that rebuild one dataset, in order: its DADA2 (or published-OTU)
# store, then its assignment grid. "pauvert" is the production mock (projects
# `dada2` and `assign_taxo`); "species" is the species section (OD 43), which
# reads every assignment store and has no DADA2 step.
dataset_projects <- function(dataset) {
  if (dataset == "species") {
    return("species")
  }
  if (dataset == "pauvert") {
    return(c("dada2", "assign_taxo"))
  }
  paste0(c("dada2_", "assign_taxo_"), dataset)
}

# Every dataset of the paper, in the order of the M2.4 run.
production_datasets <- function() {
  c(
    "pauvert",
    "korhonen_wood",
    "gssp_air",
    "rolling_mock",
    "tedersoo_illumina",
    "insilico",
    "species"
  )
}
