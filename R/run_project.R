# Build one {targets} project of _targets.yaml, then prune the store objects
# that no longer belong to it (targets of databases or parameters removed from
# config.R). Sourcing this file runs nothing, unlike make.R:
#   source("R/run_project.R")
#   run_project("assign_taxo_mini")
run_project <- function(name) {
  withr::with_envvar(c(TAR_PROJECT = name), {
    message("== targets project: ", name, " (store: ", targets::tar_config_get("store"), ")")
    targets::tar_make()
    targets::tar_prune()
  })
}
