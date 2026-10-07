### Writes today's figures in front of ETAT.md (R/etat_helpers.R). Run from
### the project root; run_bio_night.R calls it at the end of every night.
###
###   Rscript tools/etat.R
###
### Takes a few minutes: tar_outdated() sources every pipeline.

library("here")
here::i_am("tools/etat.R")
source(here("R/etat_helpers.R"))

projects_yaml <- yaml::read_yaml(here("_targets.yaml"))
date <- format(Sys.Date())

projects <- lapply(names(projects_yaml), function(p) {
  message("etat: ", p)
  project_counts(p, here(projects_yaml[[p]]$store))
}) |>
  dplyr::bind_rows()

units <- lapply(names(projects_yaml), function(p) {
  input_units(p, here(projects_yaml[[p]]$store))
}) |>
  dplyr::bind_rows()

block <- etat_block(date, projects, units, test_counts(here("tests")))

path <- here("ETAT.md")
lines <- if (file.exists(path)) {
  readLines(path)
} else {
  c(
    "# ETAT",
    "",
    "Dated figures of the project, newest first. Written by `tools/etat.R`",
    "(`R/etat_helpers.R`), never by hand: README, CONTEXT, AGENTS and",
    "ROADMAP state the decisions and point here for the counts.",
    ""
  )
}
writeLines(prepend_etat(lines, block, date), path)
message("etat: ", path, " updated for ", date)
