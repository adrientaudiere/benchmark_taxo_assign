# Full analysis pipeline — sourcing this file (re)builds EVERYTHING in order:
# reference databases, then every dataset of the paper (a day or more; at
# night, run_bio_night.R does the same per dataset, with a log).
#
# Pipelines are named projects in _targets.yaml (assign_taxo_mini for smoke
# tests on the mini_* databases with its own store). run_project()
# (R/run_project.R) selects one, builds it and prunes the objects that no
# longer belong to it.
#
# To run a single project, do not source this file; instead:
#   source("R/run_project.R")
#   run_project("assign_taxo_mini"); run_project("cross_val_mini")

source("R/run_project.R")

# Step 0: reference databases (docs/reference_databases.md). Downloads the
#   general FASTA releases of config.R::reference_sources, then derives the
#   dada2 / sintax formats and every database of config.R::benchmark_dbs.
#   Idempotent: skips files that already exist. Pass force = TRUE to rebuild.
source("make_databases.R")
download_reference_sources()
derive_all_variants()

# Step 1: per dataset, its DADA2 (or published-OTU) store → d_asv, then its
# assignment grid; the species section (OD 43) last. Cross-validation is
# parked (OD 28) and not run here.
for (project in unlist(lapply(production_datasets(), dataset_projects))) {
  run_project(project)
}
