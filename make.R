# Full analysis pipeline — sourcing this file (re)builds EVERYTHING in order:
# reference databases, then the three production projects (hours).
#
# Pipelines are named projects in _targets.yaml (dada2, assign_taxo, cross_val,
# plus assign_taxo_mini / cross_val_mini for smoke tests on the mini_* databases
# with their own stores). run_project() (R/run_project.R) selects one, builds it
# and prunes the objects that no longer belong to it.
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

# Step 1: DADA2 denoising → d_asv (and OTU variants).
run_project("dada2")

# Step 2: Taxonomic assignment (all method × db combinations in parallel).
run_project("assign_taxo")

# Step 3: Cross-validation.
# WARNING: full run (cv_fold_number folds) takes many hours.
# The cross_val_mini project uses 2 folds / 200 sequences (config.R).
run_project("cross_val")
