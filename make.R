# Full analysis pipeline — run this file to (re)build everything in order.
#
# Pipelines are named projects in _targets.yaml (dada2, assign_taxo, cross_val,
# plus assign_taxo_mini / cross_val_mini for smoke tests on the mini_* databases
# with their own stores). run_project() selects one, builds it and prunes the
# objects that no longer belong to it (stale target names from older versions).
#
# Smoke test:   run_project("assign_taxo_mini"); run_project("cross_val_mini")
# Single step:  Sys.setenv(TAR_PROJECT = "assign_taxo"); targets::tar_make()

run_project <- function(name) {
  withr::with_envvar(c(TAR_PROJECT = name), {
    message("== targets project: ", name, " (store: ", targets::tar_config_get("store"), ")")
    targets::tar_make()
    targets::tar_prune()
  })
}

# Step 0: derive reference databases from the three source files in
#   data/data_raw/refseq/ (Unite.fasta, Euk_ITS_v2.fasta, Euk_SSU_v2.fasta).
#   Idempotent: skips files that already exist. Pass force = TRUE to rebuild.
source("make_databases.R")
derive_all_variants()

# Step 1: DADA2 denoising → d_asv (and OTU variants).
run_project("dada2")

# Step 2: Taxonomic assignment (all method × db combinations in parallel).
run_project("assign_taxo")

# Step 3: Cross-validation.
# WARNING: full run (cv_fold_number folds) takes many hours.
# Lower cv_fold_tested / cv_max_seq in config.R for a smoke test.
run_project("cross_val")
