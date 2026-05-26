# Full analysis pipeline — run this file to (re)build everything in order.
#
# Step 0: derive reference databases from the three source files in
#   data/data_raw/refseq/ (Unite.fasta, Euk_ITS_v2.fasta, Euk_SSU_v2.fasta).
#   Idempotent: skips files that already exist. Pass force = TRUE to rebuild.
source("make_databases.R")
derive_all_variants()

# Step 1: DADA2 denoising → d_asv (and OTU variants).
targets::tar_make(script = "script_dada2.R", store = "store_dada2")

# Step 2: Taxonomic assignment (all method × db combinations in parallel).
targets::tar_make(script = "script_assign_taxo_parallel.R", store = "store_assign_taxo")

# Step 3: Cross-validation (independent pipeline, run separately).
# WARNING: full run (cv_fold_number folds) takes many hours.
# Keep cv_fold_tested = 2 in config.R for smoke-testing.
targets::tar_make(script = "script_cross_val.R", store = "store_cross_val")
