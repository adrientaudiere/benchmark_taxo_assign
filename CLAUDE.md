# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A benchmark of taxonomic assignment methods (dada2 / sintax / lca / blastn, plus optional idtaxa) against several reference databases (UNITE, Eukaryome ITS/SSU, with/without cutadapt-trimmed and Fungi-only variants) for ITS fungal metabarcoding data. **It is not an R package** — it is a `{targets}` pipeline project plus a Quarto analysis notebook. Do not run `devtools::check()`, `R CMD check`, or any R-package release tooling against this directory.

The whole pipeline is built around `phyloseq` objects and depends on functions from sister projects in the pqverse workspace (see *External code dependencies* below).

## Running the pipelines

Three `{targets}` pipelines live at the project root, each with its own store. `make.R` calls the first two in sequence:

```r
tar_make(script = "script_dada2.R",                store = "store_dada2")
tar_make(script = "script_assign_taxo_parallel.R", store = "store_assign_taxo")
```

The cross-validation pipeline is run independently:

```r
tar_make(script = "script_cross_val.R", store = "store_cross_val")
```

Useful per-pipeline commands (replace `<script>`/`<store>` accordingly):

```r
tar_manifest(script = "script_assign_taxo_parallel.R")
tar_visnetwork(script = "script_assign_taxo_parallel.R", targets_only = TRUE)
tar_make(script = "script_assign_taxo_parallel.R", store = "store_assign_taxo")
tar_destroy("all", store = "store_assign_taxo/")
# In a second R session at the project root:
tar_poll(store = "store_assign_taxo")
```

Pipeline progress is also written to `data/data_final/autometric_log_dada2.txt` and `data/data_final/autometric_log_assign_taxo.txt` via the `autometric` package (warn the user: full assign-taxo runs are long — hours, sometimes days, depending on the database set).

The analysis lives in `analysis/benchmark.qmd` and reads targets from both stores via `tar_read(..., store = "store_assign_taxo")`. `In_silico_simulation.qmd` is a separate exploratory notebook for simulated communities.

## Project dependencies

The project no longer uses `renv` — `.Rprofile` is empty and the on-disk `renv/` / `renv.lock` are legacy artifacts. Packages are loaded from the user's main R library. The non-CRAN dependencies are the sister pqverse projects listed below (sourced/`load_all`'d by absolute path) plus `cutadapt` in a conda env named `cutadaptenv`.

## External code dependencies

The benchmark consumes three sister pqverse packages:

- **`dbpq`** — `library(dbpq)`. Provides `cutadapt_rm_primers_db`, `filter_db`, `format2sintax`, `format2dada2`, `get_file_extension`, `find_vsearch`, `is_vsearch_installed`. `make_databases.R` and `R/cross_val.R` call it via `dbpq::`. Functions that the benchmark still implements locally are listed in `proposals_for_dbpq.md` (candidates for upstreaming).
- **`MiscMetabar`** — still loaded by `devtools::load_all("~/Nextcloud/IdEst/Projets/MiscMetabar/")` and `library("MiscMetabar")` in the pipeline scripts. Provides `add_new_taxonomy_pq`, `assign_*`, `cutadapt_remove_primers`, `mumu_pq`, `asv2otu`, `sam_data_matching_names`, `track_wkflow`, `add_external_seq_pq`, `add_shuffle_seq_pq`, `add_dna_to_phyloseq`, etc.
- **`comparpq`** — still sourced by absolute path in `script_assign_taxo_parallel.R` and `analysis/benchmark.qmd` (`source("~/.../comparpq/R/compare_taxo.R")` and friends). Provides `tc_metrics_mock`, `resolve_taxo_conflict`, `taxtab_replace_pattern_by_NA`, `rename_ranks_pq`, `simplify_taxo`.

`tidypq` is *not* used by the benchmark — none of the verbs it exports (sample/taxa/occurrence/tree filters) match a current call site.

Note the inconsistency: `script_dada2.R` sources MiscMetabar's `R/` directly with `lapply(list.files(..., full.names = TRUE), source)`, while `script_assign_taxo_parallel.R` and the qmd use `devtools::load_all()`. Either is acceptable in this project, but you cannot assume MiscMetabar is installed system-wide.

## Pipeline architecture

### `script_dada2.R` → `store_dada2`

Standard paired-end DADA2 flow for the ITS-1F/ITS2 primer pair (Pauvert et al. 2018), with `n_threads = 4` and `seq_len_min = 200`. Key stages:

1. Primer removal via `cutadapt_remove_primers` — **requires the `cutadaptenv` conda env**; the call prepends `source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv &&` to the cutadapt command. If cutadapt isn't installed in that env on the current machine, this target will fail.
2. `filterAndTrim` → `derepFastq` → `learnErrors` → `dada` → `mergePairs` → `makeSequenceTable` → chimera removal (`chimera_removal_vs`) → length filter.
3. Sample renaming via `sam_data_matching_names` using `data/data_raw/metadata/sam_data.csv` and the `samp_` prefix.
4. Initial taxonomy assignment with `assignTaxonomy` against `data/data_raw/refseq/Unite_Fungi.fasta` (this is the seed taxonomy; the real benchmark in `script_assign_taxo.R` overwrites/augments it).
5. Three parallel post-clustering variants: ASV, vsearch OTU, idtaxa/clusterize OTU — each also passed through `mumu_pq` cleaning. Final targets: `d_asv`, `d_vs`, `d_idtaxa`, `d_asv_mumu`, `d_vs_mumu`, `d_idtaxa_mumu`, plus a `track_df` summary.

### `script_assign_taxo.R` → `store_assign_taxo`

Reads `d_asv` from `store_dada2` and runs a factorial of taxonomic-assignment methods × databases as a *sequential* chain:

- `methods` × `db` combinations are built in `values_map` (a tibble). Each row carries a unique `full_name`, a `previous_target` symbol pointing to the prior row, and a `db_path` (dada2 format for `method == "dada2"`, sintax format otherwise).
- The chain is realized by `tarchetypes::tar_eval(...)`. Each target calls `add_new_taxonomy_pq(previous_target, method, ref_fasta = db_path, suffix = paste0("_", full_name), ...)`, so each step adds new taxonomy columns to the phyloseq object produced by the previous step. **The chain is intentionally sequential**: the last target (suffixed `_all_taxo`) contains all assignments stacked together. Inserting or reordering rows in `values_map` shifts every downstream target.
- Before the chain, fake sequences are injected: `add_shuffle_seq_pq(d_asv, prop_fake = 0.5)` then `add_external_seq_pq(..., "data/data_raw/fake_ref/fake_ref_asv_100.fasta")`. These provide the negative-control taxa used by `tc_metrics_mock` in the analysis (taxa names matching `^fake_` or `^external_`).
- Two `if (FALSE)` blocks at the top are manual one-shot helpers (cutadapt all sintax-format DBs; build mini-DBs by taking the first 10 000 lines). Toggle to `TRUE` only when you actually need to regenerate those files.

### `analysis/benchmark.qmd`

Reads the final `*_all_taxo` target, cleans taxonomy strings (`taxtab_replace_pattern_by_NA` with patterns like `.*_sp$`, `.*_incertae_sedis`, `unclassified.*`), builds `Gen_sp_*` columns by concatenating Genus + Species per method, then computes consensus taxonomies with `resolve_taxo_conflict` (methods: `rel_majority`, `preference`, `unanimity`). Performance metrics (TP/FP/FN/FDR/TPR/PPV/F1/TN/MCC/ACC, à la Hleap et al. 2021) are computed by `tc_metrics_mock` against `data/data_raw/metadata/taxo_mock.csv`.

## Reference databases

`data/data_raw/refseq/` holds two parallel directory trees — `dada2_format/` and `sintax_format/` — with one fasta per database. The parallel pipeline chooses the format by method: dada2 uses `dada2_format/`, everything else (sintax/lca/blastn) uses `sintax_format/`. Every database variant is derived by `make_databases.R`, which delegates to `dbpq` wherever possible (`dbpq::format2sintax`, `dbpq::filter_db`, `dbpq::cutadapt_rm_primers_db`, `dbpq::find_vsearch`). The handful of derivations dbpq does not cover yet (parens stripping, inverse-pattern filtering, head-by-records, vsearch clustering, balanced subsetting) stay local and are listed in `proposals_for_dbpq.md`. Source `make_databases.R`, then call `derive_all_variants()`; every step is idempotent. The legacy `Archives/some_bash_script` is kept for historical reference only.

To work with a faster subset, switch `script_assign_taxo.R` to the mini DBs by uncommenting `mutate(db = paste0("mini_", db))` in the `values_map` pipeline and adjusting the `preference_pattern` in `benchmark.qmd` accordingly.

## Cross-validation tooling

`R/cross_val.R` defines `cross_val()` and `cross_val_param()` — k-fold cross-validation of a single assignment method on a reference fasta. `script_cross_val.R` wires it into a targets pipeline (→ `store_cross_val`) covering 4 methods × 9 DBs × 2 variants (remove_tested TRUE/FALSE). `R/cv_to_tidy.R` converts single-bootstrap `cross_val()` output into a long tibble. Worked examples live in `R/examples_cross_val.R`. The `dada2_2steps` branch in `cross_val()` still `stop()`s — never wired up.

Three constants in `config.R` control CV runs:
- `cv_fold_number` — total folds (default 10; publication value).
- `cv_fold_tested` — folds actually run (default 2 for smoke-testing; set to `cv_fold_number` for publication).
- `cv_max_seq` — cap on sequences sampled from each DB before folding (default 100 for fast smoke-tests; set to `NULL` for the full DB).

Known limitation: the `dada2` branch in `cross_val()` calls `assignTaxonomy` with `minBoot = 0` and never applies the `min_bootstrap` threshold post-hoc, so the bootstrap filter has no effect for dada2 cross-validation.

## Parallel pipeline

`make.R` now runs `script_dada2.R` → `store_dada2` and `script_assign_taxo_parallel.R` → `store_assign_taxo`. The sequential predecessor lives in `Archives/script_assign_taxo.R` (kept for reference; not on `make.R`'s path). In the parallel layout each (method, db) consumes `d_asv_for_assignation` directly instead of `previous_target`, and the final `d_all_taxo` is built by `tarchetypes::tar_combine` calling `combine_taxo_assignments(base, !!!.x)`.

- `script_assign_taxo_parallel.R` — adds a `crew_controller_local` (workers from `config.R::n_workers`) and a `benchmark_costs` target that aggregates `data/data_final/autometric_log_assign_taxo.txt` per target (phase tagged via `autometric::log_phase_set(full_name)` inside each assignment).
- `R/combine_taxo_assignments.R` — the combine helper. Detects new columns by `setdiff(colnames(new_tt), colnames(base_tt))`, so it does not need to know per-row suffixes.
- `tests/test_combine_taxo_assignments.R` — testthat fixture proving the combine output matches a `Reduce()`-style chain accumulator. Run with `Rscript tests/test_combine_taxo_assignments.R` from the project root.
- `config.R` — shared constants (primers, `n_threads`, `seq_len_min`, `prop_fake`, paths, conda prelude, `n_workers`, `targets_seed`, `cv_fold_number`, `cv_fold_tested`, `cv_max_seq`). Sourced by `script_assign_taxo_parallel.R` and `script_cross_val.R`; `script_dada2.R` still has inline copies and should be migrated next time that pipeline is touched.

When pointing `analysis/benchmark.qmd` at the parallel store, three edits are needed (the `_all_taxo` suffix no longer exists): read `tar_read(d_all_taxo, ...)` directly, drop the `gsub("_all_taxo", ...)` calls, and remove the `rename_ranks_pq(..., gsub("_all_taxo", "", ...))` block.

## Other notes

- `Archives/` holds files that are no longer wired into any pipeline: the sequential `script_assign_taxo.R`, the broken `some_bash_script`, the exploratory `archive_benchmarking_taxonomy.R`, and `R/cutadapt_rm_primers_db.R` (now superseded by `dbpq::cutadapt_rm_primers_db`). Keep them as historical reference; do not edit unless explicitly asked.
- `R/functions.R` used to mix three concerns (helper definition, fake-ref derivation as top-level side effects, exploratory `cross_val` runs at the bottom). It has been split: helper stays in `R/functions.R`, derivation moved to `make_databases.R` (which delegates to dbpq), `cross_val` to `R/cross_val.R`, exploratory runs to `R/examples_cross_val.R`. The file is now side-effect-free on source.
- `proposals_for_dbpq.md` lists six functions still implemented locally in `make_databases.R` and `R/functions.R` that fit dbpq's scope. Read it before adding a new helper here — the upstream candidate may already be enumerated.
- `tar_option_set(seed = 22)` is set in both pipelines; preserve it when modifying targets so re-runs stay reproducible.
- `n_threads = 4` is hard-coded in `script_dada2.R` (still inline) and in `config.R` (used by the parallel script).
