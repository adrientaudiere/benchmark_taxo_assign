# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository. It describes what is on `make.R`'s path **today**; history lives in `Archives/` and in the ROADMAP done log.

## What this project is

A benchmark of taxonomic assignment methods (dada2 / sintax / lca / blastn) against six reference-database variants (UNITE and Eukaryome ITS/SSU, with `_Fungi`-filtered and cutadapt-trimmed `_cut` variants) for ITS fungal metabarcoding, evaluated on a mock community, by cross-validation, in silico and (planned) on a biological community. **It is not an R package** — it is a set of `{targets}` pipelines plus a Quarto analysis project. Do not run `devtools::check()`, `R CMD check`, or any R-package release tooling against this directory.

Everything is built around `phyloseq` objects and the sister pqverse packages (see *External code dependencies*).

## Layout

```
_targets.yaml        # five named {targets} projects (TAR_PROJECT)
config.R             # shared constants: primers, threads, db_list, db_meta, CV knobs, mini_db
make.R               # runs DB derivation + the three production projects in order
make_databases.R     # idempotent derivation of every reference-database variant (dbpq)
pipelines/           # dada2.R, assign_taxo.R, cross_val.R
R/                   # load_pqverse.R, values_map.R, autometric_helpers.R,
                     # combine_taxo_assignments.R, cross_val.R, cv_to_tidy.R,
                     # create_fake_pq_from_refseq.R, examples_cross_val.R
analysis/            # Quarto project: 00_setup.R, 01..06 chapters, in_silico_simulation.qmd, sandbox/
tests/               # Rscript-runnable testthat files
figures/             # notebook output (gitignored)
data/                # see data/README.md
store_dada2/ store_assign_taxo/ store_cross_val/   # targets stores (meta/ tracked)
docs/                # design diagram, external references
critique_fable/      # structural review of 2026-09-10 (evidence for ROADMAP section 1)
Archives/            # superseded scripts, kept for reference
```

## Running the pipelines

Pipelines are named projects in `_targets.yaml`; select one with `TAR_PROJECT` and call the targets verbs without `script =` / `store =`:

```r
Sys.setenv(TAR_PROJECT = "assign_taxo")
targets::tar_manifest()
targets::tar_visnetwork(targets_only = TRUE)
targets::tar_make()
targets::tar_prune()
# In a second R session at the project root, with the same TAR_PROJECT:
targets::tar_poll()
```

| Project | Script | Store | Purpose |
|---|---|---|---|
| `dada2` | `pipelines/dada2.R` | `store_dada2` | DADA2 denoising → `d_asv` |
| `assign_taxo` | `pipelines/assign_taxo.R` | `store_assign_taxo` | methods × DBs assignments → `d_all_taxo`, `benchmark_costs` |
| `cross_val` | `pipelines/cross_val.R` | `store_cross_val` | k-fold CV → `cv_results` |
| `assign_taxo_mini`, `cross_val_mini` | same scripts | `store_*_mini` (gitignored) | smoke tests on the `mini_*` databases |

`make.R` runs `derive_all_variants()` then the three production projects through `run_project()` (`tar_make()` + `tar_prune()`). `config.R` sets `mini_db <- grepl("_mini$", TAR_PROJECT)` and prints its value: the `*_mini` projects use the `mini_*` fasta files and write to their own stores, so a smoke test can never overwrite a production store. Target names are the same in both.

Warn the user before launching: a full `assign_taxo` run takes hours, `cross_val` many hours (`cv_fold_number`, `cv_fold_tested`, `cv_max_seq` in `config.R` are at their publication values; lower them for a smoke test). `dada2` needs cutadapt in the `cutadaptenv` conda env.

Every store also has a `session_info` target (`pqverse_versions()` + `sessioninfo::session_info()`) recording which checkout commits produced the run.

## Analysis (Quarto project in `analysis/`)

```
quarto render analysis                      # everything in _quarto.yml's render list
quarto render analysis/02_q1_methods.qmd    # one chapter
```

- `00_setup.R` — sourced first by every chapter: packages, `load_pqverse()`, `values_map`, constants (`tax_order`, `single_methods`, `preference_pattern`), helpers (`filter_default_settings()`, `save_fig()`, `cache_write()` / `cache_read()`, `plot_tc_metrics_mock()`).
- `01_load_and_clean.qmd` — reads `d_all_taxo` from `store_assign_taxo`, cleans taxonomy strings, builds `Gen_sp_*` and the consensus columns, computes `res_comp_tax` with `tc_metrics_mock()`, then writes `analysis/_cache/*.rds`. **Render it first**; the other chapters read the cache.
- `02_q1_methods`, `03_q2_databases`, `04_q3_consensus`, `05_m2_costs`, `06_cross_validation` — one manuscript question each; figures go to `figures/` via `save_fig()`.
- `in_silico_simulation.qmd` — D1b (InSilicoSeq + miaSim), exploratory.
- `sandbox/` — unfinished material, not rendered by default.
- Every chapter starts with `here::i_am("analysis/<file>.qmd")`: `_quarto.yml` is a `here` root criterion, so without it `here()` would resolve to `analysis/`.
- `execute-dir: project` means the working directory is `analysis/` during render; use `here()` for every path.

## Project dependencies

No `renv`. Packages come from the user's main R library; the pqverse packages come from their development checkouts (decision 10 in the ROADMAP), loaded by `R/load_pqverse.R`:

```r
source(here("R/load_pqverse.R"))
load_pqverse(c("MiscMetabar", "comparpq"))     # + "dbpq", "greenAlgoR", "tidypq" in the notebook
```

`PQVERSE_PKG_DIR` (env var, default `~/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg`) is the only path outside the project. Do not add `library("MiscMetabar")` or absolute `load_all()` paths; do not assume the packages are installed system-wide. Inside crew workers each target calls `load_pqverse()` again because checkouts are not installed packages (`tar_option_set(packages = ...)` lists only CRAN/Bioconductor packages).

External tools: `vsearch` and BLAST+ on `PATH`, cutadapt in the `cutadaptenv` conda env (`config.R::cutadapt_conda_prelude`), InSilicoSeq for D1b-A.

## External code dependencies

- **`MiscMetabar`** — `add_new_taxonomy_pq`, `assign_*`, `cutadapt_remove_primers`, `filter_trim`, `chimera_removal_vs`, `mumu_pq`, `asv2otu`, `sam_data_matching_names`, `track_wkflow`, `add_external_seq_pq`, `add_shuffle_seq_pq`, `add_dna_to_phyloseq`, `simplify_taxo`, `tc_bar`, `tc_circle`. dada2 is an *Import* of MiscMetabar, so `pipelines/dada2.R` attaches `dada2` itself.
- **`comparpq`** — `tc_metrics_mock`, `resolve_taxo_conflict`, `taxtab_replace_pattern_by_NA`, `rename_ranks_pq`.
- **`dbpq`** — `cutadapt_rm_primers_db`, `filter_db`, `format2sintax`, `format2dada2`, `find_vsearch`, `is_vsearch_installed` (used by `make_databases.R`). Helpers still local to this project and worth upstreaming are listed in `proposals_for_dbpq.md`; read it before adding a new helper.
- **`greenAlgoR`** — `ga_footprint()` for the CO₂eq figures (M2.2).
- **`tidypq`** — being adopted (ROADMAP decision 13, items S5.1–S5.2): use `tidypq::pq_to_tidy()` instead of `psmelt()`, and `tax_table_to_df()` / `mutate_taxa_pq()` / `filter_taxa_pq()` / `rename_samples_pq()` instead of editing phyloseq slots by hand. Call-site inventory in `critique_fable/06_tidypq_adoption.md`.

## Pipeline architecture

### `pipelines/dada2.R` → `store_dada2`

Paired-end DADA2 flow for the ITS-1F/ITS2 primer pair (Pauvert et al. 2018), `n_threads` and `seq_len_min` from `config.R`. Stages: cutadapt primer removal (conda env), `filter_trim` → `derepFastq` → `learnErrors` → `dada` → `mergePairs` → `makeSequenceTable` → `chimera_removal_vs` → length filter → sample renaming (`sam_data_matching_names`, `samp_` prefix) → seed taxonomy with `assignTaxonomy` on `refseq/dada2_format/Unite_Fungi.fasta` (rank `Phylum`, not `Phyla`) → `d_asv`, plus vsearch / clusterize OTU variants and their `mumu_pq` versions, `track_df`, `benchmark_costs_dada2`, `session_info`. Heavy phases are wrapped in `with_autometric()`.

### `pipelines/assign_taxo.R` → `store_assign_taxo`

- `d_asv_file` (`format = "file"` on `store_dada2/objects/d_asv`) → `d_asv`, so a rebuilt DADA2 store invalidates everything downstream.
- Negative controls: `d_asv_shuffled` (`add_shuffle_seq_pq`, prefix `fake_`) → `d_asv_for_assignation` (`add_external_seq_pq` with `data/data_raw/fake_ref/fake_ref_asv_100.fasta`, prefix `external_`). Taxa matching `^fake_|^external_` are the TN denominator of `tc_metrics_mock()`.
- `values_map <- build_values_map(dbs = db_list, mini_db = mini_db)` (`R/values_map.R`): 12 method/parameter rows × 6 DBs = 72 assignment targets, run in parallel by `tarchetypes::tar_eval` on a `crew_controller_group` (`dada2_ctrl` with one worker, `fast_ctrl` with `n_workers`). Each target calls `load_pqverse()` then `with_autometric(full_name, add_new_taxonomy_pq(...), dir = autometric_dir_assign)`.
- **`full_name` strings are the target names** (`dada2__Unite___0.5`, `blastn__EUK_ITS_v2___0.5...rel_majority...100`, …). `tests/test_values_map.R` pins them; changing `full_name_for()` orphans every store object.
- `d_all_taxo` = `tar_combine` of all assignments through `combine_taxo_assignments(d_asv_for_assignation, !!!.x)` (`R/combine_taxo_assignments.R`, detects new columns by `setdiff()` so it does not need to know suffixes).
- `benchmark_costs` aggregates the newest autometric log per target (`R/autometric_helpers.R`: `read_autometric_dir()` + `summarise_autometric_costs()`) and joins `values_map`.

### `pipelines/cross_val.R` → `store_cross_val`

`cv_values_map <- build_cv_values_map(...)`: 4 methods × 6 DBs × `remove_tested ∈ {TRUE, FALSE}` (standard / leaked). Each target runs `cross_val()` (`R/cross_val.R`) through `run_cv()` then `cv_to_tidy()`; `cv_results` binds everything. Known limitation: the `dada2` branch of `cross_val()` calls `assignTaxonomy` with `minBoot = 0` and never applies `min_bootstrap`; the `dada2_2steps` branch `stop()`s.

## Reference databases

`data/data_raw/refseq/` holds `dada2_format/` and `sintax_format/` trees with one fasta per variant; dada2 reads the former, sintax/lca/blastn the latter (`db_path_for()` in `R/values_map.R`). All variants derive from three source files (`Unite.fasta`, `Euk_ITS_v2.fasta`, `Euk_SSU_v2.fasta`) through `make_databases.R::derive_all_variants()`, idempotent (`force = TRUE` rebuilds). The benchmarked set is `db_list` in `config.R` (six DBs; decision 7 keeps only `EUK_SSU_v2_Fungi_cut` among the SSU variants); `db_meta` there maps each DB to its simplification type for the Q2 figures. `mini_*` files (first 10 000 records) exist for every variant and are used by the `*_mini` projects only (decision 4: never in the Results).

## Cost layer (M2)

`autometric::log_start()` only logs the process it runs in, so logging must start *inside* the crew worker: `with_autometric(phase, expr, dir)` writes `data/data_final/autometric/<pipeline>/<phase>__<timestamp>.txt` per target and per run; the aggregation keeps the newest file per phase. The legacy single-file logs (`data/data_final/autometric_log_*.txt`) are no longer written. The M2 figures in the store were computed before this fix (ROADMAP S1.1) and must be regenerated after the next production run.

## Tests

```
Rscript tests/test_combine_taxo_assignments.R   # combine == sequential chain accumulator
Rscript tests/test_values_map.R                 # target names identical to the legacy formula
```

Run from the project root. Add a test next to these for any helper under `R/`.

## Other notes

- `tar_option_set(seed = targets_seed)` (22) in every pipeline; `derive_fake_ref()` uses the same seed. Preserve it.
- `n_threads = 4`, `n_workers = 3` in `config.R`; the parallel pipeline can use up to `n_workers * n_threads + n_threads` cores.
- `store_*/meta/meta` are tracked on purpose (decision 9); `tar_prune()` in `make.R` keeps them small. `figures/` is not tracked (decision 8).
- `Archives/` holds files no longer wired into any pipeline (sequential `script_assign_taxo.R`, `some_bash_script`, `archive_benchmarking_taxonomy.R`, `R/cutadapt_rm_primers_db.R`). Do not edit unless asked.
- `ROADMAP.md` section 1 lists the remaining structural work (tests, docs, tidypq adoption); `critique_fable/` holds the evidence. `CONTEXT.md` is the glossary.
