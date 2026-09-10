# Context — benchmar_assign_taxo

Ubiquitous language for this project: the shared vocabulary used by the
developer and AI agents. Each entry maps a domain term or piece of jargon to a
plain definition, so conversations and code stay consistent and concise.

Idea from https://github.com/mattpocock/skills.

## Glossary

**benchmark** — Systematic comparison of four taxonomic-assignment methods (dada2, sintax, lca, blastn) across six reference-database variants (`db_list` in `config.R`) for ITS fungal metabarcoding. Performance is evaluated on a mock community, cross-validation, in silico simulations, and a real biological community (Bokulich 2020 four-dataset approach).

**assignment method** — Algorithm used to assign a taxonomic label to an ASV/OTU sequence. The four methods in scope: `dada2` (RDP naive Bayes classifier, `assignTaxonomy`), `sintax` (vsearch SINTAX bootstrap), `lca` (Lowest Common Ancestor), `blastn` (BLAST nucleotide search with a vote algorithm).

**reference database** — FASTA file providing labelled sequences against which query ASVs are classified. Main databases used: UNITE (ITS fungi), Eukaryome ITS, Eukaryome SSU. Each source is derived into multiple variants (full, Fungi-only `_Fungi`, primer-trimmed `_cut`, 99%-clustered) by `make_databases.R`.

**`values_map`** — Tibble built by `build_values_map()` (`R/values_map.R`) from `db_list`, enumerating all (method, parameter, db) combinations to benchmark: 72 rows. Each row carries `full_name` (the targets name, pinned by `tests/test_values_map.R`), `db_path` (format chosen by method, `mini_*` when `mini_db`), and the crew `controller`. `pipelines/assign_taxo.R` dispatches one target per row; `d_all_taxo` is assembled by `tarchetypes::tar_combine`. `build_cv_values_map()` is the cross-validation counterpart.

**`_targets.yaml` / `TAR_PROJECT`** — Five named targets projects: `dada2`, `assign_taxo`, `cross_val` (production stores) and `assign_taxo_mini`, `cross_val_mini` (same scripts on the `mini_*` databases, own gitignored stores). `Sys.setenv(TAR_PROJECT = "assign_taxo"); tar_make()`. `make.R::run_project()` wraps `tar_make()` + `tar_prune()`.

**`load_pqverse()`** — `R/load_pqverse.R`; loads the sister packages from their development checkouts under `PQVERSE_PKG_DIR` (decision 10). Called by every pipeline, inside every crew-worker target, and by `analysis/00_setup.R`. `pqverse_versions()` records version + git commit in the `session_info` target of each store.

**analysis chapters** — Quarto project in `analysis/`: `00_setup.R` (shared setup and helpers), `01_load_and_clean` (writes `analysis/_cache/*.rds`), `02_q1_methods`, `03_q2_databases`, `04_q3_consensus`, `05_m2_costs`, `06_cross_validation`, `in_silico_simulation.qmd`, `sandbox/`. Each chapter starts with `here::i_am()` because `_quarto.yml` is a `here` root criterion.

**mock community** — Lab-made community of known fungal taxa (Pauvert et al. 2019) used as the primary ground truth for computing TP/FP/FN/TN/MCC/F1/TAR/TDR. The known taxonomy is in `data/data_raw/metadata/taxo_mock.csv`.

**fake taxa** — Negative-control sequences injected before assignment: shuffled ASV sequences (`add_shuffle_seq_pq`, prefix `fake_`) and sequences from an external FASTA (`add_external_seq_pq`, prefix `external_`). Taxa matching `^fake_|^external_` are the TN denominator in `tc_metrics_mock()`.

**`tc_metrics_mock()`** — Function from `comparpq` that computes classification performance metrics (TP, FP, FN, TN, MCC, ACC, F1, PPV/TAR, TPR/TDR) for each (method × db × taxonomic rank) combination against a mock truth table.

**consensus voting / `resolve_taxo_conflict()`** — `comparpq` function that combines taxonomy assignments from multiple methods/databases into a single consensus call. Five strategies: `unanimity`, `consensus`, `abs_majority`, `rel_majority`, `preference`. Controlled by `strict` (boolean) and `nb_agree_threshold` (integer: 1, 2, or 3).

**cross-validation (CV)** — k-fold held-out evaluation of each (method, db) pair using the reference database itself as the source of truth. Implemented in `R/cross_val.R` (`cross_val()` and `cross_val_param()`), run via `pipelines/cross_val.R` → `store_cross_val`. Key parameters in `config.R`: `cv_fold_number`, `cv_fold_tested`, `cv_max_seq`.

**`remove_tested`** — Cross-validation parameter (`remove_tested_sequences = TRUE/FALSE`) that controls whether the sequences being classified are excluded from the reference database (TRUE = standard, non-leaked; FALSE = leaked, upper-bound estimate).

**in silico simulation** — Two sub-branches (D1b): InSilicoSeq generates simulated fastqs from a curated FASTA (sequencing errors included); miaSim generates community matrices from neutral/niche ecological models (no sequencing errors). Both yield a known ground truth for metric computation.

**`store_dada2` / `store_assign_taxo` / `store_cross_val`** — Three independent `{targets}` stores at the project root, one per production project of `_targets.yaml`; `make.R` runs all three in order. `store_assign_taxo` depends on `store_dada2/objects/d_asv` as a file target. The notebooks read them via `tar_read(..., store = here("store_assign_taxo"))`. `meta/meta` is tracked in git (decision 9), `objects/` is not.

**`autometric` / `with_autometric()`** — Per-target logging of wall time, peak RAM and mean CPU. Because `autometric::log_start()` only logs its own process, `with_autometric()` (`R/autometric_helpers.R`) starts the logger *inside* the crew worker and writes one file per target and per run under `data/data_final/autometric/<pipeline>/`. `summarise_autometric_costs()` keeps the newest file per phase; `benchmark_costs` (assignments) and `benchmark_costs_dada2` (denoising) are the aggregated targets.

**`config.R`** — Single shared-constants file sourced by the three pipelines, `make_databases.R` and `analysis/00_setup.R`. Contains primers (`fw_primer_sequences`, `rev_primer_sequences`), `n_threads`, `seq_len_min`, `prop_fake`, paths, conda prelude for cutadapt, `n_workers`, `targets_seed`, the CV knobs, `db_list`, `db_meta`, and `mini_db` (derived from `TAR_PROJECT`).

**`combine_taxo_assignments()`** — Helper in `R/combine_taxo_assignments.R` used by `tarchetypes::tar_combine` to merge per-(method, db) `tax_table` columns into one wide phyloseq object. Detects new columns via `setdiff(colnames(new_tt), colnames(base_tt))`; does not need to know per-row suffixes.

**`proposals_for_dbpq.md`** — Inventory of six helper functions still implemented locally in `make_databases.R` / `R/create_fake_pq_from_refseq.R` that fit `dbpq`'s scope and are candidates for upstreaming. Consult it before adding a new helper to avoid duplication.

**M1 metric set** — Primary performance metrics used in results figures: F1 + MCC (Hleap 2021) + TAR/TDR per taxonomic rank. TP/FP/FN/TN go to supplementary material only. Decided 2026-05-19.

**M2 environmental metrics** — Carbon-footprint and compute-cost layer: `wall_time_s`, `peak_resident_mb`, `mean_cpu_pct` from `benchmark_costs`, converted to CO₂eq via `greenAlgoR`.

**`cutadaptenv`** — Conda environment name containing `cutadapt`. Required by the DADA2 pipeline's primer-removal step; the script prepends `source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv &&` before each cutadapt call.

## Key decisions

**Four assignment methods, IdTaxa excluded (2026-05-19)** — The methods set is fixed at dada2, sintax, lca, blastn. IdTaxa is excluded because it requires a separate training step, making cost comparisons unfair. The `idtaxa` block remains commented out in `values_map`; the manuscript will include a one-line justification pointing to Murali et al. 2018.

**SSU database reduced to one variant (2026-05-19)** — Of four EUK_SSU variants, only `EUK_SSU_v2_Fungi_cut` (named `EUK_SSU_v1_9_3_Fungi_cut` at decision time) is retained. ITS sequences are a subregion of 18S SSU, so comparing all SSU variants would be redundant with the ITS database axis. The Fungi-filtered, primer-trimmed variant is the most directly comparable.

**Mini databases are smoke-test only (2026-05-19)** — `derive_mini()` in `make_databases.R` produces fast 10 000-line database subsets for iteration. They are excluded from all Results figures; use `dplyr::filter(!startsWith(db, "mini_"))` before any publication figure.

**Parallel pipeline over sequential** — `pipelines/assign_taxo.R` replaces the archived `Archives/script_assign_taxo.R`. Each (method, db) reads `d_asv_for_assignation` directly and runs in parallel via `crew_controller_local`; the sequential predecessor had each target depend on `previous_target`, creating a bottleneck.

**Biological validation dataset: Taudière et al. 2018 (2026-05-19)** — Tree endophyte ITS data from Taudière et al. (doi:10.1016/j.funeco.2018.07.008), with `site` and `height` as sample modalities. No objective ground truth exists; metrics are agreement-only (Jaccard, Bray-Curtis, richness). Pending: fastq availability check and SRA/Dryad/ENA fetch.

**Structural decisions (2026-09-10)** — Figures are not versioned (`figures/` gitignored); `store_*/meta/meta` stay tracked; the reference code is the dev checkout of every pqverse package (`load_pqverse()`); SRR30413326 joins as a second mock (D1d); tidypq is adopted in the notebook and the pipelines (ROADMAP S5.x); smoke tests run in the `*_mini` targets projects, never in a production store. Details in ROADMAP section 3b.

**`tar_option_set(seed = 22)`** — Fixed random seed set in both `pipelines/dada2.R` and `pipelines/assign_taxo.R`. Must be preserved when modifying targets so pipeline re-runs stay reproducible.
