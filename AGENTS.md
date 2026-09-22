# AGENTS.md

Guidance for AI coding agents (and for any human arriving in this folder). It
describes what is on `make.R`'s path **today**; history lives in `Archives/` and
in [`HISTORY.md`](HISTORY.md). `CLAUDE.md` is a pointer to this file.

**Read order.** This file (how to work here) → [`CONTEXT.md`](CONTEXT.md)
(glossary) → [`ROADMAP.md`](ROADMAP.md) (the todo, one line per item) →
`docs/objectives_design.md` (the design, and the only place a design fact is
*defined*). `HISTORY.md` holds every ticked item, dated narrative and past
decision: read it for the measurements behind a ROADMAP line.

**Two review documents of 2026-09-22**, both addressed but not applied:
[`critique_2026-09-22.md`](critique_2026-09-22.md) (aims and code, with the
developer's framing answer in its addendum) and
[`literature_2026-09-22.md`](literature_2026-09-22.md) (positioning against the
six benchmarks this project must cite). Their proposals live as `C22-…` items in
`ROADMAP.md` sections 2 and 4. **Do not apply a critique item on its own
authority**: the developer triages them.

## What this project is

A benchmark of taxonomic assignment methods (dada2 / sintax / lca / blastn) against nine reference databases (`config.R::benchmark_dbs`: UNITE with and without singletons and EUKARYOME ITS, each whole, Fungi-only and Fungi plus non-fungal representatives) for ITS fungal metabarcoding, evaluated on two mock communities (Pauvert et al. 2019, assigned; Tedersoo et al. 2026, denoised in `store_dada2_tedersoo_illumina`, not assigned yet), by cross-validation, in silico and (planned) on biological datasets. **It is not an R package** — it is a set of `{targets}` pipelines plus a Quarto analysis project. Do not run `devtools::check()`, `R CMD check`, or any R-package release tooling against this directory.

Everything is built around `phyloseq` objects and the sister pqverse packages (see *External code dependencies*).

**Scope (developer, 2026-09-17): taxonomic assignment of Fungi, and nothing else.** No community ecology: the miaSim / Hubbell branch is dropped, and biological datasets are read only through `NA_real`, `NA_fake` and `ext_fungi`.

## Design of 2026-09-17: what is built, what is not

The design was redrawn on 2026-09-17 and is **implemented**: items 0.1–0.10 are done (`HISTORY.md` §0) — the nine-database grid, the `fungi + rep` databases, the OTU input, the lca `lca_cutoff` and blastn `min_id` sweeps, the per-dataset ITSx region, then the per-unit truth of the mocks (`R/mock_truth.R`), the Hleap scoring in comparpq (`tc_metrics_unit()`), the chapters reading it and the documents. **Not built** (ROADMAP 0.11): the cross-validation databases, the final zoom and the discussion notes. **No production run has been made on the new grid**: `store_assign_taxo` still holds the 2026-09 results of the former seven-database grid, so every figure must be regenerated. Before changing `config.R`, `R/values_map.R`, a pipeline or a chapter, read:

- **`docs/objectives_design.md`** — the source of truth for the design: inputs, controls, databases, methods, parameters and metrics per objective; its §4 numbers the developer's decisions 1–27 (not the same numbering as the decisions of `HISTORY.md` §3–3d: always write "objectives_design decision N"), and a new decision is written there first;
- **`docs/hleap_2021_metrics.md`** §5–6 — per-unit truth and Hleap-style scoring;
- **`HISTORY.md` §0** — the implementation plan, item by item, with what each one measured.

`docs/experimental_design.md` was updated on 2026-09-22 (item 0.10) and now reads 96 computations; **its measured costs of §3 were read on the former 44-computation grid** and describe less than half the current work — kept as an order of magnitude because they are the only figures anyone has.

In short: databases `[EUK_ITS_v2.1 | Unite_all_20250219 | Unite_s_all_20250219] × [all | fungi | fungi + rep]`, no `_cut` database in the grid (`config.R::zoom_dbs` keeps `EUK_ITS_v2.1_Fungi_cut` for a final zoom on the best combinations, Pauvert only); inputs ASV, OTU (`d_vs` assigned directly, no MUMU, no ITSx) and ITSx ASV (fungi and fungi + rep only), i.e. 96 computations per dataset (36 + 36 + 24) and 504 rows; lca `lca_cutoff` ∈ {0.8, 0.9, 1} from one stored vsearch hit table; blastn `min_id` ∈ {90, 92, 95, 97}; per-unit truth built from the Sanger sequences of the Pauvert and Tedersoo mocks (`R/mock_truth.R`), one confusion-matrix cell per unit and rank as in Hleap et al. 2021 (`comparpq::tc_metrics_unit()`).

Chapter 01 now scores three inputs against the per-unit truth (ROADMAP 0.8a, 2026-09-22): `unit` = `ASV` (`d_all_taxo`, raw and ITSx columns, plus the consensus columns), `OTU` (`d_all_taxo_otu`, the OTUs assigned directly) and `OTU_derived` (`resolve_otu_taxonomy()` on the ASV assignments, kept so Q5 compares "assign then cluster" with "cluster then assign"). It writes one long `res_comp_tax` with `dataset`, `unit` and `external_scoring` columns, plus `res_abundance` (the M1.5 sweep). **Chapters 02 to 05 must read it through `default_view()`** (`00_setup.R`), which keeps one dataset, one input and one reading of the consensus columns — without it they average the three inputs and count every consensus twice, silently.

## Framing of the manuscript (2026-09-22)

Decided by the developer and recorded in `HISTORY.md` §3d (decisions 37–40) and
`literature_2026-09-22.md` §4. Useful whenever a change touches what the paper
claims:

- **Orsholm et al. 2026 is the foundation, not a competitor.** For fungal ITS the algorithm-family question is settled (composition-based classifiers win), so the four methods are a **deliberate zoom**, not a gap. Honest nuance: `lca` is search-based like `blastn`, so the zoom is 2 composition-based (dada2, sintax) × 2 similarity-search (lca, blastn), with blastn as the abstaining contrast.
- **The complexity goes where the literature has none**: the composition of the reference database (Q2) and voting across classifications (Q3). **Q3 is a primary axis**, not a supplementary one.
- **A vote is only as good as its voters are independent.** `Unite_all` ⊂ `Unite_s_all`, and lca and blastn read the same reference, so a naive consensus counts shared evidence twice (ROADMAP C22-P6, C22-P21). Never present a Q3 gain as a voting result without that check.
- **The harness and part of the benchmarked software share an author** (MiscMetabar, comparpq). The engines (vsearch, BLAST+, dada2) are third-party; commits are pinned by the `session_info` targets; the paper declares it (ROADMAP C22-N1).


## Layout

```
_targets.yaml        # seven named {targets} projects (TAR_PROJECT)
config.R             # shared constants: primers, threads, blastn_min_cover, ITSx, db_list, db_meta, CV knobs, mini_db
make.R               # sourcing it runs DB derivation + the three production projects (hours)
make_databases.R     # idempotent derivation of every reference-database variant (dbpq)
pipelines/           # dada2.R, dada2_bio.R (dataset projects), assign_taxo.R, cross_val.R
R/                   # load_pqverse.R, run_project.R, values_map.R, assign_compute.R, itsx.R, autometric_helpers.R,
                     # combine_taxo_assignments.R, cross_val.R, cv_queries.R, cv_to_tidy.R,
                     # otu_taxonomy.R, mock_truth.R, score_assignments.R,
                     # create_fake_pq_from_refseq.R, examples_cross_val.R (sourced by nothing, C22-S5)
analysis/            # Quarto project: 00_setup.R, 01..06 chapters, in_silico_simulation.qmd, sandbox/
tests/               # Rscript-runnable testthat files
figures/             # notebook output (gitignored)
data/                # see data/README.md
store_dada2/ store_assign_taxo/ store_cross_val/   # targets stores (meta/ tracked)
docs/                # objectives_design.md (decided design, read first), hleap_2021_metrics.md (scoring and
                     # per-unit truth), experimental_design.md (ITS regions; its measured costs predate the grid),
                     # reference_databases.md, diagnostics_cv_and_databases.md,
                     # README.md (index), design diagram, external references
critique_fable/      # structural review of 2026-09-10 (evidence for the S-items of ROADMAP section 5)
critique_2026-09-22.md   # review of the aims and of the code size; its addendum holds the developer's framing
literature_2026-09-22.md # the six benchmarks to position against, and where this project sits
ROADMAP.md HISTORY.md    # the todo (one line per item) and the archive (done items, narratives, decisions)
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
| `dada2_tedersoo_illumina` | `pipelines/dada2_bio.R` | `store_dada2_tedersoo_illumina` | DADA2 on the Tedersoo et al. 2026 Illumina mock (built: 238 ASVs, `d_vs` 180 OTUs) |
| `assign_taxo_tedersoo_illumina` | `pipelines/assign_taxo.R` | `store_assign_taxo_tedersoo_illumina` | same grid as `assign_taxo`, ITSx region `full` (`bio_datasets$itsx_region`); not run yet, and no truth table until ROADMAP 0.6 |

The dataset projects come in pairs (`dada2_<dataset>`, `assign_taxo_<dataset>`, ROADMAP D1c); `config.R::bio_datasets` holds each dataset's fastq folder, sample file and primers, and `config.R` derives the dataset from `TAR_PROJECT`.

`make.R` runs `derive_all_variants()` then the three production projects through `run_project()` (`R/run_project.R`: `tar_make()` + `tar_prune()`). To build one project, `source("R/run_project.R"); run_project("assign_taxo_mini")`: never `source("make.R")` for that, it runs everything. `config.R` sets `mini_db <- grepl("_mini$", TAR_PROJECT)` and prints its value: the `*_mini` projects use the `mini_*` fasta files and write to their own stores, so a smoke test can never overwrite a production store. Target names are the same in both.

Warn the user before launching: `derive_all_variants()` takes about 40 min when it rebuilds the `fungi + rep` databases (`force = TRUE`, mostly the vsearch pass that keeps the representatives away from the controls), a full `assign_taxo` run takes hours, `cross_val` many hours (`config.R`: 10 folds and `cv_max_seq = 5000` for the production projects; the `*_mini` projects switch to 2 folds and 200 sequences automatically through `mini_db`). `dada2` and `cross_val` (query trimming) need cutadapt in the `cutadaptenv` conda env.

Every store also has a `session_info` target (`pqverse_versions()` + `sessioninfo::session_info()`) recording which checkout commits produced the run.

## Analysis (Quarto project in `analysis/`)

```
quarto render analysis                      # everything in _quarto.yml's render list
quarto render analysis/02_q1_methods.qmd    # one chapter
```

- `00_setup.R` — sourced first by every chapter: packages, `load_pqverse()`, `values_map`, constants (`tax_order`, `single_methods`, `preference_pattern`), helpers (`filter_default_settings()`, `save_fig()`, `cache_write()` / `cache_read()`, `plot_tc_metrics_mock()`).
- `01_load_and_clean.qmd` — reads `d_all_taxo` from `store_assign_taxo`, cleans taxonomy strings, builds `Gen_sp_*`, the OTU object `d_otu_taxo` (Q5: `resolve_otu_taxonomy()` with the `otu_clusters` target) and the consensus columns, scores the three inputs against the per-unit truth with `R/score_assignments.R` (which calls `comparpq::tc_metrics_unit()` once per `external_scoring` group), and writes `res_comp_tax` and `res_abundance` to `analysis/_cache/*.rds`. **Render it first**; the other chapters read the cache.
- `02_q1_methods`, `03_q2_databases`, `04_q3_consensus`, `05_m2_costs`, `06_cross_validation` — one manuscript question each; figures go to `figures/` via `save_fig()`.
- `in_silico_simulation.qmd` — D1b (InSilicoSeq), exploratory: the input fasta drawn from `config.R::insilico_source`, the Docker `iss generate` recipe, the targets projects still to create, and the per-record truth parsed from the UNITE headers. The miaSim / Hubbell branch was removed on 2026-09-22 (ROADMAP 0.9).
- `sandbox/` — unfinished material, not rendered by default.
- Every chapter starts with `here::i_am("analysis/<file>.qmd")`: `_quarto.yml` is a `here` root criterion, so without it `here()` would resolve to `analysis/`.
- `execute-dir: project` means the working directory is `analysis/` during render; use `here()` for every path.

## Project dependencies

No `renv`. Packages come from the user's main R library; the pqverse packages come from their development checkouts (decision 10 in the ROADMAP), loaded by `R/load_pqverse.R`:

```r
source(here("R/load_pqverse.R"))
load_pqverse(c("MiscMetabar", "comparpq"))     # + "dbpq", "greenAlgoR", "tidypq" in the notebook
```

`PQVERSE_PKG_DIR` (env var, default `~/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg`) is the only path outside the project. Do not add `library("MiscMetabar")` or absolute `load_all()` paths; do not assume the packages are installed system-wide.

**Rule for targets that run on crew workers** (`pipelines/assign_taxo.R`, `pipelines/cross_val.R`): workers only attach the CRAN/Bioconductor packages listed in `tar_option_set(packages = ...)` (`here`, `phyloseq`, `dplyr`, `tidyr`, `tibble`, …), never the pqverse checkouts. So a target either (a) calls `load_pqverse()` as its first statement (the `tar_eval` assignment and CV targets do), or (b) is declared `deployment = "main"` and runs in the main process where the checkouts are already loaded (every other target: store loading, negative controls, `tar_combine`, `benchmark_costs`, `session_info`). A "could not find function" error from `tar_make()` means a target broke this rule; if the missing function is from a CRAN package, add it to `packages`.

External tools: `vsearch` and BLAST+ on `PATH`, cutadapt in the `cutadaptenv` conda env (`config.R::cutadapt_conda_prelude`), ITSx 1.1.3 in the `itsxenv` conda env (`config.R::itsx_conda_prelude`; `conda create -n itsxenv -c conda-forge -c bioconda itsx`), InSilicoSeq for D1b-A.

## External code dependencies

- **`MiscMetabar`** — `add_new_taxonomy_pq`, `assign_*`, `cutadapt_remove_primers`, `filter_trim`, `chimera_removal_vs`, `mumu_pq`, `asv2otu`, `sam_data_matching_names`, `track_wkflow`, `add_dna_to_phyloseq`, `simplify_taxo`, `tc_bar`, `tc_circle`. dada2 is an *Import* of MiscMetabar, so `pipelines/dada2.R` attaches `dada2` itself and `R/cross_val.R` calls `dada2::assignTaxonomy()` (crew workers never attach it).
- **`comparpq`** — `tc_metrics_mock`, `tc_metrics_unit` / `tc_metrics_unit_vec` (per-unit scoring as in Hleap et al. 2021, added 2026-09-22 for item 0.7; **this is what the chapters call** — `tc_metrics_mock()` is kept for a mock that only comes with a list of expected taxa), `resolve_taxo_conflict`, `taxtab_replace_pattern_by_NA`, `rename_ranks_pq`, `add_shuffle_seq_pq`, `add_external_seq_pq` (the two negative-control injectors live here, not in MiscMetabar; checked 2026-09-16 in `comparpq/R/fake_creation.R`).
- **`dbpq`** — `cutadapt_rm_primers_db` (reverse-complements `primer_rev` by default since 2026-09-15; also used by `R/cv_queries.R` on the CV queries), `filter_db`, `format2sintax`, `format2dada2`, `find_vsearch`, `is_vsearch_installed` (used by `make_databases.R`). Helpers still local to this project and worth upstreaming are listed in `proposals_for_dbpq.md`; read it before adding a new helper.
- **`greenAlgoR`** — `ga_footprint()` for the CO₂eq figures (M2.2).
- **`tidypq`** — being adopted (ROADMAP decision 13, items S5.1–S5.2): use `tidypq::pq_to_tidy()` instead of `psmelt()`, and `tax_table_to_df()` / `mutate_taxa_pq()` / `filter_taxa_pq()` / `rename_samples_pq()` instead of editing phyloseq slots by hand. Call-site inventory in `critique_fable/06_tidypq_adoption.md`.

## Pipeline architecture

### `pipelines/dada2.R` → `store_dada2`

Paired-end DADA2 flow for the ITS-1F/ITS2 primer pair (Pauvert et al. 2018), `n_threads` and `seq_len_min` from `config.R`. Stages: cutadapt primer removal (conda env), `filter_trim` → `derepFastq` → `learnErrors` → `dada` → `mergePairs` → `makeSequenceTable` → `chimera_removal_vs` → length filter → sample renaming (`sam_data_matching_names`, `samp_` prefix) → seed taxonomy with `assignTaxonomy` on `refseq/dada2_format/<seed_taxonomy_db>.fasta` (`config.R`; rank `Phylum`, not `Phyla`) → `d_asv`, plus vsearch / clusterize OTU variants and their `mumu_pq` versions, `track_df`, `benchmark_costs_dada2`, `session_info`. Heavy phases are wrapped in `with_autometric()`.

### `pipelines/assign_taxo.R` → `store_assign_taxo`

- `d_asv_file` (`format = "file"` on `store_dada2/objects/d_asv`) → `d_asv`, so a rebuilt DADA2 store invalidates everything downstream.
- Negative controls: `d_asv_shuffled` (`add_shuffle_seq_pq`, prefix `fake_`) → `d_asv_for_assignation` (`add_external_seq_pq` with `data/data_raw/fake_ref/fake_ref_asv_100.fasta`, prefix `external_`; drawn since 2026-09-17 only from the 13 non-fungal kingdoms retained in every release, `config.R::rep_kingdom_min_records`, `external_control_kingdoms()`). The fasta is the file target `fake_ref_file`: regenerating it with `derive_fake_ref(force = TRUE)` invalidates every assignment (`derive_all_variants()` alone keeps an existing file; docs §6). Taxa matching `^fake_|^external_` are the controls; `external_scoring` decides whether the `external_` ones enter the confusion matrix (`_Fungi` databases) or are reported through `ext_fungi` / `ext_correct` (whole releases and `_Fungi_rep`).
- ITSx input (ROADMAP Q4, decision 19): `itsx_asv` (`run_itsx()` in `R/itsx.R`, ITSx in the `itsxenv` conda env) → `itsx_dropped_taxa` (`itsx_duplicated_taxa()`: ASVs whose ITS1 equals that of a more abundant ASV, decision 20) → `d_asv_common` (the ASVs of both inputs) → `d_asv_itsx` (the ITSx region replaces each detected ASV sequence: `config.R::itsx_region`, ITS1 for Pauvert, the full ITS for a dataset whose `bio_datasets$itsx_region` is `full`, read from ITSx's `<prefix>.full.fasta`; undetected ASVs keep theirs; ITSx runs with all organism profiles because it only trims flanks, it must not act as a Fungi filter) → `d_asv_itsx_shuffled` → `d_asv_itsx_for_assignation`. Both shuffles run under `withr::with_seed(targets_seed, …)`, so the two inputs have the same `fake_1..n` and `external_*` taxa, which `combine_taxo_assignments()` requires.
- OTU input (objectives_design decision 4, ROADMAP 0.4): `d_vs_file` (`format = "file"` on `<dada2_store>/objects/d_vs`, the 97 % vsearch OTUs, no MUMU) → `d_vs` → `d_vs_shuffled` → `d_vs_for_assignation`, with the same `prop_fake`, seed and external sequences as the ASV inputs, but its own `fake_` taxa.
- `values_map <- build_values_map(dbs = db_list, mini_db = mini_db)` (`R/values_map.R`): 21 method/parameter rows (dada2 and sintax `min_bootstrap` ∈ {0.4, 0.5, 0.6}; lca `lca_cutoff` ∈ {0.8, 0.9, 1}; blastn 3 votes × `min_id` ∈ {90, 92, 95, 97}; `build_methods_grid()`) × (9 DBs with the raw ASVs + 9 with the OTUs + the 6 Fungi-filtered DBs of `config.R::itsx_db_list` with the ITSx input; `preprocess`: `none`, `otu`, `itsx`; column `input_pq`) = 504 assignment targets, but only 96 computations (ROADMAP S8.2, `R/assign_compute.R`): one `compute_<method>__<db>` (or `itsx_compute_…`, `otu_compute_…`) target per method × database × input runs `compute_assignment()` on a `crew_controller_group` (`dada2_ctrl`: one worker with `assign_threads_dada2` threads; `fast_ctrl`: `n_workers` workers with `assign_threads_fast`), inside `load_pqverse()` and `with_autometric(compute_name, …)`; each `full_name` target then runs `derive_assignment()` in the main process (dada2 / sintax threshold; lca at its `lca_cutoff` from the stored vsearch hits through `MiscMetabar::assign_vsearch_lca(hits_table = )`, identical to vsearch `--lcaout` at that cutoff; blastn `min_id`, `min_cover` and vote through `MiscMetabar::assign_blastn(blast_table = )`). sintax and lca are called with `clean_pq = FALSE`, otherwise they skip the negative controls (B19). Each computation reads its reference through a file target (`ref_file_targets`: one `ref_sintax__<db>` / `ref_dada2__<db>` per fasta, 18 targets; same in `pipelines/cross_val.R`). `d_all_taxo` combines the ASV and ITSx rows on `d_asv_for_assignation`, `d_all_taxo_otu` the OTU rows on `d_vs_for_assignation`.
- Post-clustering from member ASVs (Q5, decision 25), kept beside the directly assigned OTUs so Q5 compares "assign then cluster" with "cluster then assign" (`unit = "OTU_derived"`): `otu_clusters` (`otu_membership()` in `R/otu_taxonomy.R`: 97 % vsearch clusters of `d_asv_common`, with the most abundant member as archetype); chapter 01 merges `d_all_taxo` into OTUs and keeps, per assignment column, the value all member ASVs agree on (`resolve_otu_taxonomy()`, `unanimity`, NA ignored).
- **`full_name` strings are the target names** (`dada2__Unite_s_all_20250219___0.5`, `lca__EUK_ITS_v2.1___0.8`, `blastn__EUK_ITS_v2.1___0.5...rel_majority...100...95`, `itsx_sintax__EUK_ITS_v2.1_Fungi___0.5`, `otu_lca__Unite_all_20250219___1`, …). `tests/test_values_map.R` pins them; changing `full_name_for()` orphans every store object.
- `d_all_taxo` = `tar_combine` of all assignments through `combine_taxo_assignments(d_asv_for_assignation, !!!.x)` (`R/combine_taxo_assignments.R`, detects new columns by `setdiff()` so it does not need to know suffixes).
- `benchmark_costs` aggregates the newest autometric log per computation (`R/autometric_helpers.R`: `read_autometric_dir()` + `summarise_autometric_costs()`) and joins `values_map` on `compute_name`: the three rows of one computation share its cost, and `phase` holds the `full_name` the chapters join on.

### `pipelines/cross_val.R` → `store_cross_val`

`cv_values_map <- build_cv_values_map(dbs = cv_db_list, ...)`: 4 methods × the 7 DBs of `config.R::cv_db_list` (the former grid, kept until the CV databases are chosen from the mock results, objectives_design decision 27) × `remove_tested ∈ {TRUE, FALSE}` (standard / leaked) = 56 targets. The queries are trimmed to the ITS1F–ITS2 amplicon (decision 24, `R/cv_queries.R`): `cross_val()` draws `cv_oversample × cv_max_seq` records, `trim_cv_queries()` cuts them with `dbpq::cutadapt_rm_primers_db()` (reverse-primer site required, `primer_min_overlap`), and `cv_select_queries()` keeps the drawn records up to the `cv_max_seq`-th one holding the site; records without it are never queried but stay in the training part, which always holds the records of the benchmarked database. A `_Fungi_cut` database no longer holds the primer sites: its queries are drawn and trimmed from its `_Fungi` source (`query_fasta`, columns `query_db` / `query_db_path` / `query_ref_file` of `cv_values_map`) and paired by name with the `_cut` records. dada2-format headers carry no identifier: `cross_val()` names the dada2 records with the identifiers of the sintax-format file of the same database (`id_fasta`, `query_id_fasta`, `record_ids()`; same records in the same order) and keeps the dada2 headers for the truth table and the training fasta. Before 2026-09-15 its name deduplication left dada2 one record per taxonomy string (ROADMAP B22). The metrics compare the assignment and truth tables row by row: `cv_align_rows()` puts the sintax and lca rows in the query order (vsearch `--sintax` with several threads does not keep it; every sintax CV result before 2026-09-15 is invalid, B23) and `cross_val()` stops when a method returns another number of rows than queries. Each target runs `cross_val()` (`R/cross_val.R`) through `run_cv()` then `cv_to_tidy()`; `cv_results` binds everything. **The reference is the drawn pool, not the database**: `cross_val()` keeps `cv_oversample × cv_max_seq` records (with the production settings, 5 413 of the 156 820 of `Unite_all_20250219` and 6 284 of the 1 877 003 of `EUK_ITS_v2.1`: the pool size follows the ITS2-site rate of the database, not its size), which is the documented behaviour of `max_seq` but makes every CV target search a small reference — blastn abstains (87–95 % NA), lca and sintax answer with the nearest available taxon (ROADMAP B24). `cross_val()` applies `min_bootstrap` to dada2 (`minBoot = 100 * min_bootstrap`) and to sintax (passed to `assign_sintax()`); before 2026-09-11 dada2 CV used no threshold and sintax always 0.5. Each CV target uses `config.R::cv_threads` threads (`nproc` for sintax / lca / blastn, `multithread` for dada2; keep `n_workers * cv_threads` within the cores). The `dada2_2steps` branch `stop()`s. ITSx does not apply to cross-validation (queries are reference sequences).

## Reference databases

Full procedure, header formats and how to move to a new release: **`docs/reference_databases.md`**. Summary:

- `config.R::reference_sources` lists the general FASTA releases with versioned names (`Unite_s_all_20250219` = UNITE 19.02.2025 eukaryotes with singletons, `Unite_all_20250219` = without singletons, `EUK_ITS_v2.1`), their URL and DOI. The SSU release was dropped on 2026-09-15 (decision 7 revised: ITS1 is not part of the 18S gene). `config.R::benchmark_dbs` lists the benchmarked databases (source × simplification: full, `_Fungi`, `_Fungi_rep`), `config.R::zoom_dbs` the databases built for the final zoom only (`EUK_ITS_v2.1_Fungi_cut`); `db_list`, `cv_db_list`, `db_meta`, `seed_taxonomy_db`, `fake_ref_source` and `preference_db` come from there. **Never hard-code a database name elsewhere**: pipelines, chapters and tests read `config.R`.
- `make_databases.R` step 1, `download_reference_sources()`: `dbpq::download_unite_db(url, extract = TRUE)` / `dbpq::download_eukaryome_db(url)` → `data/data_raw/refseq/sources/<source>.fasta` + `sources/manifest.csv` (URL, DOI, md5, date).
- Step 2, `derive_all_variants()`: `derive_dada2()` / `derive_sintax()` → `dada2_format/<source>.fasta` and `sintax_format/<source>.fasta` (both go through `general_headers_fixed()`: the UNITE `|k__` separator that makes dbpq drop the kingdom, and EUKARYOME name qualifiers such as `Lactarius(Fungi)`, `(Candida)` or `Mortierella.s.str` that break sintax bootstraps and mock matching), `derive_kingdom_only()` (kingdom exactly Fungi, anchored pattern), `derive_cutadapted()` (trimmed at ITS1F and at the reverse complement of ITS2 when found, records without primer sites kept, records shorter than `cut_min_length` after trimming dropped: `config.R::cut_discard_untrimmed`, `primer_min_overlap`, `cut_min_length`; ROADMAP S6.9), fake reference, `derive_fungi_rep()` (`_Fungi` records + `config.R::rep_share` × the fungal record count of non-fungal records, rounded up, drawn from the full release of the same source over the kingdoms of `retained_kingdoms()` by `water_fill()`, never a record within `config.R::rep_max_identity_to_external` (0.97) of an `external_` control, excluded from the candidates by one vsearch pass, candidates as queries and controls as database, which costs about 35 min for the three databases; same draw written in both formats; list in `sources/<db>_reps.csv`; its `mini_*` files hold 10 000 Fungi records + 500 representatives), `mini_*` subsets. Idempotent; `force = TRUE` rebuilds. Each reference fasta is a file target (`ref_dada2__<db>` / `ref_sintax__<db>`, column `ref_file` of `values_map`), so a database rebuilt in place reruns exactly the targets that read it (ROADMAP S7.9).
- dada2 reads `dada2_format/`, sintax/lca/blastn read `sintax_format/` (`db_path_for()` in `R/values_map.R`). `mini_*` files are used by the `*_mini` projects only (decision 4: never in the Results).

**Legacy files deleted (2026-09-11).** The unversioned reference files (`data/data_raw/refseq/{Unite,Unite_RefS,Euk_ITS_v2,Euk_SSU_v2}.fasta` and the unversioned files of `dada2_format/` / `sintax_format/`) lost the UNITE kingdom in sintax headers (`…|k__Kingdom;tax=p:…`, ranks shifted by one) and let non-Fungi records into the `_Fungi` files; do not recreate them. Every store built between 2026-05-26 and 2026-09-10 also ran all targets on the first database (`ifelse()` recycling in `values_map`, fixed). See `critique_fable/07_results_validity.md`.

## Cost layer (M2)

`autometric::log_start()` only logs the process it runs in, so logging must start *inside* the crew worker: `with_autometric(phase, expr, dir)` writes `data/data_final/autometric/<pipeline>/<phase>__<timestamp>.txt` per target and per run; the aggregation keeps the newest file per phase. The legacy single-file logs (`data/data_final/autometric_log_*.txt`) are no longer written. The M2 figures in the store were computed before this fix (ROADMAP S1.1) and must be regenerated after the next production run.

Measured on the production store (2026-09-16, `docs/experimental_design.md` §3): 44 computations, **2.15 h of target time** (1 h 34 wall), of which dada2 is 1.14 h (53 %) with a **53.1 GB** peak against 1.2–1.3 GB for sintax, lca and blastn. **The per-database split of that same table is not usable**: because the aggregation keeps the newest log per phase and the 2026-09-14 and 2026-09-15 runs are mixed, `EUK_ITS_v2.1_Fungi` comes out at 0.06 h against 1.19 h for its own trimmed variant, which holds the same records shortened. Totals are reliable, per-database figures are not, until a clean run is remeasured — so do not build a Q2.5 or M2 figure on them.

## Tests

```
Rscript tests/run_all.R                          # all twelve files (testthat::test_dir)
Rscript tests/test_values_map.R                  # one file: target names identical to the legacy formula
```

Files: `test_combine_taxo_assignments.R` (combine == sequential chain, empty-assignment guard), `test_values_map.R`, `test_cv_to_tidy.R`, `test_create_fake_pq_from_refseq.R`, `test_make_databases_helpers.R` (awk/sed helpers, `derive_fake_ref()`, `water_fill()` and `derive_fungi_rep()` on tiny fastas), `test_assign_compute.R` (derived lca and blastn rows equal the direct assignment; skipped without vsearch / BLAST+), `test_autometric_helpers.R`, `test_itsx.R` (refseq replacement; the real ITSx runs, ITS1 and full ITS, are skipped without the `itsxenv` env), `test_cv_queries.R` (query subsample rule, row alignment of the assignments; the cutadapt run is skipped without `cutadaptenv`), `test_otu_taxonomy.R` (OTU taxonomy from member ASVs; the vsearch run is skipped without vsearch), `test_mock_truth.R` (per-unit truth of a mock: name helpers, lineage cache, matching rule and ties; the vsearch runs are skipped without vsearch), `test_score_assignments.R` (which `external_scoring` each column is given, and the scoring call it wraps). Run from the project root; each file is self-contained. Add a test next to these for any helper under `R/` or `make_databases.R`.

## Other notes

- `tar_option_set(seed = targets_seed)` (22) in every pipeline; `derive_fake_ref()` uses the same seed. Preserve it.
- `n_threads = 4`, `n_workers = 3` in `config.R`; the parallel pipeline can use up to `n_workers * n_threads + n_threads` cores.
- `store_*/meta/meta` are tracked on purpose (decision 9); `tar_prune()` in `make.R` keeps them small. `figures/` is not tracked (decision 8).
- `Archives/` holds files no longer wired into any pipeline (sequential `script_assign_taxo.R`, `some_bash_script`, `archive_benchmarking_taxonomy.R`, `R/cutadapt_rm_primers_db.R`). Do not edit unless asked.
- `ROADMAP.md` is the todo: section 1 blocks the first production run, section 2 is the manuscript framing, section 4 the code debt of the 2026-09-22 critique, section 5 the remaining structural work of 2026-09-10 (`critique_fable/` holds its evidence). `HISTORY.md` is the archive, `CONTEXT.md` the glossary.
