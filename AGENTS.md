# AGENTS.md

Guidance for AI coding agents (and for any human arriving in this folder). It
describes what is on `make.R`'s path **today**; history lives in `Archives/` and
in [`HISTORY.md`](HISTORY.md). `CLAUDE.md` is a pointer to this file.

**Read order.** This file (how to work here) → [`CONTEXT.md`](CONTEXT.md)
(glossary) → [`ROADMAP.md`](ROADMAP.md) (the todo, one line per item) →
`docs/objectives_design.md` (the design, and the only place a design fact is
*defined*). `HISTORY.md` holds every ticked item, dated narrative and past
decision: read it for the measurements behind a ROADMAP line.

**Figures live in [`ETAT.md`](ETAT.md)** (targets, outdated and orphan
targets, computations, units per input, tests; newest first), written by
`tools/etat.R` at the end of every `run_bio_night.R`, never by hand. When a
count is needed in a document, point to `ETAT.md` instead of writing it.

**Review documents live in `reviews/`**; the latest is
[`reviews/code_review_2026-10-05.md`](reviews/code_review_2026-10-05.md).

**Two review documents of 2026-09-22**, both addressed but not applied:
[`reviews/critique_2026-09-22.md`](reviews/critique_2026-09-22.md) (aims and code, with the
developer's framing answer in its addendum) and
[`reviews/literature_2026-09-22.md`](reviews/literature_2026-09-22.md) (positioning against the
six benchmarks this project must cite). Their framing proposals live as `C22-…`
items in `docs/manuscript_framing.md`; their code-debt items were deleted on
2026-09-30 as stale, to be replaced by a new code review at the end of the
project (`ROADMAP.md` §4). **Do not apply a critique item on its own
authority**: the developer triages them.

## What this project is

A benchmark of taxonomic assignment methods (dada2 / sintax / lca / blastn) against nine reference databases (`config.R::benchmark_dbs`: UNITE with and without singletons and EUKARYOME ITS, each whole, Fungi-only and Fungi plus non-fungal representatives) for ITS fungal metabarcoding, evaluated on three mock communities (Pauvert et al. 2019, Tedersoo et al. 2026, Rolling et al. 2022), an in silico community (InSilicoSeq) and two biological datasets read through their controls only (GSSP air, Korhonen wood). Cross-validation is parked (OD 28) and its code archived (`Archives/cross_val/`). **It is not an R package** — it is a set of `{targets}` pipelines plus a Quarto analysis project. Do not run `devtools::check()`, `R CMD check`, or any R-package release tooling against this directory.

Everything is built around `phyloseq` objects and the sister pqverse packages (see *External code dependencies*).

**Scope (developer, 2026-09-17): taxonomic assignment of Fungi, and nothing else.** No community ecology: the miaSim / Hubbell branch is dropped, and biological datasets are read only through `NA_real`, `NA_fake` and `ext_fungi`.

## Design of 2026-09-17: what is built, what is not

The design was redrawn on 2026-09-17 and is **implemented**: items 0.1–0.10 are done (`HISTORY.md` §0) — the nine-database grid, the `fungi + rep` databases, the OTU input, the lca `lca_cutoff` and blastn `min_id` sweeps, the per-dataset ITSx region, then the per-unit truth of the mocks (`R/mock_truth.R`), the Hleap scoring in comparpq (`tc_metrics_unit()`), the chapters reading it and the documents. **Not built**: the cross-validation databases (0.11a, `docs/parked.md`), the final zoom (0.11b, ROADMAP) and the discussion notes (0.11c, `docs/manuscript_framing.md`). **Since OD 34 (2026-09-27) every project assigns its ITSx inputs only** (region × ASV / OTU, input names `<region>_<unit>`, e.g. `ITS1_asv`) against the nine databases, one computation at a time (`config.R::single_worker`). What has run and what is outdated: `ROADMAP.md` (status) and `ETAT.md` (counts). Before changing `config.R`, `R/values_map.R`, a pipeline or a chapter, read:

- **`docs/objectives_design.md`** — the source of truth for the design: inputs, controls, databases, methods, parameters and metrics per objective; its §4 numbers the developer's decisions 1–27 (not the same numbering as the decisions of `HISTORY.md` §3–3d: always write "objectives_design decision N"), and a new decision is written there first;
- **`docs/hleap_2021_metrics.md`** §5–6 — per-unit truth and Hleap-style scoring;
- **`HISTORY.md` §0** — the implementation plan, item by item, with what each one measured.

`docs/experimental_design.md` §3 holds costs measured on the grid of 2026-09-16, before OD 34: an order of magnitude only, until the M2.4 run (OD 31) replaces them.

In short: databases `[EUK_ITS_v2.1 | Unite_all_20250219 | Unite_s_all_20250219] × [all | fungi | fungi + rep]`, no `_cut` database in the grid (`config.R::zoom_dbs` keeps `EUK_ITS_v2.1_Fungi_cut` for a final zoom on the best combinations, Pauvert only); inputs = ITSx region(s) of the dataset × ASV / OTU (OD 34; regions and units per dataset in `config.R::bio_datasets`); the method rows come from `build_values_map()` (`R/values_map.R`: dada2 and sintax `min_bootstrap`, lca `lca_cutoff` ∈ {0.51, 0.8, 1} from one stored vsearch hit table, blastn votes × `min_id`), one computation per method × database × input; per-unit truth built from the Sanger (or expected) sequences of the mocks (`R/mock_truth.R`), one confusion-matrix cell per unit and rank as in Hleap et al. 2021 (`comparpq::tc_metrics_unit()`). Grid sizes per project: `ETAT.md`. Primary endpoints: F1 and PPV at **Genus** (OD 43); Species is read by exact and near match in its own section (chapter 09).

Chapter 01 scores every input against the per-unit truth: `unit` = `ASV` (`d_all_taxo_<region>_asv`, plus the consensus columns), `OTU` (`d_all_taxo_<region>_otu`, the OTUs assigned directly) and `OTU_derived` (`resolve_otu_taxonomy()` on the ASV assignments, kept so Q5 compares "assign then cluster" with "cluster then assign"). It writes one long `res_comp_tax` with `dataset`, `unit` and `external_scoring` columns, plus `res_abundance` (the M1.5 sweep). **Chapters 02 to 05 must read it through `default_view()`** (`00_setup.R`), which keeps one dataset, one input and one reading of the consensus columns — without it they average the three inputs and count every consensus twice, silently.

## Framing of the manuscript (2026-09-22)

Decided by the developer and recorded in `HISTORY.md` §3d (decisions 37–40) and
`reviews/literature_2026-09-22.md` §4. Useful whenever a change touches what the paper
claims:

- **Orsholm et al. 2026 is the foundation, not a competitor.** For fungal ITS the algorithm-family question is settled (composition-based classifiers win), so the four methods are a **deliberate zoom**, not a gap. Honest nuance: `lca` is search-based like `blastn`, so the zoom is 2 composition-based (dada2, sintax) × 2 similarity-search (lca, blastn), with blastn as the abstaining contrast.
- **The complexity goes where the literature has none**: the composition of the reference database (Q2) and voting across classifications (Q3). **Q3 is a primary axis**, not a supplementary one.
- **A vote is only as good as its voters are independent.** `Unite_all` ⊂ `Unite_s_all`, and lca and blastn read the same reference, so a naive consensus counts shared evidence twice (ROADMAP C22-P6, C22-P21). Never present a Q3 gain as a voting result without that check.
- **The harness and part of the benchmarked software share an author** (MiscMetabar, comparpq). The engines (vsearch, BLAST+, dada2) are third-party; the `session_info` targets record each checkout's commit and whether it was dirty (`pqverse_versions()`, `project_version()`); the paper declares it (ROADMAP C22-N1).


## Layout

```
_targets.yaml        # the named {targets} projects (TAR_PROJECT)
config.R             # shared constants: primers, threads, pinned search settings, ITSx, databases, datasets, mini_db
make.R               # sourcing it runs DB derivation + every dataset's projects (a day or more)
run_bio_night.R      # the night runner: per dataset, its DADA2 store then its grid; log in logs/
make_databases.R     # idempotent, atomic derivation of every reference-database variant (dbpq)
pipelines/           # dada2.R (Pauvert), dada2_bio.R and load_bio.R (other datasets), assign_taxo.R, species.R
R/                   # one file per concern; run_project.R (projects of a dataset), values_map.R (the grid),
                     # assign_compute.R, itsx_inputs.R, mock_truth.R, score_assignments.R, species_*.R,
                     # autometric_helpers.R + sample_tree.R (costs), load_pqverse.R, etat_helpers.R, ...
analysis/            # Quarto project: 00_setup.R, chapters 01–05 and 07–09, in_silico_simulation.qmd, sandbox/
tests/               # Rscript-runnable testthat files (tests/run_all.R)
tools/etat.R         # writes ETAT.md
tools/approach_schema.R  # docs/approach_schema.csv (relations) -> docs/approach_schema.svg + approach_schema_use<k>.svg (Graphviz dot)
tools/mock_scoring_schema.R  # docs/mock_scoring_{matrix,metrics,table}.csv -> docs/mock_scoring.svg (truth matrix, metrics, extended mock table)
figures/             # notebook output (gitignored)
data/                # see data/README.md
store_*/             # targets stores (meta/ tracked)
logs/                # run logs (gitignored)
docs/                # objectives_design.md (decided design, read first), hleap_2021_metrics.md (scoring and
                     # per-unit truth), experimental_design.md, reference_databases.md, manuscript_framing.md,
                     # parked.md, diagnostics_cv_and_databases.md
reviews/             # review documents (code review 2026-10-05, critiques and literature of 2026-09-22, ...)
ROADMAP.md HISTORY.md ETAT.md   # the todo, the archive, the dated figures
Archives/            # superseded code, kept for reference (cross_val/, sandbox_pre_od34/, ...)
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
| `dada2` | `pipelines/dada2.R` | `store_dada2` | DADA2 denoising of the Pauvert mock → `d_asv` |
| `assign_taxo` | `pipelines/assign_taxo.R` | `store_assign_taxo` | the grid on the Pauvert ITSx inputs → `d_all_taxo_<input>`, `benchmark_costs` |
| `assign_taxo_mini` | same script | `store_assign_taxo_mini` (gitignored) | smoke test on the `mini_*` databases |
| `dada2_<dataset>` | `pipelines/dada2_bio.R` (fastq) or `pipelines/load_bio.R` (published OTUs, OD 30) | `store_dada2_<dataset>` | `d_asv` of `tedersoo_illumina`, `insilico`, `rolling_mock`, `gssp_air`, `korhonen_wood` |
| `assign_taxo_<dataset>` | `pipelines/assign_taxo.R` | `store_assign_taxo_<dataset>` | the same grid on that dataset's ITSx inputs (regions and units from `config.R::bio_datasets`) |
| `species`, `species_mini` | `pipelines/species.R` | `store_species`, `store_species_mini` | species by exact and near match (OD 43), reads every assignment store: run it last |

The dataset projects come in pairs (`dada2_<dataset>`, `assign_taxo_<dataset>`, ROADMAP D1c); `config.R::bio_datasets` holds each dataset's fastq folder, sample file and primers, and `config.R` derives the dataset from `TAR_PROJECT`.

`make.R` runs `derive_all_variants()` then, for every dataset of `production_datasets()`, its DADA2 project and its grid through `run_project()` (`R/run_project.R`: `tar_make()`, then `tar_prune()` only when nothing is left outdated). The nights use `Rscript run_bio_night.R <dataset> ...` (same order per dataset, each project pruned after a build that left nothing outdated, log in `logs/`, `ETAT.md` refreshed at the end). To build one project, `source("R/run_project.R"); run_project("assign_taxo_mini")`: never `source("make.R")` for that, it runs everything. `config.R` sets `mini_db <- grepl("_mini$", TAR_PROJECT)` and prints its value: the `*_mini` projects use the `mini_*` fasta files and write to their own stores, so a smoke test can never overwrite a production store. Target names are the same in both.

Warn the user before launching: `derive_all_variants()` takes about 40 min when it rebuilds the `fungi + rep` databases (`force = TRUE`, mostly the vsearch pass that keeps the representatives away from the controls), a dataset's grid takes hours, and dada2 × full `EUK_ITS_v2.1` needs about 55 GB: the machine must be free. The DADA2 projects need cutadapt in the `cutadaptenv` conda env.

**dada2 is the patched build** (developer, 2026-09-30): Bioconductor dada2 1.40.0 scores every taxonomy beyond the 65 536th of a reference with the wrong k-mer profile (HISTORY.md §2). The project uses the build of `patches/install_dada2_patched.R` in `r_lib_patched/` (git-ignored), which `.Rprofile` puts first for any R session started at the project root; `assert_dada2_patched()` stops the pipelines otherwise. After a fresh clone, or after an update of R, run `Rscript patches/install_dada2_patched.R` first. A session started outside the project root, or with `--vanilla` / `--no-init-file`, gets the Bioconductor build.

Every store also has a `session_info` target (`pqverse_versions()`, `project_version()`, `sessioninfo::session_info()`) recording which commits produced the run and whether each checkout was dirty, with a hash of its diff (code review §1.7). Commit before a production run.

## Analysis (Quarto project in `analysis/`)

```
quarto render analysis                      # everything in _quarto.yml's render list
quarto render analysis/02_q1_methods.qmd    # one chapter
```

- `00_setup.R` — sourced first by every chapter: packages, `load_pqverse()`, `values_map`, constants (`tax_order`, `single_methods`, `preference_pattern`), helpers (`filter_default_settings()`, `save_fig()`, `cache_write()` / `cache_read()`, `plot_tc_metrics_mock()`).
- `01_load_and_clean.qmd` — reads `d_all_taxo` from `store_assign_taxo`, cleans taxonomy strings, builds `Gen_sp_*`, the OTU object `d_otu_taxo` (Q5: `resolve_otu_taxonomy()` with the `otu_clusters` target) and the consensus columns, scores the three inputs against the per-unit truth with `R/score_assignments.R` (which calls `comparpq::tc_metrics_unit()` once per `external_scoring` group), and writes `res_comp_tax`, `res_abundance`, `res_subsample` (Q5: ASVs subsampled to the OTU count, 100 draws × 2 schemes, about 15 of its ~21 min of render) `res_paired` (Q4: each region on the ASVs trimmed in every region of its dataset, `ASV_common`, and on every common ASV, `ASV_common_all`, OD 38) and `itsx_trimming` (trimmed / untrimmed units per region) to `analysis/_cache/*.rds`. **Render it first**; the other chapters read the cache.
- `02_q1_methods`, `03_q2_databases`, `04_q3_consensus`, `05_m2_costs` — one manuscript question each; figures go to `figures/` via `save_fig()`. Their summaries read Genus only (OD 43). Chapter 06 (cross-validation) is in `Archives/cross_val/`.
- `07_q4_its_region` — Q4: ITS1 vs ITS2 vs full on the datasets with several regions and a truth (Tedersoo, in silico), ASVs, one tab per dataset; reads both `res_comp_tax` and `res_insilico` (written by `in_silico_simulation.qmd`, which `_quarto.yml` therefore renders first).
- `08_q5_units` — Q5: ASV vs OTU vs OTU_derived per mock and ITSx region (the views of 02–05), single methods, default settings; reads `res_comp_tax`.
- `09_species_exact` — the species section (OD 43): the exact-match calls of the `species` store (`assignSpecies`, vsearch, blastn; `R/species_scoring.R` turns hits into calls with the unanimity rule) scored at Species with chapter 01's functions and truth, beside the main grid's `Gen_sp` rows; agreement of the three methods, precision and recall, why few species are named, controls, cost (S5), and the near matches of vsearch and blastn (S6: one difference, or 90 % covered at 100 % with the record stopping early). Reads `res_comp_tax` and the `species` store, one tab per mock input present in it.
- `in_silico_simulation.qmd` — D1b (InSilicoSeq): the input fasta drawn from `config.R::insilico_source`, the Docker `iss generate` recipe, the per-record truth parsed from the UNITE headers, and the scoring of `assign_taxo_insilico` (positive control, OD 29); writes the `res_insilico` cache chapter 07 reads.
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

**Rule for targets that run on crew workers** (`pipelines/assign_taxo.R`): workers only attach the CRAN/Bioconductor packages listed in `tar_option_set(packages = ...)` (`here`, `phyloseq`, `dplyr`, `tidyr`, `tibble`, …), never the pqverse checkouts. So a target either (a) calls `load_pqverse()` as its first statement (the `tar_eval` assignment and CV targets do), or (b) is declared `deployment = "main"` and runs in the main process where the checkouts are already loaded (every other target: store loading, negative controls, `tar_combine`, `benchmark_costs`, `session_info`). A "could not find function" error from `tar_make()` means a target broke this rule; if the missing function is from a CRAN package, add it to `packages`.

External tools: `vsearch` and BLAST+ on `PATH`, cutadapt in the `cutadaptenv` conda env (`config.R::cutadapt_conda_prelude`), ITSx 1.1.3 in the `itsxenv` conda env (`config.R::itsx_conda_prelude`; `conda create -n itsxenv -c conda-forge -c bioconda itsx`), InSilicoSeq for D1b-A.

## External code dependencies

- **`MiscMetabar`** — `add_new_taxonomy_pq`, `assign_*`, `cutadapt_remove_primers`, `filter_trim`, `chimera_removal_vs`, `mumu_pq`, `asv2otu`, `sam_data_matching_names`, `track_wkflow`, `add_dna_to_phyloseq`, `simplify_taxo`, `tc_bar`, `tc_circle`. `itsx_pq()` (with `add_detected = TRUE`, column `ITSx_detected`, since 2026-10-07). dada2 is an *Import* of MiscMetabar, so the DADA2 pipelines attach `dada2` themselves and `R/assign_compute.R` calls `dada2::assignTaxonomy()` (crew workers never attach it).
- **`comparpq`** — `tc_metrics_mock`, `tc_metrics_unit` / `tc_metrics_unit_vec` (per-unit scoring as in Hleap et al. 2021, added 2026-09-22 for item 0.7; **this is what the chapters call** — `tc_metrics_mock()` is kept for a mock that only comes with a list of expected taxa), `resolve_taxo_conflict`, `taxtab_replace_pattern_by_NA`, `rename_ranks_pq`, `add_shuffle_seq_pq`, `add_external_seq_pq` (the two negative-control injectors live here, not in MiscMetabar; checked 2026-09-16 in `comparpq/R/fake_creation.R`).
- **`dbpq`** — `cutadapt_rm_primers_db` (reverse-complements `primer_rev` by default since 2026-09-15), `filter_db`, `format2sintax`, `format2dada2`, `find_vsearch`, `is_vsearch_installed` (used by `make_databases.R`). Helpers still local to this project and worth upstreaming are listed in `reviews/proposals_for_dbpq.md`; read it before adding a new helper.
- **`greenAlgoR`** — `ga_footprint()` for the CO₂eq figures (M2.2).
- **`tidypq`** — being adopted (ROADMAP decision 13, items S5.1–S5.2): use `tidypq::pq_to_tidy()` instead of `psmelt()`, and `tax_table_to_df()` / `mutate_taxa_pq()` / `filter_taxa_pq()` / `rename_samples_pq()` instead of editing phyloseq slots by hand. Call-site inventory in `reviews/critique_fable/06_tidypq_adoption.md`.

## Pipeline architecture

### `pipelines/dada2.R` → `store_dada2`

Paired-end DADA2 flow for the ITS-1F/ITS2 primer pair (Pauvert et al. 2018), `n_threads` and `seq_len_min` from `config.R`. Stages: cutadapt primer removal (conda env), `filter_trim` → `derepFastq` → `learnErrors` → `dada` → `mergePairs` → `makeSequenceTable` → `chimera_removal_vs` → length filter → sample renaming (`sam_data_matching_names`, `samp_` prefix) → seed taxonomy with `assignTaxonomy` on `refseq/dada2_format/<seed_taxonomy_db>.fasta` (`config.R`; rank `Phylum`, not `Phyla`) → `d_asv`, `track_df`, `benchmark_costs_dada2`, `session_info`. Heavy phases are wrapped in `with_autometric()`. No OTU or mumu variant since 2026-10-07: the grid clusters the ITSx ASVs itself (OD 34). `pipelines/dada2_bio.R` is the same flow for the other fastq datasets (`config.R::bio_datasets`), `pipelines/load_bio.R` builds `d_asv` from published OTUs (OD 30).

### `pipelines/assign_taxo.R` → `store_assign_taxo`

- `d_asv_file` (`format = "file"` on `store_dada2/objects/d_asv`) → `d_asv`, so a rebuilt DADA2 store invalidates everything downstream.
- Negative controls (`R/itsx_inputs.R`): `fake_` shuffled from the trimmed units (`add_shuffle_seq_pq`) and `external_` from `data/data_raw/fake_ref/fake_ref_asv_100.fasta` (`add_external_seq_pq`; one record per phylum first, then drawn from the non-fungal kingdoms retained in every release, `config.R::rep_kingdom_min_records`, `external_control_kingdoms()`; redrawn by OD 39). The fasta is the file target `fake_ref_file`: regenerating it with `derive_fake_ref(force = TRUE)` invalidates every assignment (`derive_all_variants()` alone keeps an existing file; docs §6). Taxa matching `^fake_|^external_` are the controls; `external_scoring` decides whether the `external_` ones enter the confusion matrix (`_Fungi` databases) or are reported through `ext_fungi` / `ext_correct` (whole releases and `_Fungi_rep`).
- ITSx inputs, the only inputs since 2026-09-27 (OD 34): per region of `config.R::itsx_regions`, `d_<region>_asv_real` (`MiscMetabar::itsx_pq()` on `d_asv`: ITSx in the `itsxenv` conda env, all organism profiles because it only trims flanks, undetected ASVs keep their sequence, ASVs identical once trimmed merged into the most abundant, origin in `ITSx_origin`, detection in `ITSx_detected`) → `itsx_controls_<region>` (the external controls trimmed alone, `R/itsx_inputs.R`) → `d_<region>_asv_for_assignation` (`with_controls()`: `fake_` shuffled from the trimmed units under `withr::with_seed(targets_seed, …)`, then the trimmed `external_`, their origin written back). With an OTU input (`itsx_units`), `d_<region>_otu_real` (`asv2otu()` on the trimmed ASVs), `otu_clusters_<region>` and `d_<region>_otu_for_assignation`. One `d_all_taxo_<input>` per input.
- `values_map <- build_values_map(dbs = db_list, mini_db = mini_db)` (`R/values_map.R`, `build_methods_grid()`): the method/parameter rows × the nine databases × the inputs, one `<input>_<method>__<db>___<parameters>` row target per line but one computation per method × database × input (`<input>_compute_<method>__<db>`, ROADMAP S8.2, `R/assign_compute.R`). A computation runs `compute_assignment()` on one crew worker when `config.R::single_worker` is TRUE (the M2.4 setting: one computation at a time), with `config.R::assign_threads[[method]]` threads (all cores, sintax 1) and the pinned search settings (`lca_*`, `blastn_args`, `sintax_cmd_args`, `sintax_seed`), after the `dada2_build_on_worker` check of the patched dada2, inside `load_pqverse()` and `with_autometric(compute_name, …)` (which also writes the memory and CPU of the process tree, read by `add_tree_costs()` into `benchmark_costs`; logs under `data/data_final/autometric/<TAR_PROJECT>`); each row target then runs `derive_assignment()` in the main process (dada2 / sintax threshold; lca at its `lca_cutoff` from the stored vsearch hits through `MiscMetabar::assign_vsearch_lca(hits_table = )`, identical to vsearch `--lcaout` at that cutoff; blastn `min_id`, `min_cover` and vote through `MiscMetabar::assign_blastn(blast_table = )`). sintax and lca are called with `clean_pq = FALSE`, otherwise they skip the negative controls (B19). Each computation reads its reference through a file target (`ref_file_targets`: one `ref_sintax__<db>` / `ref_dada2__<db>` per fasta). Sizes per project: `ETAT.md`.
- Post-clustering from member ASVs (Q5, decision 25), kept beside the directly assigned OTUs so Q5 compares "assign then cluster" with "cluster then assign" (`unit = "OTU_derived"`): `otu_clusters_<region>` (`otu_membership()` in `R/otu_taxonomy.R`: 97 % vsearch clusters of the trimmed ASVs, with the most abundant member as archetype); chapter 01 merges `d_all_taxo_<region>_asv` into OTUs and keeps, per assignment column, the value all member ASVs agree on (`resolve_otu_taxonomy()`, `unanimity`, NA ignored).
- **`full_name` strings are the target names** (`ITS1_asv_dada2__Unite_s_all_20250219___0.5`, `ITS1_otu_lca__EUK_ITS_v2.1___0.8`, `ITS2_asv_blastn__EUK_ITS_v2.1___0.5...rel_majority...100...95`, …). `tests/test_values_map.R` pins them; changing `full_name_for()` orphans every store object.
- `d_all_taxo_<input>` = `tar_combine` of the assignments of one input through `combine_taxo_assignments(d_<input>_for_assignation, !!!.x)` (`R/combine_taxo_assignments.R`, detects new columns by `setdiff()` so it does not need to know suffixes).
- `benchmark_costs` aggregates the newest autometric log per computation (`R/autometric_helpers.R`: `read_autometric_dir()` + `summarise_autometric_costs()`) and joins `values_map` on `compute_name`: the rows of one computation share its cost, and `phase` holds the `full_name` the chapters join on.

### `pipelines/species.R` → `store_species`

Species by exact and near match (OD 43, chapter 09): the inputs of every assignment store pooled into distinct queries (`pool_queries()`), then per database `dada2::assignSpecies()` on `species_format/`, vsearch exact (`--id 1 --query_cov 1`, kept when the query is an exact substring), vsearch near candidates (`--id 0.9 --query_cov 0.9`) and blastn's own search (megablast, up to 2 000 records per query, records indexed by position because BLAST+ rewrites some ids). Near hits carry `diffs` (differences over the whole query) and `exact_cover` (100 % over ≥ 90 % of the query, the record stopping early), computed from the sequences the same way for both methods (`R/species_exact.R`); the call rule is a scoring choice (`R/species_scoring.R`). `SPECIES_DATASETS` restricts the datasets. Costs logged like the grid's.

### Cross-validation (archived)

Parked by OD 28; code, chapter and night script in `Archives/cross_val/` (its README lists the known defects); `store_cross_val` untouched at the root.

## Reference databases

Full procedure, header formats and how to move to a new release: **`docs/reference_databases.md`**. Summary:

- `config.R::reference_sources` lists the general FASTA releases with versioned names (`Unite_s_all_20250219` = UNITE 19.02.2025 eukaryotes with singletons, `Unite_all_20250219` = without singletons, `EUK_ITS_v2.1`), their URL and DOI. The SSU release was dropped on 2026-09-15 (decision 7 revised: ITS1 is not part of the 18S gene). `config.R::benchmark_dbs` lists the benchmarked databases (source × simplification: full, `_Fungi`, `_Fungi_rep`), `config.R::zoom_dbs` the databases built for the final zoom only (`EUK_ITS_v2.1_Fungi_cut`); `db_list`, `cv_db_list`, `db_meta`, `seed_taxonomy_db`, `fake_ref_source` and `preference_db` come from there. **Never hard-code a database name elsewhere**: pipelines, chapters and tests read `config.R`.
- `make_databases.R` step 1, `download_reference_sources()`: `dbpq::download_unite_db(url, extract = TRUE)` / `dbpq::download_eukaryome_db(url)` → `data/data_raw/refseq/sources/<source>.fasta` + `sources/manifest.csv` (URL, DOI, md5, date).
- Step 2, `derive_all_variants()`: `derive_dada2()` / `derive_sintax()` → `dada2_format/<source>.fasta` and `sintax_format/<source>.fasta` (both go through `general_headers_fixed()`: the UNITE `|k__` separator that makes dbpq drop the kingdom, and EUKARYOME name qualifiers such as `Lactarius(Fungi)`, `(Candida)` or `Mortierella.s.str` that break sintax bootstraps and mock matching), `derive_kingdom_only()` (kingdom exactly Fungi, anchored pattern), `derive_cutadapted()` (trimmed at ITS1F and at the reverse complement of ITS2 when found, records without primer sites kept, records shorter than `cut_min_length` after trimming dropped: `config.R::cut_discard_untrimmed`, `primer_min_overlap`, `cut_min_length`; ROADMAP S6.9), fake reference, `derive_fungi_rep()` (`_Fungi` records + `config.R::rep_share` × the fungal record count of non-fungal records, rounded up, drawn from the full release of the same source over the kingdoms of `retained_kingdoms()` by `water_fill()`, never a record within `config.R::rep_max_identity_to_external` (0.97) of an `external_` control, excluded from the candidates by one vsearch pass, candidates as queries and controls as database, which costs about 35 min for the three databases; same draw written in both formats; list in `sources/<db>_reps.csv`; its `mini_*` files hold 10 000 Fungi records + 500 representatives), `mini_*` subsets, `derive_species()`. Idempotent and atomic (each file written to `<file>.part.fasta`, then renamed); an output older than its input is reported (`STALE` warning), not rebuilt, since the stores hash these files; `force = TRUE` rebuilds. Each reference fasta is a file target (`ref_dada2__<db>` / `ref_sintax__<db>`, column `ref_file` of `values_map`), so a database rebuilt in place reruns exactly the targets that read it (ROADMAP S7.9).
- dada2 reads `dada2_format/`, sintax/lca/blastn read `sintax_format/` (`db_path_for()` in `R/values_map.R`). `mini_*` files are used by the `*_mini` projects only (decision 4: never in the Results).

**Legacy files deleted (2026-09-11).** The unversioned reference files (`data/data_raw/refseq/{Unite,Unite_RefS,Euk_ITS_v2,Euk_SSU_v2}.fasta` and the unversioned files of `dada2_format/` / `sintax_format/`) lost the UNITE kingdom in sintax headers (`…|k__Kingdom;tax=p:…`, ranks shifted by one) and let non-Fungi records into the `_Fungi` files; do not recreate them. Every store built between 2026-05-26 and 2026-09-10 also ran all targets on the first database (`ifelse()` recycling in `values_map`, fixed). See `reviews/critique_fable/07_results_validity.md`.

## Cost layer (M2)

`autometric::log_start()` only logs the process it runs in, so logging starts *inside* the crew worker: `with_autometric(phase, expr, dir)` writes `data/data_final/autometric/<TAR_PROJECT>/<phase>__<timestamp>.txt` per target and per run, plus, since 2026-10-07, the RSS of the worker and its child processes (vsearch, BLAST+; `R/sample_tree.R`, `__tree.csv`) and their CPU time (`__cpu.csv`). `summarise_autometric_costs()` keeps the newest log per phase; `add_tree_costs()` adds `wall_s`, `cpu_self_s`, `cpu_child_s` and `peak_tree_mb`. Costs logged before M2.4 are provisional (OD 31): read them as orders of magnitude. blastn uses about one core whatever its threads: `makeblastdb` (single-threaded) dominates each computation.

## Tests

```
Rscript tests/run_all.R                          # every file (testthat::test_dir); counts in ETAT.md
Rscript tests/test_values_map.R                  # one file: target names identical to the legacy formula
```

One `tests/test_<file>.R` per helper file of `R/` (and `test_make_databases_helpers.R`); the tests that need vsearch, BLAST+, cutadapt or ITSx skip without them. `test_assign_compute.R` checks that every derived row equals the direct call (lca, blastn, sintax) and pins dada2's left-packing rule; `test_values_map.R` pins the target names. Run from the project root; each file is self-contained. Add a test next to these for any helper under `R/` or `make_databases.R`.

## Other notes

- `tar_option_set(seed = targets_seed)` (22) in every pipeline; `derive_fake_ref()` uses the same seed. Preserve it.
- `config.R::single_worker` (TRUE) runs one computation at a time with `assign_threads[[method]]` threads; `n_workers` only applies when it is FALSE. `n_threads` serves the DADA2 pipelines and ITSx.
- `store_*/meta/meta` are tracked on purpose (decision 9); `run_project(prune = TRUE)` keeps them small; the night runs prune too, so an object that is no longer a target is gone after a clean night: read or save it before. `figures/` is not tracked (decision 8).
- `Archives/` holds files no longer wired into any pipeline (the cross-validation in `cross_val/`, sandbox scripts of the pre-OD 34 inputs in `sandbox_pre_od34/`, older scripts). Do not edit unless asked.
- `ROADMAP.md` is the todo, one short line per item, in the order the work happens (runs, analyses, developer decisions, end of project). `docs/manuscript_framing.md` holds the manuscript-framing proposals, `docs/parked.md` the unplanned items (cross-validation, upstream reports, side questions). `HISTORY.md` is the archive (§6: everything closed when the ROADMAP was simplified on 2026-09-30), `CONTEXT.md` the glossary.
