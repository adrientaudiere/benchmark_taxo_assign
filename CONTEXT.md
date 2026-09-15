# Context — benchmar_assign_taxo

Ubiquitous language for this project: the shared vocabulary used by the
developer and AI agents. Each entry maps a domain term or piece of jargon to a
plain definition, so conversations and code stay consistent and concise.

Idea from https://github.com/mattpocock/skills.

## Glossary

**benchmark** — Systematic comparison of four taxonomic-assignment methods (dada2, sintax, lca, blastn) across seven reference databases (`db_list` in `config.R`, from `benchmark_dbs`) for ITS fungal metabarcoding. Performance is evaluated on a mock community, cross-validation, in silico simulations, and a real biological community (Bokulich 2020 four-dataset approach).

**assignment method** — Algorithm used to assign a taxonomic label to an ASV/OTU sequence. The four methods in scope: `dada2` (RDP naive Bayes classifier, `assignTaxonomy`), `sintax` (vsearch SINTAX bootstrap), `lca` (Lowest Common Ancestor), `blastn` (BLAST nucleotide search with a vote algorithm).

**reference database** — FASTA file providing labelled sequences against which query ASVs are classified. Built by `make_databases.R` from a *reference source*, in dada2 and sintax header formats, as a full, `_Fungi` (kingdom exactly Fungi) or `_Fungi_cut` variant (then trimmed with cutadapt at ITS1F and at the reverse complement of ITS2 when found; records without primer sites are kept whole, decision 23). The benchmarked set is `config.R::benchmark_dbs`. Procedure: `docs/reference_databases.md`.

**reference source** — A general FASTA release downloaded as-is (`dbpq::download_unite_db()` / `dbpq::download_eukaryome_db()`) into `data/data_raw/refseq/sources/<source>.fasta`, listed in `config.R::reference_sources` with its URL and DOI. Source names carry the release (`Unite_s_all_20250219` = UNITE 19.02.2025 eukaryotes with singletons, `Unite_all_20250219` = without singletons, `EUK_ITS_v2.1`), so a new release yields new file and target names. Provenance in `sources/manifest.csv`.

**general FASTA** — The provider's reference release with taxonomy in `k__`/`p__` form (UNITE: `Name|Acc|SH|type|k__…`; EUKARYOME: `ID;k__…`), before conversion by `derive_dada2()` / `derive_sintax()`.

**EUKARYOME name qualifiers** — Annotations EUKARYOME adds to rank values: the kingdom or order of a homonym (`g__Lactarius(Fungi)`), bracketed genera (`g__(Candida)`), `.s.str` (sensu stricto), `.nom.prov` (provisional name). `make_databases.R::general_headers_fixed()` keeps the name only before conversion, so names match the mock and sintax bootstraps parse; placeholders (`Xxx.gen05`, `.fam.incertae.sedis`) are kept.

**`values_map`** — Tibble built by `build_values_map()` (`R/values_map.R`) from `db_list`, enumerating all (method, parameter, db) combinations to benchmark: 12 method/parameter rows × (`length(db_list)` databases with the raw ASVs + `length(itsx_db_list)` Fungi-filtered databases with the ITSx input; `preprocess` = `none` or `itsx`), i.e. 132 rows with the seven databases of 2026-09-15 (ROADMAP S8.5, S6.9). `input_pq` names the phyloseq target an assignment reads; ITSx rows have `itsx_`-prefixed `full_name`s. `compute_name` groups the rows computed together: one computation per method × database × input (`compute_<method>__<db>`, 44 in total) serves three rows, whose thresholds, votes or copies are derived afterwards (`R/assign_compute.R`, ROADMAP S8.2). Each row carries `full_name` (the targets name, pinned by `tests/test_values_map.R`), `db_path` (format chosen by method, `mini_*` when `mini_db`), `ref_file` (name of the file target tracking `db_path`: `ref_sintax__<db>` / `ref_dada2__<db>`), and the crew `controller`. `pipelines/assign_taxo.R` dispatches one target per row; `d_all_taxo` is assembled by `tarchetypes::tar_combine`. `build_cv_values_map()` is the cross-validation counterpart.

**`_targets.yaml` / `TAR_PROJECT`** — Five named targets projects: `dada2`, `assign_taxo`, `cross_val` (production stores) and `assign_taxo_mini`, `cross_val_mini` (same scripts on the `mini_*` databases, own gitignored stores). `Sys.setenv(TAR_PROJECT = "assign_taxo"); tar_make()`. `run_project()` (`R/run_project.R`) wraps `tar_make()` + `tar_prune()`; sourcing `make.R` runs every step, so source `R/run_project.R` to build a single project.

**`load_pqverse()`** — `R/load_pqverse.R`; loads the sister packages from their development checkouts under `PQVERSE_PKG_DIR` (decision 10). Called by every pipeline, inside every crew-worker target, and by `analysis/00_setup.R`. `pqverse_versions()` records version + git commit in the `session_info` target of each store.

**analysis chapters** — Quarto project in `analysis/`: `00_setup.R` (shared setup and helpers), `01_load_and_clean` (writes `analysis/_cache/*.rds`), `02_q1_methods`, `03_q2_databases`, `04_q3_consensus`, `05_m2_costs`, `06_cross_validation`, `in_silico_simulation.qmd`, `sandbox/`. Each chapter starts with `here::i_am()` because `_quarto.yml` is a `here` root criterion.

**tidypq** — pqverse package of tidyverse-style verbs for phyloseq objects (samples, taxa, occurrences, tree scales). Adopted in this project (decision 13): `tax_table_to_df()` to read a tax_table as a tibble, `mutate_taxa_pq()` to write columns back, `pq_to_tidy()` instead of `psmelt()` (unassigned ranks become `"Unknown"`, zero abundances are dropped). Phyloseq slots (`@tax_table`, `@otu_table`) are no longer edited by hand where a verb exists.

**ITSx input (`d_asv_itsx`)** — Copy of `d_asv` whose sequences are the ITS1 extracted by ITSx (`R/itsx.R`, `itsxenv` conda env), without the 18S flank and the 5.8S start of the ITS1F–ITS2 ASVs; ASVs ITSx does not detect keep their full sequence. ITSx only trims the flanks here: all organism profiles are used and its ability to keep Fungi only is deliberately not used. It receives the same negative controls as `d_asv` and feeds the `itsx_<method>` assignments (question Q4, decision 19).

**mock community** — Lab-made community of known fungal taxa (Pauvert et al. 2019) used as the primary ground truth for computing TP/FP/FN/TN/MCC/F1/TAR/TDR. The known taxonomy is in `data/data_raw/metadata/taxo_mock.csv`.

**fake taxa** — Negative-control sequences injected before assignment: shuffled ASV sequences (`add_shuffle_seq_pq`, prefix `fake_`) and sequences from an external FASTA (`add_external_seq_pq`, prefix `external_`): `fake_ref_asv_100.fasta`, 100 non-Fungi records of `config.R::fake_ref_source` drawn by `derive_fake_ref()`, tracked by the file target `fake_ref_file` so a regenerated file invalidates every assignment. Taxa matching `^fake_|^external_` are the TN denominator in `tc_metrics_mock()`.

**`tc_metrics_mock()`** — Function from `comparpq` that computes classification performance metrics (TP, FP, FN, TN, MCC, ACC, F1, PPV/TAR, TPR/TDR) for each (method × db × taxonomic rank) combination against a mock truth table.

**consensus voting / `resolve_taxo_conflict()`** — `comparpq` function that combines taxonomy assignments from multiple methods/databases into a single consensus call. Five strategies: `unanimity`, `consensus`, `abs_majority`, `rel_majority`, `preference`. Controlled by `strict` (boolean) and `nb_agree_threshold` (integer: 1, 2, or 3).

**cross-validation (CV)** — k-fold held-out evaluation of each (method, db) pair using the reference database itself as the source of truth. Implemented in `R/cross_val.R` (`cross_val()` and `cross_val_param()`), run via `pipelines/cross_val.R` → `store_cross_val`. Key parameters in `config.R`: `cv_fold_number`, `cv_fold_tested`, `cv_max_seq` (10 folds and 5000 sequences in production, 2 folds and 200 sequences when `mini_db`), `cv_oversample`.

**trimmed CV queries** — In cross-validation, the test sequences are cut to the ITS1F–ITS2 amplicon with cutadapt, like the reads of an Illumina run (`R/cv_queries.R`, decision 24); the training part keeps full-length records. Records without the ITS2 site are never queried but stay in the training part; `cv_oversample` records are drawn per wanted query so that each database gets `cv_max_seq` queries. The queries of a `_Fungi_cut` database, whose records no longer hold the primer sites, are drawn and trimmed from its `_Fungi` source (`query_fasta`) and paired by name with the `_cut` records. Records are named by their sintax-format identifiers for every method, also for dada2, whose headers are the taxonomy only (`id_fasta`, ROADMAP B22).

**OTU taxonomy (`d_otu_taxo`)** — Q5 object built in chapter 01: `d_all_taxo` merged into the 97 % vsearch clusters of the ASVs (`otu_clusters` target, `otu_membership()`), each assignment column keeping the value all member ASVs assigned at that rank agree on (`resolve_otu_taxonomy()`, `unanimity`, NA ignored). No assignment is rerun on OTUs; the negative controls are not clustered (decision 25). Its metrics are `res_comp_tax_otu`.

**`remove_tested`** — Cross-validation parameter (`remove_tested_sequences = TRUE/FALSE`) that controls whether the sequences being classified are excluded from the reference database (TRUE = standard, non-leaked; FALSE = leaked, upper-bound estimate).

**in silico simulation** — Two sub-branches (D1b): InSilicoSeq generates simulated fastqs from a curated FASTA (sequencing errors included); miaSim generates community matrices from neutral/niche ecological models (no sequencing errors). Both yield a known ground truth for metric computation.

**`store_dada2` / `store_assign_taxo` / `store_cross_val`** — Three independent `{targets}` stores at the project root, one per production project of `_targets.yaml`; `make.R` runs all three in order. `store_assign_taxo` depends on `store_dada2/objects/d_asv` as a file target. The notebooks read them via `tar_read(..., store = here("store_assign_taxo"))`. `meta/meta` is tracked in git (decision 9), `objects/` is not.

**`autometric` / `with_autometric()`** — Per-target logging of wall time, peak RAM and mean CPU. Because `autometric::log_start()` only logs its own process, `with_autometric()` (`R/autometric_helpers.R`) starts the logger *inside* the crew worker and writes one file per target and per run under `data/data_final/autometric/<pipeline>/`. `summarise_autometric_costs()` keeps the newest file per phase; `benchmark_costs` (assignments) and `benchmark_costs_dada2` (denoising) are the aggregated targets.

**`config.R`** — Single shared-constants file sourced by the three pipelines, `make_databases.R` and `analysis/00_setup.R`. Contains primers (`fw_primer_sequences`, `rev_primer_sequences`), `n_threads`, `seq_len_min`, `prop_fake`, paths, conda prelude for cutadapt, the primer-trimming settings (`primer_min_overlap`, `cut_discard_untrimmed`), `n_workers`, `targets_seed`, the CV knobs, `db_list`, `db_meta`, and `mini_db` (derived from `TAR_PROJECT`).

**`combine_taxo_assignments()`** — Helper in `R/combine_taxo_assignments.R` used by `tarchetypes::tar_combine` to merge per-(method, db) `tax_table` columns into one wide phyloseq object. Detects new columns via `setdiff(colnames(new_tt), colnames(base_tt))`; does not need to know per-row suffixes.

**`proposals_for_dbpq.md`** — Inventory of six helper functions still implemented locally in `make_databases.R` / `R/create_fake_pq_from_refseq.R` that fit `dbpq`'s scope and are candidates for upstreaming. Consult it before adding a new helper to avoid duplication.

**M1 metric set** — Primary performance metrics used in results figures: F1 + MCC (Hleap 2021) + TAR/TDR per taxonomic rank. TP/FP/FN/TN go to supplementary material only. Decided 2026-05-19.

**M2 environmental metrics** — Carbon-footprint and compute-cost layer: `wall_time_s`, `peak_resident_mb`, `mean_cpu_pct` from `benchmark_costs`, converted to CO₂eq via `greenAlgoR`.

**`cutadaptenv`** — Conda environment name containing `cutadapt`. Required by the DADA2 pipeline's primer-removal step, the `_Fungi_cut` database derivation and the trimming of the CV queries; the scripts prepend `source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv &&` before each cutadapt call.

## Key decisions

**Four assignment methods, IdTaxa excluded (2026-05-19)** — The methods set is fixed at dada2, sintax, lca, blastn. IdTaxa is excluded because it requires a separate training step, making cost comparisons unfair. The `idtaxa` block remains commented out in `values_map`; the manuscript will include a one-line justification pointing to Murali et al. 2018.

**No SSU database (2026-09-15, revises 2026-05-19)** — The 2026-05-19 decision kept only `EUK_SSU_v2_Fungi_cut`, on the belief that ITS lies inside the 18S gene. It does not: trimming an SSU record at ITS1F leaves the 46 bp of 18S that start the ASVs, so the SSU release is no longer benchmarked (ROADMAP decision 7, S6.9).

**Mini databases are smoke-test only (2026-05-19)** — `derive_mini()` in `make_databases.R` produces fast 10 000-line database subsets for iteration. They are excluded from all Results figures; use `dplyr::filter(!startsWith(db, "mini_"))` before any publication figure.

**Parallel pipeline over sequential** — `pipelines/assign_taxo.R` replaces the archived `Archives/script_assign_taxo.R`. Each (method, db) reads `d_asv_for_assignation` directly and runs in parallel via `crew_controller_local`; the sequential predecessor had each target depend on `previous_target`, creating a bottleneck.

**Biological validation dataset: Taudière et al. 2018 (2026-05-19)** — Tree endophyte ITS data from Taudière et al. (doi:10.1016/j.funeco.2018.07.008), with `site` and `height` as sample modalities. No objective ground truth exists; metrics are agreement-only (Jaccard, Bray-Curtis, richness). Pending: fastq availability check and SRA/Dryad/ENA fetch.

**Structural decisions (2026-09-10)** — Figures are not versioned (`figures/` gitignored); `store_*/meta/meta` stay tracked; the reference code is the dev checkout of every pqverse package (`load_pqverse()`); SRR30413326 joins as a second mock (D1d); tidypq is adopted in the notebook and the pipelines (ROADMAP S5.x); smoke tests run in the `*_mini` targets projects, never in a production store. Details in ROADMAP section 3b.

**blastn query cover 80 % (2026-09-11)** — `config.R::blastn_min_cover`. BLAST `qcovs` is the share of the ASV covered by the alignment with a reference, whatever the reference length. The ITS1F–ITS2 ASVs start with 30–45 bp of 18S that ITS references lack, so the `assign_blastn()` default of 95 rejected most true hits (ROADMAP decision 18, critique 07 B18). ITS extraction with ITSx is the proposed alternative (Q4).

**`tar_option_set(seed = 22)`** — Fixed random seed set in both `pipelines/dada2.R` and `pipelines/assign_taxo.R`. Must be preserved when modifying targets so pipeline re-runs stay reproducible.
