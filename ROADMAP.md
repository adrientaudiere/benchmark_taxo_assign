# Benchmark plan — answering the manuscript questions

Working plan derived from `Coffre_principal/PROJETS/Projets IdEst/Dev packages R/MiscMetabar/autour_de_MiscMetabar/benchmark_taxo_assign/` (`Taxonomic assignation manuscript.md` + `Taxonomic assignation.md`). Use this as the running todo; tick items as you go.

## Little side questions to explore

### Software parameters

What is the effect of set --maxaccepts 16 in sintax ?
What is the effect of set --dbmask none in sintax ?

### Databases building

What is the effect of taxing the Unite database with singletons as RefSeq only (Unite_RefS.fasta) ? 
What is the effect of adding only some sequences to a FUNGI filtered database ?

## Actions requiring manual execution (in temporal order)

These steps require you to actually run commands in a terminal or edit external files. They are ordered by the suggested execution sequence; later steps may depend on earlier ones completing.

### ✅ Phase 1 — Complete the mock-community Results section

| #   | Action                                                                                                                    | Command / location                                                                              | Blocker for                              | Status   |
| --- | ------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------- | ---------------------------------------- | -------- |
| 1   | Re-run assign_taxo pipeline with v2 databases (renamed from v1_9_3 in Q2.3)                                               | `tar_make(script="script_assign_taxo_parallel.R", store="store_assign_taxo")` from project root | Q2 figures with Fungi-only ITS databases | ✅ Done  |
| 2   | Add IdTaxa exclusion note to the methods table in `Taxonomic assignation.md` (external Markdown file in Coffre_principal) | Edit manually — mark IdTaxa as "evaluated elsewhere, not in this study"                         | Q1.5b                                    | ✅ Done  |

### ✅ Phase 2 — Cross-validation

| #   | Action                                                                                                                      | Command / location                                                                                                               | Blocker for                               | Status  |
| --- | --------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------- | ------- |
| 3   | **(Smoke test first)** Run CV pipeline with `cv_fold_tested = 2L` (already set in `config.R`) to verify the pipeline wiring | `tar_make(script="script_cross_val.R", store="store_cross_val")`                                                                 | Confirming D1a wiring before the long run | ✅ Done |
| 4   | Set `cv_fold_tested <- cv_fold_number` in `config.R` for the publication run, then re-run                                   | Edit `config.R`, then `tar_make(script="script_cross_val.R", store="store_cross_val")` (hours-scale)                            | D1a CV figures                            | ✅ Done |


### ✅ Phase 3 — Database derivation verification

| #   | Action                                                                                        | Command / location                                                                                              | Blocker for | Status  |
| --- | --------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------- | ----------- | ------- |
| 5   | Verify `make_databases.R::derive_all_variants()` end-to-end with new v2 source files         | `source("make_databases.R"); derive_all_variants()` — sources Unite.fasta, Euk_ITS_v2.fasta, Euk_SSU_v2.fasta  | C1 sign-off | ✅ Done |

### 🔄 Phase 4 — In silico simulations

| #   | Action                                                                                        | Command / location                                                                                                                                                     | Blocker for           |
| --- | --------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- |
| 6   | Curate input taxon list from UNITE Fungi (50–200 species, one ref seq each, spanning 7 phyla) | Manual selection from `data/data_raw/refseq/`                                                                                                                          | D1b-A InSilicoSeq run |
| 7   | Run InSilicoSeq to generate simulated fastqs (3 replicate seeds)                              | `iss generate --genomes mini_Unite.fasta --sequence_type amplicon --n_reads 10000 --abundance zero_inflated_lognormal --model MiSeq` (see `Archives/some_bash_script`) | D1b-A pipeline        |
| 8   | Run DADA2 + assignment on InSilicoSeq fastqs                                                  | `tar_make(script="script_dada2.R", store="store_dada2_insilico")` then `tar_make(script="script_assign_taxo_parallel.R", store="store_assign_taxo_insilico")`          | D1b-A analysis        |

### Phase 5 — Biological community (Taudière 2018)

| #   | Action                                                                                                                                   | Command / location                                                                                                                                              | Blocker for   |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| 9   | Check data availability in the paper (doi:10.1016/j.funeco.2018.07.008) — locate SRA / Dryad / ENA accession for raw fastqs or OTU table | Manual check of paper's data availability statement                                                                                                             | D1c setup     |
| 10  | Fetch raw fastqs (if available) into `data/data_raw/rawseq_endophyte/`                                                                   | `fastq-dump` / `wget` from accession                                                                                                                            | D1c DADA2 run |
| 11  | Confirm primers match ITS-1F/ITS2; if not, update `config.R` for this run only                                                           | Manual primer check against paper's Methods                                                                                                                     | D1c DADA2 run |
| 12  | Run DADA2 + assignment for biological community                                                                                          | `tar_make(script="script_dada2.R", store="store_dada2_endophyte")` then `tar_make(script="script_assign_taxo_parallel.R", store="store_assign_taxo_endophyte")` | D1c analysis  |

---

## Main research questions (from the manuscript Results outline)

1. **Q1** — Effect of classification algorithm and its parameters
2. **Q2** — Effect of reference database and database simplification (Fungi-only, `_cut`, clustering)
3. **Q3** — Effect of consensus voting across methods/databases/parameters

Plus the methodological scaffolding the manuscript announces:

- **D1** — Datasets: cross-validation, in silico, mock community, biological community (Bokulich 2020 four-data approach)
- **M1** — Performance metrics: TP/FP/FN/TN/MCC/ACC/F1/Precision/Recall (Hleap 2021 set)
- **M2** — Environmental metrics: CPU / wall time / CO₂eq

## Snapshot — what is already produced by the pipeline

| Piece                                                                    | Where                                                                              | Status                                                            |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| ASV phyloseq `d_asv`                                                     | `store_dada2`, from `script_dada2.R`                                               | ✅ runs end-to-end                                                 |
| Fake taxa injection (shuffle + external) → `d_asv_for_assignation`       | `script_assign_taxo_parallel.R`                                                    | ✅ TN material is in place                                         |
| 4 methods × 7 DBs × bootstraps → `d_all_taxo` (combined `tax_table`)     | `store_assign_taxo`, via `tar_combine` + `combine_taxo_assignments`                | ✅ produced                                                        |
| Per-target runtime / memory (`benchmark_costs`)                          | `store_assign_taxo`                                                                | ✅ wired, not yet plotted                                          |
| Mock truth table `taxo_mock`                                             | `data/data_raw/metadata/taxo_mock.csv`                                             | ✅                                                                 |
| NA cleanup + `Gen_sp` construction                                       | `analysis/benchmark.qmd` § "Import value from store_assign_taxo"                   | ✅                                                                 |
| `tc_metrics_mock()` (TP/FP/FN/MCC/ACC/F1) per (method × db × rank)       | `benchmark.qmd` via `comparpq`                                                     | ✅ produced as `res_comp_tax`                                      |
| Three consensus strategies (unanimity, rel_majority, preference) applied | `benchmark.qmd` § "Create column using consensus..."                               | ✅ partial                                                         |
| Per-rank NA proportion plots                                             | `benchmark.qmd` § "Proportion of NA"                                               | ✅ exploratory plots                                               |
| Cross-validation helper `cross_val()` + targets pipeline                 | `R/cross_val.R` + `script_cross_val.R` → `store_cross_val`                        | ✅ pipeline run (phases 2–3 done)                                  |
| In silico notebook                                                       | `In_silico_simulation.qmd`                                                         | 🔄 D1b-B (miaSim) in progress                                     |
| Biological community dataset                                             | Taudière et al. 2018, doi:10.1016/j.funeco.2018.07.008 (endophytes; site × height) | ⚠️ chosen, fastq fetch pending                                    |
| `benchmark.qmd` ↔ parallel store wiring                                  | `benchmark.qmd` still grep's `_all_taxo`                                           | ❌ needs the three small edits (see CLAUDE.md "Parallel pipeline") |

## Q1 — Algorithm & parameters

Already producible from `res_comp_tax`. What's missing is the **systematic figure set** and a stats layer.

- [x] **Q1.1** — Re-point `benchmark.qmd` at `tar_read(d_all_taxo, store = "store_assign_taxo")` (drop the `_all_taxo` grep, drop the `rename_ranks_pq(...)` block). Required before anything else in the qmd works on the parallel store.
- [x] **Q1.2** — Lock the **per-rank × method** figure: ACC, F1, MCC across `c("dada2", "sintax", "lca", "blastn")` faceted by rank K→S. `plot_tc_metrics_mock()` already builds this; pin a final version and save as `figures/fig_q1_methods.{pdf,png}`.
- [x] **Q1.3** — **min_bootstrap sweep figure** for dada2 / sintax: F1 or MCC as a function of `min_bootstrap ∈ {0.4, 0.5, 0.6}`, per rank. Aggregate from `res_comp_tax` (`bootstrap` column already extracted there).
- [x] **Q1.4** — **Vote-algorithm sweep figure** for blastn: F1 or MCC across `c("rel_majority", "abs_majority", "unanimity")` × DB × rank. The data is already in `res_comp_tax`; needs a dedicated faceted panel.
- [x] **Q1.5** — ~~Add IdTaxa~~ **Decided out** (decision 1). Action: keep the `idtaxa` block commented out in `script_assign_taxo_parallel.R::values_map`, and add a one-paragraph justification in the Methods section ("IdTaxa requires a training step and is therefore excluded; see Murali et al. 2018 for an IdTaxa-focused comparison").
- [x] **Q1.5b** — Update the methods table in `Taxonomic assignation.md` to mark IdTaxa as "evaluated elsewhere, not in this study" so reviewers see the decision was deliberate.
- [x] **Q1.6** — Pairwise method comparison **table**: per-rank ACC/F1 means (±sd) across DBs via `knitr::kable`. **Win matrix**: `slice_max` per (rank, db) → `count(tax_level, algo)` → pivoted table. Both added to `benchmark.qmd` § "Q1.6".
- [x] **Q1.7** — Statistical test: `lme4::lmer(values ~ algo + (1|db))` vs null per rank. LRT via `anova(m_null, m_full)`; χ², df, p-value table with significance stars. Added to `benchmark.qmd` § "Q1.7". Requires `lme4` on search path.

## Q2 — Database & simplification

The DB axis spans six labels (Unite, Unite_Fungi, EUK_ITS_v1_9_3, EUK_ITS_v1_9_3_Fungi, EUK_ITS_v1_9_3_Fungi_cut, EUK_SSU_v1_9_3_Fungi_cut) (decision 7). The qmd currently has DB on the y-axis of the NA plots but no systematic accuracy view.

- [x] **Q2.1** — **Per-DB accuracy figure**: all 4 methods faceted (`fig_q2_db.{pdf,png}`). Update `method_for_q2` to the best method from Q1.2 before finalising. Color = simplification type (full/Fungi/cut/Fungi+cut) via `db_meta` tibble.
- [x] **Q2.2** — **Simplification-effect table** + figure: `db_pair_delta()` helper computes Δ F1/MCC for three pairs; `knitr::kable` summary + `fig_q2_simplification.{pdf,png}`.
- [x] **Q2.3** — `EUK_ITS_v1_9_3_Fungi` and `EUK_ITS_v1_9_3_Fungi_cut` added to `values_map` in both `script_assign_taxo_parallel.R` and `values_map_for_qmd.R` (files confirmed present in both `dada2_format/` and `sintax_format/`). **Re-run `tar_make` to compute the new targets.** Note: `EUK_SSU_99_*` variants exist only in `sintax_format/` — needs a separate values_map block for non-dada2 methods if desired.
- [x] **Q2.4** — ~~Decide on mini DBs~~ **Decided out** (decision 4). `derive_mini()` stays in `make_databases.R` as a smoke-test tool only. Make sure no Results figure or table references `mini_*` rows. Filter them out of `values_map` selections when finalising figures (`dplyr::filter(!startsWith(db, "mini_"))`).
- [x] **Q2.5** — Discussion-side: tie DB-size to runtime via `benchmark_costs` (does a smaller DB pay for itself in accuracy/cost?). Figure code added to `benchmark.qmd` § "Q2.5 — DB size vs compute cost and accuracy" (`fig_q2_5_size_cost.{pdf,png}`).

## Q3 — Consensus voting

The qmd already builds three consensus columns. The manuscript table and `comparpq::resolve_taxo_conflict` document **five strategies × `strict` flag × `nb_agree_threshold`**.

- [x] **Q3.1** — Apply all five strategies × `strict ∈ {FALSE, TRUE}`: `unanimity`, `consensus`, `abs_majority`, `rel_majority`, `preference`. Loop in benchmark.qmd (extra_consensus_configs) adds 9 new columns covering all missing combinations.
- [x] **Q3.2** — **Consensus vs single-method figure**: F1/MCC per rank with consensus strategies plotted alongside the best single methods (conservatism order: preference→rel_majority→abs_majority→consensus→unanimity). `figures/fig_q3_consensus_vs_single.{pdf,png}`.
- [x] **Q3.3** — `nb_agree_threshold` sweep on `rel_majority` with values **1, 2, 3** (decision 6). `rel_majority_nb2_consensus` and `rel_majority_nb3_consensus` columns added; `figures/fig_q3_nb_threshold.{pdf,png}`.
- [x] **Q3.4** — `all_consensus_suffixes` vector drives the `ranks_df` loop; `tc_metrics_mock` now runs on all 12 consensus columns. `is_consensus`, `consensus_strict`, `nb_agree` columns added to `res_comp_tax` for downstream filtering.
- [x] **Q3.5** — **Method×DB matrix for `preference`**: `geom_tile` heatmap of F1 (Gen_sp) with method on y and DB on x; the preferred cell (sintax × EUK_ITS_v1_9_3) highlighted with a red border; preference-consensus F1 annotated. `figures/fig_q3_heatmap.{pdf,png}`.

## D1 — Multi-dataset coverage

Mock is done. The other three Bokulich-style datasets are partly missing.

| Dimension                   | Mock community                                                           | Cross-validation                                 | In silico                                                   | Biological community                                  |
| --------------------------- | ------------------------------------------------------------------------ | ------------------------------------------------ | ----------------------------------------------------------- | ----------------------------------------------------- |
| **Dataset ID**              | D1                                                                       | D1a                                              | D1b                                                         | D1c                                                   |
| **Data origin**             | Lab-made community of known taxa                                         | Reference database itself (held-out folds)       | Computationally simulated reads                             | Real environmental samples (Taudière 2018 endophytes) |
| **Ground truth**            | Known taxonomy from culture collection                                   | DB labels (held-out)                             | Input FASTA taxonomy                                        | None                                                  |
| **Sequencing errors**       | Real (wet-lab)                                                           | None (DB sequences)                              | Simulated (InSilicoSeq MiSeq model) / None (miaSim)         | Real (wet-lab)                                        |
| **Community structure**     | Fixed, known abundances                                                  | Uniform (one seq per taxon)                      | Controlled (zero-inflated lognormal / neutral-niche models) | Unknown, natural                                      |
| **"Novel taxon" possible?** | No — all taxa are in the DB                                              | No — sequences are drawn from the DB             | No (D1b-A) / depends on model (D1b-B)                       | Yes — environmental taxa may be absent from any DB    |
| **Metrics computable**      | Full: TP/FP/FN/TN/MCC/F1/TAR/TDR                                         | good/wrong/NA proportions per fold               | Full (same as mock)                                         | Agreement-only: Jaccard, Bray-Curtis, richness        |
| **Main analytical risk**    | Mock may not represent natural diversity                                 | DB leakage (remove_tested = TRUE/FALSE variants) | Sequencing model may not capture real error profile         | No objective benchmark possible                       |
| **Script / store**          | `script_dada2.R` + `script_assign_taxo_parallel.R` → `store_assign_taxo` | `script_cross_val.R` → `store_cross_val`         | `In_silico_simulation.qmd` (D1b-B); `store_*_insilico` for D1b-A | New `store_*_endophyte` (planned)                     |
| **Status**                  | ✅ complete                                                               | ✅ pipeline run (phases 2–3 done)                | 🔄 D1b-B (miaSim) in progress; D1b-A (InSilicoSeq) outlined | ❌ fastq fetch pending                                 |

### D1a — Cross-validation
- [x] Wire `cross_val()` into its own targets script (`script_cross_val.R` → `store_cross_val`). Inputs: each DB in `values_map$db`; each method. Outputs: a tibble per (method, db) with good/wrong/NA proportions per rank (averaged over folds via `cv_to_tidy()`). Note: dada2 single-bootstrap branch in `cross_val()` does not apply the bootstrap filter (known upstream limitation).
- [x] Add the **leaked** variant (set `remove_tested_sequences = FALSE` in `cross_val`). Both `remove_tested = TRUE` (standard) and `FALSE` (leaked) rows are in `cv_values_map`; the `leaked_suffix` column drives target naming.
- [x] Aggregation target: `cv_results` via `tarchetypes::tar_combine(bind_rows(!!!.x))`. Each per-(method,db,variant) tibble has columns `tax_level`, `mean`, `sd`, `metric`, `method`, `db`, `remove_tested`, `min_bootstrap`. Run with: `tar_make(script = "script_cross_val.R", store = "store_cross_val")`.
- [x] Smoke-test size control: `max_seq` parameter added to `cross_val()` (subsamples the DB before folding). Wired through `run_cv()` in `script_cross_val.R` via `cv_max_seq` in `config.R`. Three publication knobs in `config.R`: `cv_fold_number = 10L`, `cv_fold_tested = 2L` (raise to `cv_fold_number` for publication), `cv_max_seq = 100L` (set to `NULL` for full DB). Fix also applied: duplicate sequences within a fold are deduplicated before `create_fake_pq_from_refseq()` to avoid MiscMetabar's refseq validation error on EUK SSU databases.
- [ ] CV figures parallel to Q1 & Q2 (F1 per rank, per method × DB). Requires pipeline run first.

Exemple of a CV figure:

```{r}
cv_res <- tar_read("cv_results", store="store_cross_val")
cv_res |> dplyr::filter(metric=="good_classifications") |> arrange(desc(mean)) |> ggplot(aes(x=mean, y=factor(method), fill=db)) + geom_violin() + facet_grid(min_bootstrap~tax_level)
```



### D1b — In silico simulations (two sub-branches per decision 3)

**D1b-A — InSilicoSeq path (fastq from fasta).**
- [ ] Pick a curated input taxon list from UNITE Fungi (e.g. 50–200 species, one ref sequence each, spanning all 7 phyla).
- [ ] Run InSilicoSeq via the existing Docker recipe (already drafted in `Archives/some_bash_script` and in `Taxonomic assignation.md`): `iss generate --genomes mini_Unite.fasta --sequence_type amplicon --n_reads 10000 --abundance zero_inflated_lognormal --model MiSeq`. Reproducibly via 3 replicate seeds.
- [ ] Feed the simulated fastqs through `script_dada2.R` → `d_asv` → existing assignment loop. The "truth" table is the input taxon list.

**D1b-B — miaSim path (community matrices).**
- [ ] Generate community matrices via miaSim's neutral (Hubbell) and niche-based (Logistic, Lotka-Volterra) models. Each yields a phyloseq of *species × samples*, with known relative abundances.
- [ ] Bridge to assignment: the simulated species are real UNITE entries (so the refseq slot can be filled), but the abundances come from miaSim. Then run through the assignment loop the same way.
- [ ] Compare results between D1b-A (sequencing error included) and D1b-B (only community structure varies) — that pairing is itself a result.

- [ ] **D1b shared** — `In_silico_simulation.qmd` should grow two sections (one per sub-branch) and reuse the same `tc_metrics_mock` machinery. Truth tables differ but the metrics columns match.

### D1c — Biological communities (decision 2)

**Dataset:** Taudière et al. 2018, *Fungal Ecology*, [doi:10.1016/j.funeco.2018.07.008](https://doi.org/10.1016/j.funeco.2018.07.008). Tree endophyte ITS metabarcoding with `site` and `height` as sample modalities.

- [ ] Fetch the raw fastqs (or a derived OTU/ASV table) — check the paper's data availability statement for SRA / Dryad / ENA accession. If only OTU tables are public, this branch becomes "agreement-only" with no DADA2 rerun.
- [ ] If raw fastqs are available: drop them in `data/data_raw/rawseq_endophyte/` and run a copy of `script_dada2.R` with adapted primer / sample metadata. Stand up `store_dada2_endophyte` parallel to `store_dada2`.
- [ ] Confirm primers match the Pauvert ITS-1F/ITS2 pair — if not, update `fw_primer_sequences` / `rev_primer_sequences` in `config.R` *for this run only* (do not overwrite the canonical values).
- [ ] Run the assignment loop on `d_asv_endophyte`. No `tc_metrics_mock` (no ground truth).
- [ ] Agreement-only metrics: Jaccard / Bray-Curtis between method-pairs, per-sample richness comparison, consensus-vs-single calls. Decompose by `site` and `height`.
- [ ] **Methods-section text**: this is the "real conditions" leg of the four-data approach; emphasize that disagreement here is informative even without truth.

## M1 — Performance metrics

The set used in `tc_metrics_mock` already covers TP/FP/FN/TN/MCC/ACC/F1. The manuscript notes also mention Bokulich's TAR/TDR, Edgar's MC/OC/EPQ, and Bokulich-2018 over-/under-classification rates.

- [x] **M1.1** — Headline metric set fixed (decision 5): **F1 + MCC + TAR/TDR per rank** for mock-community results. TP/FP/FN/TN go to SM only. TAR (=PPV) and TDR (=TPR) added as `bind_rows` aliases after `res_comp_tax` creation in `benchmark.qmd`. `vote_algorithm` and `bootstrap_num` columns also added there for cleaner downstream filtering.
- [ ] **M1.2** — Add Edgar-style **misclassification rate / over-classification rate / EPQ** as alternative columns to `tc_metrics_mock` output, only for cross-validation datasets where "novel" vs "known" is well-defined.
- [x] **M1.3** — Confirm that the **fake taxa** (`add_shuffle_seq_pq` + `add_external_seq_pq`) are correctly counted as the TN denominator across methods — there is a comment in `ieauieau_tmp.R` (`fake_taxa_cond`) suggesting an in-progress reimplementation. Resolve which version is in `comparpq` now. **Resolved:** `ieauieau_tmp.R` does not exist on disk; the only implementation is in `comparpq/R/compare_taxo.R::tc_metrics_mock_vec()`, which correctly sets `fake_taxa_cond <- taxa_names(physeq) %in% fake_taxa_names` (matching `^fake_|^external_`) and counts TN as NA assignments among those taxa.

## M2 — Environmental metrics

`benchmark_costs` target already aggregates `wall_time_s`, `peak_resident_mb`, `mean_cpu_pct` per (method, db, bootstrap, vote). Missing: usage in the analysis.

- [x] **M2.1** — Plot **accuracy vs cost**: x = `wall_time_s` (from `tar_read(benchmark_costs, ...)`), y = F1 (at Gen_sp), color = method, labels = DB. `fig_m2_cost_species.{pdf,png}` (main) + `fig_m2_cost_allranks.{pdf,png}` (SM, all ranks in one row).
- [x] **M2.2** — CO₂eq via `greenAlgoR` (the pqverse package listed in the workspace CLAUDE.md). Feed `benchmark_costs` through it; one row per assignment target.
- [x] **M2.3** — Also log `benchmark_costs` for the **dada2 store** so the DADA2 preprocessing cost is attributable separately from the assignment cost (currently only `data_final/autometric_log_assign_taxo.txt` is post-processed; `autometric_log_dada2.txt` is dormant).

## Cross-cutting / housekeeping

- [ ] **C1** — Run `make_databases.R::derive_all_variants()` end-to-end on a clean machine to verify the dbpq-delegated derivations match the legacy outputs byte-for-byte (or close enough). Use a test directory to avoid clobbering current DBs.
- [x] **C2** — Update `analysis/benchmark.qmd` so it loads `dbpq` (currently sources `comparpq` files only) — needed for `dbpq::format2dada2` if any rerun is triggered from the qmd.
- [x] **C3** — `script_dada2.R` still has inline copies of the `config.R` constants. Migrate it to `source(here("config.R"))` next time you re-execute the DADA2 store.
- [ ] **C4** — Consider proposing `combine_taxo_assignments` and `cross_val` to `comparpq` (noted in `proposals_for_dbpq.md` as not-for-dbpq). Out of scope for this benchmark; nice for the broader pqverse.

## Decisions taken (2026-05-19)

1. **IdTaxa** — *not* included in the benchmark. The methods set is fixed at four: dada2, sintax, lca, blastn. The manuscript needs a one-line justification (e.g. "IdTaxa requires a separate training step and is therefore excluded from this cost-aware benchmark; see [reference] for an IdTaxa-focused comparison").
2. **Biological dataset** — endophyte dataset from **[Taudière et al. 2018, Fungal Ecology, doi:10.1016/j.funeco.2018.07.008](https://doi.org/10.1016/j.funeco.2018.07.008)**, using `site` and `height` as sample modalities.
3. **In silico tooling** — both **InSilicoSeq** (for fastq simulation from a curated fasta) **and miaSim** (for community-level structure from neutral/ecological processes). D1b has two sub-branches.
4. **Mini DBs in publication** — dropped. `mini_*` derivations stay in `make_databases.R` purely as a smoke-test tool for fast iteration on the analysis; they are excluded from the Results.
5. **Primary metric set** — **F1 + MCC** (Hleap-style) **+ TAR/TDR per rank** for mock-community results. F1 is the headline; MCC complements it (handles class imbalance); TAR/TDR are presence/absence checks for the mock.
6. **`nb_agree_threshold` values** — **1, 2, 3** (as default in Q3.3).
7. **SSU database reduction** — of the four EUK_SSU variants (full, Fungi-only, cut, Fungi+cut), only **EUK_SSU_v1_9_3_Fungi_cut** is retained. Rationale: the ITS region is a subregion of the 18S SSU gene, so ITS sequences are already present in any full SSU database — comparing all SSU variants would be redundant with the ITS database comparisons. The Fungi-filtered, primer-trimmed variant is the most directly comparable to the ITS databases and avoids inflating the DB axis with near-duplicate conditions.

## Suggested execution order (updated for the six decisions)

A path that gives you a draftable Results section fast, then expands.

1. **First pass (week-scale)** — Q1.1 → Q1.2 → Q1.3 → Q1.4 → M1.1 (TAR/TDR addition) → Q3.1 → Q3.2 → M2.1. End state: a complete mock-community Results section answering Q1 + Q3, with a cost panel and the agreed F1+MCC+TAR/TDR metric set. Manuscript can be partly drafted.
2. **Database axis (week-scale)** — Q2.1 → Q2.2 → Q2.3 → Q2.5. Adds Q2 to the draft. Remember to filter out `mini_*` rows (decision 4).
3. **Robustness via CV (week-scale)** — D1a fully. Re-runs Q1/Q2 figures on CV data; the agreement between mock and CV is itself a result.
4. **In silico** — D1b-A (InSilicoSeq) first because it reuses the existing DADA2 pipeline. Then D1b-B (miaSim) which only varies community structure.
5. **Biological** — D1c (Taudière 2018 endophytes). The fastq fetch is the first blocker; resolve it before standing up `store_dada2_endophyte`.
6. **Statistical layer** — Q1.6, Q1.7, Q2.5 (the cost-vs-DB-size analysis), M1.2 (Edgar-style MC/OC/EPQ for CV only).
7. **Polishing** — M2.2 (CO₂eq via greenAlgoR), M2.3 (DADA2 cost split), Q1.5b (IdTaxa-exclusion note in the methods table), C1–C4.

Stop at step 2 if the goal is a short methods note; go through step 5 for the full four-dataset manuscript per Bokulich 2020.
