# Benchmark plan — answering the manuscript questions

Working plan derived from `Coffre_principal/PROJETS/Projets IdEst/Dev packages R/MiscMetabar/autour_de_MiscMetabar/benchmark_taxo_assign/` (`Taxonomic assignation manuscript.md` + `Taxonomic assignation.md`). Use this as the running todo; tick items as you go.

## Main research questions (from the manuscript Results outline)

1. **Q1** — Effect of classification algorithm and its parameters
2. **Q2** — Effect of reference database and database simplification (Fungi-only, `_cut`, clustering, mini)
3. **Q3** — Effect of consensus voting across methods/databases/parameters

Plus the methodological scaffolding the manuscript announces:

- **D1** — Datasets: cross-validation, in silico, mock community, biological community (Bokulich 2020 four-data approach)
- **M1** — Performance metrics: TP/FP/FN/TN/MCC/ACC/F1/Precision/Recall (Hleap 2021 set)
- **M2** — Environmental metrics: CPU / wall time / CO₂eq

## Snapshot — what is already produced by the pipeline

| Piece | Where | Status |
|---|---|---|
| ASV phyloseq `d_asv` | `store_dada2`, from `script_dada2.R` | ✅ runs end-to-end |
| Fake taxa injection (shuffle + external) → `d_asv_for_assignation` | `script_assign_taxo_parallel.R` | ✅ TN material is in place |
| 4 methods × 7 DBs × bootstraps → `d_all_taxo` (combined `tax_table`) | `store_assign_taxo`, via `tar_combine` + `combine_taxo_assignments` | ✅ produced |
| Per-target runtime / memory (`benchmark_costs`) | `store_assign_taxo` | ✅ wired, not yet plotted |
| Mock truth table `taxo_mock` | `data/data_raw/metadata/taxo_mock.csv` | ✅ |
| NA cleanup + `Gen_sp` construction | `analysis/benchmark.qmd` § "Import value from store_assign_taxo" | ✅ |
| `tc_metrics_mock()` (TP/FP/FN/MCC/ACC/F1) per (method × db × rank) | `benchmark.qmd` via `comparpq` | ✅ produced as `res_comp_tax` |
| Three consensus strategies (unanimity, rel_majority, preference) applied | `benchmark.qmd` § "Create column using consensus..." | ✅ partial |
| Per-rank NA proportion plots | `benchmark.qmd` § "Proportion of NA" | ✅ exploratory plots |
| Cross-validation helper `cross_val()` | `R/cross_val.R` | ⚠️ source-able, no targets pipeline |
| In silico notebook | `In_silico_simulation.qmd` | ⚠️ ~70 lines, scratchpad only |
| Biological community dataset | Cregger et al. 2018, doi:10.1016/j.funeco.2018.07.008 (endophytes; site × height) | ⚠️ chosen, fastq fetch pending |
| `benchmark.qmd` ↔ parallel store wiring | `benchmark.qmd` still grep's `_all_taxo` | ❌ needs the three small edits (see CLAUDE.md "Parallel pipeline") |

## Q1 — Algorithm & parameters

Already producible from `res_comp_tax`. What's missing is the **systematic figure set** and a stats layer.

- [ ] **Q1.1** — Re-point `benchmark.qmd` at `tar_read(d_all_taxo, store = "store_assign_taxo")` (drop the `_all_taxo` grep, drop the `rename_ranks_pq(...)` block). Required before anything else in the qmd works on the parallel store.
- [ ] **Q1.2** — Lock the **per-rank × method** figure: ACC, F1, MCC across `c("dada2", "sintax", "lca", "blastn")` faceted by rank K→S. `plot_tc_metrics_mock()` already builds this; pin a final version and save as `figures/fig_q1_methods.{pdf,png}`.
- [ ] **Q1.3** — **min_bootstrap sweep figure** for dada2 / sintax: F1 or MCC as a function of `min_bootstrap ∈ {0.4, 0.5, 0.6}`, per rank. Aggregate from `res_comp_tax` (`bootstrap` column already extracted there).
- [ ] **Q1.4** — **Vote-algorithm sweep figure** for blastn: F1 or MCC across `c("rel_majority", "abs_majority", "unanimity")` × DB × rank. The data is already in `res_comp_tax`; needs a dedicated faceted panel.
- [x] **Q1.5** — ~~Add IdTaxa~~ **Decided out** (decision 1). Action: keep the `idtaxa` block commented out in `script_assign_taxo_parallel.R::values_map`, and add a one-paragraph justification in the Methods section ("IdTaxa requires a training step and is therefore excluded; see Murali et al. 2018 for an IdTaxa-focused comparison").
- [ ] **Q1.5b** — Update the methods table in `Taxonomic assignation.md` to mark IdTaxa as "evaluated elsewhere, not in this study" so reviewers see the decision was deliberate.
- [ ] **Q1.6** — Pairwise method comparison **table**: per-rank ACC/F1 means with sd across DBs, plus the **win matrix** (for each rank, which method ranks best on each DB). Useful in Discussion.
- [ ] **Q1.7** — Statistical test on whether method effects are significant once DB is controlled for: linear mixed model `F1 ~ method + (1|db)` per rank, or a permutation test. Decide which framing.

## Q2 — Database & simplification

The DB axis spans seven labels (Unite, Unite_Fungi, EUK_ITS_v1_9_3, EUK_SSU_v1_9_3, EUK_SSU_v1_9_3_Fungi, EUK_SSU_v1_9_3_cut, EUK_SSU_v1_9_3_Fungi_cut). The qmd currently has DB on the y-axis of the NA plots but no systematic accuracy view.

- [ ] **Q2.1** — **Per-DB accuracy figure** holding method fixed (best-performing method from Q1.2). One panel per rank, x-axis = DB, color = simplification type (full / Fungi / cut / Fungi_cut). Highlights whether Fungi-filtering or primer-trimming helps.
- [ ] **Q2.2** — **Simplification-effect table**: paired difference in F1/MCC between (a) `Unite` vs `Unite_Fungi`, (b) `EUK_SSU_v1_9_3` vs `EUK_SSU_v1_9_3_Fungi`, (c) `EUK_SSU_v1_9_3_Fungi` vs `EUK_SSU_v1_9_3_Fungi_cut`. Answers "does filtering / trimming help?" with paired effect sizes.
- [ ] **Q2.3** — **Add EUK_ITS_v1_9_3_Fungi and the _99_ variants** to `values_map` once `derive_clustered()` has finished producing them (it is in `make_databases.R::derive_all_variants`). The manuscript outline explicitly lists ITS Fungi and the 99% variants but they are commented out.
- [x] **Q2.4** — ~~Decide on mini DBs~~ **Decided out** (decision 4). `derive_mini()` stays in `make_databases.R` as a smoke-test tool only. Make sure no Results figure or table references `mini_*` rows. Filter them out of `values_map` selections when finalising figures (`dplyr::filter(!startsWith(db, "mini_"))`).
- [ ] **Q2.5** — Discussion-side: tie DB-size to runtime via `benchmark_costs` (does a smaller DB pay for itself in accuracy/cost?).

## Q3 — Consensus voting

The qmd already builds three consensus columns. The manuscript table and `comparpq::resolve_taxo_conflict` document **five strategies × `strict` flag × `nb_agree_threshold`**.

- [ ] **Q3.1** — Apply all five strategies × `strict ∈ {FALSE, TRUE}`: `unanimity`, `consensus`, `abs_majority`, `rel_majority`, `preference`. Currently only three (unanimity / rel_majority / preference) are wired in `benchmark.qmd`. Extend the block.
- [ ] **Q3.2** — **Consensus vs single-method figure**: F1/MCC per rank with consensus strategies plotted alongside the best single methods. Should show the conservatism gradient described in the notes (unanimity & consensus > abs_majority > rel_majority > preference).
- [ ] **Q3.3** — `nb_agree_threshold` sweep on `rel_majority` with values **1, 2, 3** (decision 6). At least one figure showing the trade-off (FN↑ vs FP↓ as threshold rises), faceted by rank.
- [ ] **Q3.4** — Run `tc_metrics_mock` on each consensus column and add them to `res_comp_tax` (currently they go through `rename_ranks_pq` then get measured implicitly; make the join explicit). Required before Q3.2 can be drawn.
- [ ] **Q3.5** — **Method×DB matrix for `preference`**: a small heatmap of accuracy with the preferred (method, db) on the diagonal. Justifies the preference-strategy claim.

## D1 — Multi-dataset coverage

Mock is done. The other three Bokulich-style datasets are partly missing.

### D1a — Cross-validation
- [ ] Wire `cross_val()` into its own targets script (`script_cross_val.R` → `store_cross_val`). Inputs: each DB in `values_map$db`; each method. Outputs: a tibble per (method, db) with TP/FP/FN/F1 per rank, plus a fold index.
- [ ] Add the **leaked** variant (set `remove_tested_sequences = FALSE` in `cross_val`). The manuscript explicitly tests this distinction.
- [ ] Aggregation target: `cv_results` joining method × db × strategy (10-fold vs leaked) × rank.
- [ ] CV figures parallel to Q1 & Q2 (F1 per rank, per method × DB).

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

**Dataset:** Cregger et al. 2018, *Fungal Ecology*, [doi:10.1016/j.funeco.2018.07.008](https://doi.org/10.1016/j.funeco.2018.07.008). Tree endophyte ITS metabarcoding with `site` and `height` as sample modalities.

- [ ] Fetch the raw fastqs (or a derived OTU/ASV table) — check the paper's data availability statement for SRA / Dryad / ENA accession. If only OTU tables are public, this branch becomes "agreement-only" with no DADA2 rerun.
- [ ] If raw fastqs are available: drop them in `data/data_raw/rawseq_endophyte/` and run a copy of `script_dada2.R` with adapted primer / sample metadata. Stand up `store_dada2_endophyte` parallel to `store_dada2`.
- [ ] Confirm primers match the Pauvert ITS-1F/ITS2 pair — if not, update `fw_primer_sequences` / `rev_primer_sequences` in `config.R` *for this run only* (do not overwrite the canonical values).
- [ ] Run the assignment loop on `d_asv_endophyte`. No `tc_metrics_mock` (no ground truth).
- [ ] Agreement-only metrics: Jaccard / Bray-Curtis between method-pairs, per-sample richness comparison, consensus-vs-single calls. Decompose by `site` and `height`.
- [ ] **Methods-section text**: this is the "real conditions" leg of the four-data approach; emphasize that disagreement here is informative even without truth.

## M1 — Performance metrics

The set used in `tc_metrics_mock` already covers TP/FP/FN/TN/MCC/ACC/F1. The manuscript notes also mention Bokulich's TAR/TDR, Edgar's MC/OC/EPQ, and Bokulich-2018 over-/under-classification rates.

- [x] **M1.1** — Headline metric set fixed (decision 5): **F1 + MCC + TAR/TDR per rank** for mock-community results. TP/FP/FN/TN go to SM only. Confirm `tc_metrics_mock` returns TAR and TDR explicitly — if not (current code returns TP/FP/FN/MCC/ACC/F1 etc.), add `TAR = TP/(TP+FP)` and `TDR = TP/(TP+FN)` either inside `comparpq::tc_metrics_mock()` or as a downstream `mutate()` in `benchmark.qmd`.
- [ ] **M1.2** — Add Edgar-style **misclassification rate / over-classification rate / EPQ** as alternative columns to `tc_metrics_mock` output, only for cross-validation datasets where "novel" vs "known" is well-defined.
- [ ] **M1.3** — Confirm that the **fake taxa** (`add_shuffle_seq_pq` + `add_external_seq_pq`) are correctly counted as the TN denominator across methods — there is a comment in `ieauieau_tmp.R` (`fake_taxa_cond`) suggesting an in-progress reimplementation. Resolve which version is in `comparpq` now.

## M2 — Environmental metrics

`benchmark_costs` target already aggregates `wall_time_s`, `peak_resident_mb`, `mean_cpu_pct` per (method, db, bootstrap, vote). Missing: usage in the analysis.

- [ ] **M2.1** — Plot **accuracy vs cost**: x = `wall_time_s`, y = F1 (at species rank), color = method, shape = DB. One panel per rank if useful. This is the manuscript's strongest selling point — it pairs the new methodological work with a practical recommendation.
- [ ] **M2.2** — CO₂eq via `greenAlgoR` (the pqverse package listed in the workspace CLAUDE.md). Feed `benchmark_costs` through it; one row per assignment target.
- [ ] **M2.3** — Also log `benchmark_costs` for the **dada2 store** so the DADA2 preprocessing cost is attributable separately from the assignment cost (currently only `data_final/autometric_log_assign_taxo.txt` is post-processed; `autometric_log_dada2.txt` is dormant).

## Cross-cutting / housekeeping

- [ ] **C1** — Run `make_databases.R::derive_all_variants()` end-to-end on a clean machine to verify the dbpq-delegated derivations match the legacy outputs byte-for-byte (or close enough). Use a test directory to avoid clobbering current DBs.
- [ ] **C2** — Update `analysis/benchmark.qmd` so it loads `dbpq` (currently sources `comparpq` files only) — needed for `dbpq::format2dada2` if any rerun is triggered from the qmd.
- [ ] **C3** — `script_dada2.R` still has inline copies of the `config.R` constants. Migrate it to `source(here("config.R"))` next time you re-execute the DADA2 store.
- [ ] **C4** — Consider proposing `combine_taxo_assignments` and `cross_val` to `comparpq` (noted in `proposals_for_dbpq.md` as not-for-dbpq). Out of scope for this benchmark; nice for the broader pqverse.

## Decisions taken (2026-05-19)

1. **IdTaxa** — *not* included in the benchmark. The methods set is fixed at four: dada2, sintax, lca, blastn. The manuscript needs a one-line justification (e.g. "IdTaxa requires a separate training step and is therefore excluded from this cost-aware benchmark; see [reference] for an IdTaxa-focused comparison").
2. **Biological dataset** — endophyte dataset from **[Cregger et al. 2018, Fungal Ecology, doi:10.1016/j.funeco.2018.07.008](https://doi.org/10.1016/j.funeco.2018.07.008)**, using `site` and `height` as sample modalities.
3. **In silico tooling** — both **InSilicoSeq** (for fastq simulation from a curated fasta) **and miaSim** (for community-level structure from neutral/ecological processes). D1b has two sub-branches.
4. **Mini DBs in publication** — dropped. `mini_*` derivations stay in `make_databases.R` purely as a smoke-test tool for fast iteration on the analysis; they are excluded from the Results.
5. **Primary metric set** — **F1 + MCC** (Hleap-style) **+ TAR/TDR per rank** for mock-community results. F1 is the headline; MCC complements it (handles class imbalance); TAR/TDR are presence/absence checks for the mock.
6. **`nb_agree_threshold` values** — **1, 2, 3** (as default in Q3.3).

## Suggested execution order (updated for the six decisions)

A path that gives you a draftable Results section fast, then expands.

1. **First pass (week-scale)** — Q1.1 → Q1.2 → Q1.3 → Q1.4 → M1.1 (TAR/TDR addition) → Q3.1 → Q3.2 → M2.1. End state: a complete mock-community Results section answering Q1 + Q3, with a cost panel and the agreed F1+MCC+TAR/TDR metric set. Manuscript can be partly drafted.
2. **Database axis (week-scale)** — Q2.1 → Q2.2 → Q2.3 → Q2.5. Adds Q2 to the draft. Remember to filter out `mini_*` rows (decision 4).
3. **Robustness via CV (week-scale)** — D1a fully. Re-runs Q1/Q2 figures on CV data; the agreement between mock and CV is itself a result.
4. **In silico** — D1b-A (InSilicoSeq) first because it reuses the existing DADA2 pipeline. Then D1b-B (miaSim) which only varies community structure.
5. **Biological** — D1c (Cregger 2018 endophytes). The fastq fetch is the first blocker; resolve it before standing up `store_dada2_endophyte`.
6. **Statistical layer** — Q1.6, Q1.7, Q2.5 (the cost-vs-DB-size analysis), M1.2 (Edgar-style MC/OC/EPQ for CV only).
7. **Polishing** — M2.2 (CO₂eq via greenAlgoR), M2.3 (DADA2 cost split), Q1.5b (IdTaxa-exclusion note in the methods table), C1–C4.

Stop at step 2 if the goal is a short methods note; go through step 5 for the full four-dataset manuscript per Bokulich 2020.
