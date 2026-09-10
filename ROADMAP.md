# ROADMAP — benchmar_assign_taxo

Working plan derived from `Coffre_principal/PROJETS/Projets IdEst/Dev packages R/MiscMetabar/autour_de_MiscMetabar/benchmark_taxo_assign/` (`Taxonomic assignation manuscript.md` + `Taxonomic assignation.md`). Use this as the running todo; tick items as you go.

**How this file is organised** (restructured 2026-09-10):

1. [Structure & maintenance](#1-structure--maintenance-critique-fable-2026-09-10) — engineering debt found by the structural review in [`critique_fable/`](critique_fable/README.md). Ordered in execution phases; each item links to the evidence.
2. [Open scientific work](#2-open-scientific-work) — every unticked manuscript item (D1a figures, D1b, D1c, M1.2, C1, C4, side questions).
3. [Decisions](#3-decisions-taken-2026-05-19) — the seven scope decisions of 2026-05-19 and the six structural decisions of 2026-09-10 (section 3b); do not reopen without a new dated entry.
4. [Reference](#4-reference) — research questions, dataset matrix, suggested execution order.
5. [Done log](#5-done-log) — every ticked item, kept verbatim for traceability.

Tags follow the workspace convention: `[Priority / facility]`, Priority = Critical → High → Medium → Low, facility = easy (< 1 day) → moderate (a few days) → hard.

---

## 1. Structure & maintenance (critique Fable, 2026-09-10)

Evidence and proposed code for every item: [`critique_fable/`](critique_fable/README.md). Ids (S0.1 …) are referenced from those files.

### Phase 0 — Quick wins (one sitting, no rerun needed)

Applied 2026-09-10. Notes: S0.1 done except the `store_*_mini` routing, which is part of S3.1 (`_targets.yaml`); `mini_db` now defaults to FALSE, so the **next `tar_make()` on `store_assign_taxo` / `store_cross_val` is a full production rerun** (hours). S0.3 changed the `file_refseq_taxo` and `filtered` commands, so `store_dada2` will also rebuild from `filtered` onwards on its next `tar_make()` (needs the `cutadaptenv` conda env). `tar_prune()` results: `store_dada2` 208 → 37 meta rows, `store_assign_taxo` 509 → 108, `store_cross_val` 138 → 82.

- [x] **S0.1** `[Critical / easy]` `mini_db` is `TRUE` in `config.R` and the store cannot tell mini from full runs. Read it from `BENCHMARK_MINI_DB` (default FALSE), print the flag + store path at pipeline start, and route smoke runs to `store_*_mini` (see S3.1). → [01 §B2](critique_fable/01_bugs_and_traps.md)
- [x] **S0.2** `[High / easy]` Create `figures/` and `dir.create()` it in the notebook setup; 26 `ggsave()` calls currently abort the render (ggplot2 4.0.3, `create.dir = FALSE`). → [01 §B3](critique_fable/01_bugs_and_traps.md)
- [x] **S0.3** `[High / easy]` `script_dada2.R`: fix `file_refseq_taxo` (`refseq/dada2_format/Unite_Fungi.fasta`), replace the dead `~/Nextcloud/IdEst/Projets/MiscMetabar/R/` source with the pqverse loader (S2.2), fix the doubled `paste(getwd(), here(...))` output paths. Required before any D1b-A / D1c DADA2 rerun. → [01 §B4](critique_fable/01_bugs_and_traps.md)
- [x] **S0.4** `[Medium / easy]` Fix the contradictory comments in `config.R` (`cv_fold_tested`, header about `script_dada2.R`) and `make.R` (step 3). → [01 §B5](critique_fable/01_bugs_and_traps.md), [05 drift table](critique_fable/05_analysis_and_docs.md)
- [x] **S0.5** `[Medium / easy]` `git rm --cached v2.31.0.tar.gz` (vsearch source tarball) and ignore it; name and move `41576_2023_679_MOESM2_ESM.ods` to `docs/references/`; move `design_analysis.svg` to `docs/`. → [03](critique_fable/03_layout_and_files.md)
- [x] **S0.6** `[Medium / easy]` Delete `.RData` (31 MB), `.Rhistory`, `renv/` (539 MB); set `RestoreWorkspace: No` / `SaveWorkspace: No` in the `.Rproj`; add `.claude/settings.local.json` to `.gitignore`. → [01 §B9](critique_fable/01_bugs_and_traps.md), [03](critique_fable/03_layout_and_files.md)
- [x] **S0.7** `[Medium / easy]` `tar_prune()` the three stores (509 meta rows in `store_assign_taxo`, three pipeline generations mixed) and add `tar_prune()` after each `tar_make()` in `make.R`. → [01 §B8](critique_fable/01_bugs_and_traps.md)

### Phase 1 — Correctness of the cost layer (M2)

Applied 2026-09-10 (`R/autometric_helpers.R`, logs under `data/data_final/autometric/<pipeline>/`, `d_asv_file` target). Not yet verified on a real run: the M2 figures stay open until the production rerun regenerates `benchmark_costs` (see M2 re-validation in section 2).

- [x] **S1.1** `[Critical / moderate]` `benchmark_costs` is built from the 2025-02-11 *sequential* run: crew workers never log a phase (no `log_start()` in the worker) and the autometric log is append-only across runs. Start a per-target logger inside the worker (`with_autometric()` helper, one file per target per run), aggregate from the directory, then **rerun and regenerate M2.1, M2.2, M2.3, Q2.5**. Those four items are ticked in the done log but must be considered open until this lands. → [01 §B1](critique_fable/01_bugs_and_traps.md)
- [x] **S1.2** `[High / easy]` Make `d_asv` a `format = "file"` dependency of `store_assign_taxo` (`store_dada2/objects/d_asv`) so a rebuilt DADA2 store invalidates the assignments. → [01 §B7](critique_fable/01_bugs_and_traps.md)
- [ ] **S1.3** `[High / easy]` All 18 blastn targets of the current (mini) store added no taxonomy column: `assign_blastn()` returns the input unchanged with only a message when no hit passes the score filters. Guard applied 2026-09-10: `combine_taxo_assignments()` now stops and names empty assignments. Open: confirm on the production rerun that blastn hits with the full databases, and propose a loud warning upstream in MiscMetabar. → [01 §B11](critique_fable/01_bugs_and_traps.md)

### Phase 2 — One source of truth for parameters and dependencies

Applied 2026-09-10. `tests/test_values_map.R` proves the new builders reproduce the legacy target names byte for byte, so no store object was orphaned by S2.1. Behaviour change to know about: crew workers now run the **dev checkout** of MiscMetabar / comparpq (`load_pqverse()` inside each target); before, they silently ran the installed CRAN MiscMetabar because `library("MiscMetabar")` put it in `tar_option_get("packages")`. `d_asv_for_assignation_fake` renamed `d_asv_shuffled` (invalidates the assignments, already forced by S0.1). The Q2.1 `db_meta` join was silently NA for every EUK database (v1_9_3 names); fixed by moving `db_meta` to `config.R`.

- [x] **S2.1** `[High / easy]` `R/values_map.R` with `build_methods_grid()`, `build_values_map()`, `build_cv_values_map()`, `full_name_for()`; `db_list` and `db_meta` move to `config.R`; delete `R/values_map_for_qmd.R`; add `tests/test_values_map.R` pinning the expected target names. → [04](critique_fable/04_pipeline_code.md)
- [x] **S2.2** `[High / easy]` `R/load_pqverse.R` + `PQVERSE_PKG_DIR` env var replacing the 15 absolute paths (plus one dead path in `script_dada2.R`) and the file-by-file `source()` of comparpq; drop `library("MiscMetabar")` where `load_all()` follows; capture `session_info()` as a target. → [02](critique_fable/02_dependencies_and_paths.md)
- [x] **S2.3** `[Medium / easy]` Rename `R/functions.R` → `R/create_fake_pq_from_refseq.R`; `Phyla` → `Phylum` in `script_dada2.R`; `d_asv_for_assignation_fake` → `d_asv_shuffled`; pass `targets_seed` to `derive_fake_ref()`; `tempfile()` instead of the fixed `test_refseq.fasta` in `cross_val()`. → [01 §B6](critique_fable/01_bugs_and_traps.md), [04](critique_fable/04_pipeline_code.md)

### Phase 3 — Layout

Applied 2026-09-10. Scripts moved to `pipelines/{dada2,assign_taxo,cross_val}.R`; `_targets.yaml` defines five projects (`TAR_PROJECT`), the `*_mini` ones own their stores and set `mini_db` (the `BENCHMARK_MINI_DB` variable of S0.1 is gone). `analysis/` is a Quarto project: `00_setup.R` + chapters `01_load_and_clean` → `06_cross_validation`, `sandbox/` for the unfinished material, `_cache/` for the inter-chapter objects; `benchmark.qmd` deleted. Stale `v1_9_3` names found and fixed on the way: `preference_pattern` (the preference consensus matched no column), Q2.2 pairs, Q3.5 preferred cell, the `tc_bar`/`tc_circle` exploratory calls. `data/README.md` documents every subfolder. `here::i_am()` is required in each chapter because `_quarto.yml` is itself a `here` root criterion.

- [x] **S3.1** `[Medium / moderate]` `_targets.yaml` with named projects (`dada2`, `assign_taxo`, `assign_taxo_mini`, `cross_val`, `cross_val_mini`); `make.R` becomes three `run_project()` calls; optionally move the three scripts to `pipelines/` and drop the `_parallel` suffix. → [03](critique_fable/03_layout_and_files.md)
- [x] **S3.2** `[Medium / easy]` Move `In_silico_simulation.qmd` to `analysis/`; write `data/README.md` describing every `data/` subfolder, including the undocumented `mock_hleap2021/`, `taxo_mock_SRR30413326.csv` and the empty `rawseq/SRR30413326/` (keep or delete?). → [01 §B10](critique_fable/01_bugs_and_traps.md), [03](critique_fable/03_layout_and_files.md)
- [x] **S3.3** `[Medium / moderate]` Turn `analysis/` into a Quarto project (`_quarto.yml`, `execute-dir: project`, `freeze: auto`) and split `benchmark.qmd` (1 381 lines) into `00_setup.R` + numbered chapters (load/clean, Q1, Q2, Q3, M2, CV) + `sandbox/` for lines 1077–1381; add a `save_fig()` helper. → [05](critique_fable/05_analysis_and_docs.md)

### Phase 4 — Tests and documentation

- [ ] **S4.1** `[Medium / easy]` `tests/run_all.R` (`testthat::test_dir`) and tests for `cv_to_tidy()`, `create_fake_pq_from_refseq()`, the awk/sed helpers of `make_databases.R`, and the cost aggregation (two-runs-same-phase fixture). → [04 §Tests](critique_fable/04_pipeline_code.md)
- [ ] **S4.2** `[High / easy]` Apply the documentation-drift table (README, CLAUDE.md, CONTEXT.md, `config.R`, `make.R`, `proposals_for_dbpq.md`): six DBs not nine, parallel script is the only one, `make.R` runs everything, CV knobs at publication values, `mini_db` documented, v2 names. Rule: CLAUDE.md describes what is on `make.R`'s path today; history goes to `Archives/README.md`. → [05 §Documentation drift](critique_fable/05_analysis_and_docs.md)
- [ ] **S4.3** `[Low / easy]` Align the project name: README title `benchmark_taxo_assign` vs folder/Rproj `benchmar_assign_taxo`. Keep the folder (Nextcloud + memory paths), fix the title. → [03](critique_fable/03_layout_and_files.md)

### Phase 5 — tidypq adoption (decision 2026-09-10, see decisions 8–13)

`tidypq` (`pqverse_pkg/tidypq`) now provides verbs at four scales (samples, taxa, occurrences, tree) plus `pq_to_tidy()` and `tax_table_to_df()`. The workspace rule is: prefer `tidypq::pq_to_tidy()` over `psmelt()` and never hand-roll the phyloseq-to-tibble pipeline. Call-site inventory and verb mapping: [06](critique_fable/06_tidypq_adoption.md).

- [ ] **S5.1** `[High / moderate]` Notebook: rewrite the `@tax_table` slot manipulations of `analysis/benchmark.qmd` (lines 58–103, 172–207, 247) with `tax_table_to_df()` / `mutate_taxa_pq()` / `rename_taxa_pq()`, the fake-taxa filters with `filter_taxa_pq()`, and the `psmelt()` at line 1090 with `pq_to_tidy()`. Do it during the S3.3 split so each chapter is rewritten once. Add `tidypq` to the loader (S2.2). → [06](critique_fable/06_tidypq_adoption.md)
- [ ] **S5.2** `[Medium / moderate]` Pipelines: `combine_taxo_assignments()` on `tax_table_to_df()` + `mutate_taxa_pq()`; `create_fake_pq_from_refseq()` and the `tax_tib` extractions in `cross_val()` on `tax_table_to_df()`; the `rename_samples()` / `otu_table()` block of `script_dada2.R` on `rename_samples_pq()`. Keep `tests/test_combine_taxo_assignments.R` green (it pins the exact column set). Invalidates all three stores: schedule with the production rerun that S0.1 already requires. → [06](critique_fable/06_tidypq_adoption.md)
- [ ] **S5.3** `[Low / easy]` Update CLAUDE.md "External code dependencies" (tidypq is now a dependency) and CONTEXT.md (`pq_to_tidy` entry).

### Answers to the open questions (2026-09-10)

- `figures/` is ignored entirely (`figures/*` + `!figures/.gitkeep`); figures live in the manuscript folder. Applied.
- `store_*/meta/meta` stay tracked; `tar_prune()` in `make.R` keeps the diff small.
- The reference MiscMetabar (and comparpq, dbpq, greenAlgoR, tidypq) is the dev checkout under `pqverse_pkg/`; S2.2 adds a `session_info` target to pin the commit.
- `docs/references/41576_2023_679_MOESM2_ESM.ods` is kept as is, source still to be filled in `docs/references/README.md`.
- `SRR30413326` is a planned second mock community (D1d below). `mock_hleap2021/` needs a separate decision: it holds Hleap 2021's *fish* mocks (12S/COI, `Chordata`), not fungal ITS.

---

## 2. Open scientific work

### Little side questions to explore

#### Software parameters

What is the effect of set --maxaccepts 16 in sintax ?
What is the effect of set --dbmask none in sintax ?

#### Databases building

What is the effect of taxing the Unite database with singletons as RefSeq only (Unite_RefS.fasta) ?
What is the effect of adding only some sequences to a FUNGI filtered database ?

### D1a — Cross-validation (remaining)

- [ ] CV figures parallel to Q1 & Q2 (F1 per rank, per method × DB). Pipeline has run (phases 2–3 done). Target home after S3.3: `analysis/06_cross_validation.qmd`.

Exemple of a CV figure:

```{r}
cv_res <- tar_read("cv_results", store="store_cross_val")
cv_res |> dplyr::filter(metric=="good_classifications") |> arrange(desc(mean)) |> ggplot(aes(x=mean, y=factor(method), fill=db)) + geom_violin() + facet_grid(min_bootstrap~tax_level)
```

- [ ] Known limitation carried from `R/cross_val.R`: the `dada2` branch calls `assignTaxonomy` with `minBoot = 0` and never applies `min_bootstrap` post-hoc; the `dada2_2steps` branch `stop()`s. Decide: fix, or document as "dada2 CV = no bootstrap filter" in Methods.
- [ ] From the TODO block in `R/cross_val.R`: add fake sequences to CV folds to obtain a true-negative rate (trade-off analysis); consider `nperm` permutations on top of k-fold (probably too slow). `[Low / moderate]`

### D1d — Second mock community, SRR30413326 (decision 12)

Truth table already on disk: `data/data_raw/metadata/taxo_mock_SRR30413326.csv` (16 fungal taxa: Russula ×3, Lactifluus ×3, Sydowia, Sporobolomyces, Phallus, …). `data/data_raw/rawseq/SRR30413326/` exists but is empty.

- [ ] Fetch the SRR30413326 fastqs into `data/data_raw/rawseq/SRR30413326/` (`fasterq-dump SRR30413326`), and record the study / primers / read layout in `data/README.md` (S3.2).
- [ ] Confirm the primer pair; if it differs from ITS-1F/ITS2, parameterise `script_dada2.R` (primers and sample-data file come from `config.R`; a `dataset` env var or a `_targets.yaml` project per dataset, see S3.1) rather than copying the script.
- [ ] Stand up `store_dada2_srr30413326` and `store_assign_taxo_srr30413326`, same `values_map` as the main run.
- [ ] Run `tc_metrics_mock()` against the second truth table; add the dataset as a facet/colour to the Q1 and Q2 figures (agreement between the two mocks is itself a result, as for CV).
- [ ] Decide what to do with `data/data_raw/mock_hleap2021/Mocks/` (five Hleap 2021 fish mocks, 32–387 sequences each): out of scope for a fungal ITS benchmark unless a cross-kingdom sanity check is wanted; otherwise delete.

### 🔄 Phase 4 — In silico simulations (manual actions)

| #   | Action                                                                                        | Command / location                                                                                                                                                     | Blocker for           |
| --- | --------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------- |
| 6   | Curate input taxon list from UNITE Fungi (50–200 species, one ref seq each, spanning 7 phyla) | Manual selection from `data/data_raw/refseq/`                                                                                                                          | D1b-A InSilicoSeq run |
| 7   | Run InSilicoSeq to generate simulated fastqs (3 replicate seeds)                              | `iss generate --genomes mini_Unite.fasta --sequence_type amplicon --n_reads 10000 --abundance zero_inflated_lognormal --model MiSeq` (see `Archives/some_bash_script`) | D1b-A pipeline        |
| 8   | Run DADA2 + assignment on InSilicoSeq fastqs                                                  | add `dada2_insilico` / `assign_taxo_insilico` projects to `_targets.yaml` (same scripts, own stores), then `TAR_PROJECT=dada2_insilico` and `TAR_PROJECT=assign_taxo_insilico` `tar_make()`          | D1b-A analysis        |

### D1b — In silico simulations (two sub-branches per decision 3)

**D1b-A — InSilicoSeq path (fastq from fasta).**
- [ ] Pick a curated input taxon list from UNITE Fungi (e.g. 50–200 species, one ref sequence each, spanning all 7 phyla).
- [ ] Run InSilicoSeq via the existing Docker recipe (already drafted in `Archives/some_bash_script` and in `Taxonomic assignation.md`): `iss generate --genomes mini_Unite.fasta --sequence_type amplicon --n_reads 10000 --abundance zero_inflated_lognormal --model MiSeq`. Reproducibly via 3 replicate seeds.
- [ ] Feed the simulated fastqs through `script_dada2.R` → `d_asv` → existing assignment loop. The "truth" table is the input taxon list.

**D1b-B — miaSim path (community matrices).**
- [ ] Generate community matrices via miaSim's neutral (Hubbell) and niche-based (Logistic, Lotka-Volterra) models. Each yields a phyloseq of *species × samples*, with known relative abundances.
- [ ] Bridge to assignment: the simulated species are real UNITE entries (so the refseq slot can be filled), but the abundances come from miaSim. Then run through the assignment loop the same way.
- [ ] Compare results between D1b-A (sequencing error included) and D1b-B (only community structure varies) — that pairing is itself a result.
- [ ] **[Medium / moderate]** Loop to create multiple communities at the last timestep of miaSim Hubbell model (multiple independent simulations for robustness). (source: `In_silico_simulation.qmd`:69)

- [ ] **D1b shared** — `In_silico_simulation.qmd` should grow two sections (one per sub-branch) and reuse the same `tc_metrics_mock` machinery. Truth tables differ but the metrics columns match.

### Phase 5 — Biological community (Taudière 2018) (manual actions)

| #   | Action                                                                                                                                   | Command / location                                                                                                                                              | Blocker for   |
| --- | ---------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------- |
| 9   | Check data availability in the paper (doi:10.1016/j.funeco.2018.07.008) — locate SRA / Dryad / ENA accession for raw fastqs or OTU table | Manual check of paper's data availability statement                                                                                                             | D1c setup     |
| 10  | Fetch raw fastqs (if available) into `data/data_raw/rawseq_endophyte/`                                                                   | `fastq-dump` / `wget` from accession                                                                                                                            | D1c DADA2 run |
| 11  | Confirm primers match ITS-1F/ITS2; if not, update `config.R` for this run only                                                           | Manual primer check against paper's Methods                                                                                                                     | D1c DADA2 run |
| 12  | Run DADA2 + assignment for biological community                                                                                          | add `dada2_endophyte` / `assign_taxo_endophyte` projects to `_targets.yaml` (same scripts, own stores), then `TAR_PROJECT=dada2_endophyte` and `TAR_PROJECT=assign_taxo_endophyte` `tar_make()` | D1c analysis  |

### D1c — Biological communities (decision 2)

**Dataset:** Taudière et al. 2018, *Fungal Ecology*, [doi:10.1016/j.funeco.2018.07.008](https://doi.org/10.1016/j.funeco.2018.07.008). Tree endophyte ITS metabarcoding with `site` and `height` as sample modalities.

- [ ] Fetch the raw fastqs (or a derived OTU/ASV table) — check the paper's data availability statement for SRA / Dryad / ENA accession. If only OTU tables are public, this branch becomes "agreement-only" with no DADA2 rerun.
- [ ] If raw fastqs are available: drop them in `data/data_raw/rawseq_endophyte/` and run a copy of `script_dada2.R` with adapted primer / sample metadata. Stand up `store_dada2_endophyte` parallel to `store_dada2`.
- [ ] Confirm primers match the Pauvert ITS-1F/ITS2 pair — if not, update `fw_primer_sequences` / `rev_primer_sequences` in `config.R` *for this run only* (do not overwrite the canonical values).
- [ ] Run the assignment loop on `d_asv_endophyte`. No `tc_metrics_mock` (no ground truth).
- [ ] Agreement-only metrics: Jaccard / Bray-Curtis between method-pairs, per-sample richness comparison, consensus-vs-single calls. Decompose by `site` and `height`.
- [ ] **Methods-section text**: this is the "real conditions" leg of the four-data approach; emphasize that disagreement here is informative even without truth.

### M1 — Performance metrics (remaining)

- [ ] **M1.2** — Add Edgar-style **misclassification rate / over-classification rate / EPQ** as alternative columns to `tc_metrics_mock` output, only for cross-validation datasets where "novel" vs "known" is well-defined.
- [ ] **[Medium / moderate]** Finalize distance-to-true-community metric (currently the "En chantier" section of `benchmark.qmd`, lines 1077–1200; moves to `analysis/sandbox/` with S3.3). (source: `analysis/benchmark.qmd`:1077)

### M2 — Environmental metrics (reopened)

- [ ] **M2.1 / M2.2 / M2.3 / Q2.5 re-validation** — blocked by S1.1. The figures exist but were computed from the 2025-02-11 sequential log. Re-tick in the done log once regenerated from a fresh parallel run.

### Cross-cutting / housekeeping (remaining)

- [ ] **C1** — Run `make_databases.R::derive_all_variants()` end-to-end on a clean machine to verify the dbpq-delegated derivations match the legacy outputs byte-for-byte (or close enough). Use a test directory to avoid clobbering current DBs. (Phase 3 verified the *current* machine; C1 is the clean-machine check.)
- [ ] **C4** — Consider proposing `combine_taxo_assignments` and `cross_val` to `comparpq` (noted in `proposals_for_dbpq.md` as not-for-dbpq). Out of scope for this benchmark; nice for the broader pqverse.
- [ ] **[Low / easy]** Optionally add "fake" taxa to benchmark for TRUE-negative / trade-off analysis. Not urgent. (source: `R/functions.R`:103 — note: that file is 50 lines; the TODO now lives in `R/cross_val.R`:22.)

---

## 3. Decisions taken (2026-05-19)

1. **IdTaxa** — *not* included in the benchmark. The methods set is fixed at four: dada2, sintax, lca, blastn. The manuscript needs a one-line justification (e.g. "IdTaxa requires a separate training step and is therefore excluded from this cost-aware benchmark; see [reference] for an IdTaxa-focused comparison").
2. **Biological dataset** — endophyte dataset from **[Taudière et al. 2018, Fungal Ecology, doi:10.1016/j.funeco.2018.07.008](https://doi.org/10.1016/j.funeco.2018.07.008)**, using `site` and `height` as sample modalities.
3. **In silico tooling** — both **InSilicoSeq** (for fastq simulation from a curated fasta) **and miaSim** (for community-level structure from neutral/ecological processes). D1b has two sub-branches.
4. **Mini DBs in publication** — dropped. `mini_*` derivations stay in `make_databases.R` purely as a smoke-test tool for fast iteration on the analysis; they are excluded from the Results.
5. **Primary metric set** — **F1 + MCC** (Hleap-style) **+ TAR/TDR per rank** for mock-community results. F1 is the headline; MCC complements it (handles class imbalance); TAR/TDR are presence/absence checks for the mock.
6. **`nb_agree_threshold` values** — **1, 2, 3** (as default in Q3.3).
7. **SSU database reduction** — of the four EUK_SSU variants (full, Fungi-only, cut, Fungi+cut), only **EUK_SSU_v2_Fungi_cut** (named `EUK_SSU_v1_9_3_Fungi_cut` at decision time) is retained. Rationale: the ITS region is a subregion of the 18S SSU gene, so ITS sequences are already present in any full SSU database — comparing all SSU variants would be redundant with the ITS database comparisons. The Fungi-filtered, primer-trimmed variant is the most directly comparable to the ITS databases and avoids inflating the DB axis with near-duplicate conditions.

## 3b. Decisions taken (2026-09-10)

8. **Figures are not versioned** — `figures/` is gitignored (only `.gitkeep`); the manuscript folder holds the reference copies.
9. **targets meta stays in git** — `store_*/meta/meta` remain tracked for provenance; `tar_prune()` runs after every `tar_make()`.
10. **Reference code = dev checkouts** — MiscMetabar, comparpq, dbpq, greenAlgoR and tidypq are loaded from `pqverse_pkg/` (`PQVERSE_PKG_DIR`, S2.2), not from CRAN; a `session_info` target pins the versions used for the manuscript.
11. **`41576_2023_679_MOESM2_ESM.ods` kept** in `docs/references/`, source to be documented.
12. **Second mock community SRR30413326** joins the benchmark as dataset D1d; the Hleap 2021 fish mocks await a separate decision.
13. **tidypq adopted in the notebook and the pipelines** (S5.1, S5.2): phyloseq slots are no longer manipulated by hand where a tidypq verb exists; `pq_to_tidy()` replaces `psmelt()`.

---

## 4. Reference

### Main research questions (from the manuscript Results outline)

1. **Q1** — Effect of classification algorithm and its parameters
2. **Q2** — Effect of reference database and database simplification (Fungi-only, `_cut`, clustering)
3. **Q3** — Effect of consensus voting across methods/databases/parameters

Plus the methodological scaffolding the manuscript announces:

- **D1** — Datasets: cross-validation, in silico, mock community, biological community (Bokulich 2020 four-data approach)
- **M1** — Performance metrics: TP/FP/FN/TN/MCC/ACC/F1/Precision/Recall (Hleap 2021 set)
- **M2** — Environmental metrics: CPU / wall time / CO₂eq

### Snapshot — what is already produced by the pipeline

| Piece                                                                    | Where                                                                              | Status                                                            |
| ------------------------------------------------------------------------ | ---------------------------------------------------------------------------------- | ----------------------------------------------------------------- |
| ASV phyloseq `d_asv`                                                     | `store_dada2`, from `script_dada2.R`                                               | ✅ runs end-to-end (store dated 2024-12-19; rerun blocked by S0.3) |
| Fake taxa injection (shuffle + external) → `d_asv_for_assignation`       | `script_assign_taxo_parallel.R`                                                    | ✅ TN material is in place                                         |
| 4 methods × 6 DBs × bootstraps → `d_all_taxo` (combined `tax_table`)     | `store_assign_taxo`, via `tar_combine` + `combine_taxo_assignments`                | ✅ produced (2026-05-21, with `mini_db = TRUE` — see S0.1)         |
| Per-target runtime / memory (`benchmark_costs`)                          | `store_assign_taxo`                                                                | ❌ stale (S1.1)                                                    |
| Mock truth table `taxo_mock`                                             | `data/data_raw/metadata/taxo_mock.csv`                                             | ✅                                                                 |
| NA cleanup + `Gen_sp` construction                                       | `analysis/benchmark.qmd` § "Import value from store_assign_taxo"                   | ✅                                                                 |
| `tc_metrics_mock()` (TP/FP/FN/MCC/ACC/F1) per (method × db × rank)       | `benchmark.qmd` via `comparpq`                                                     | ✅ produced as `res_comp_tax`                                      |
| Three consensus strategies (unanimity, rel_majority, preference) applied | `benchmark.qmd` § "Create column using consensus..."                               | ✅ partial                                                         |
| Per-rank NA proportion plots                                             | `benchmark.qmd` § "Proportion of NA"                                               | ✅ exploratory plots                                               |
| Cross-validation helper `cross_val()` + targets pipeline                 | `R/cross_val.R` + `script_cross_val.R` → `store_cross_val`                        | ✅ pipeline run (phases 2–3 done)                                  |
| In silico notebook                                                       | `In_silico_simulation.qmd`                                                         | 🔄 D1b-B (miaSim) in progress                                     |
| Biological community dataset                                             | Taudière et al. 2018, doi:10.1016/j.funeco.2018.07.008 (endophytes; site × height) | ⚠️ chosen, fastq fetch pending                                    |
| `benchmark.qmd` ↔ parallel store wiring                                  | reads `tar_read(d_all_taxo, ...)` directly                                         | ✅ (Q1.1)                                                          |

### D1 — Multi-dataset coverage

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

### Suggested execution order (updated 2026-09-10)

A path that gives you a draftable Results section fast, then expands.

0. **Structure Phase 0 + S1.1** (one or two sittings) — fixes the traps that would bite every later step (mini DB flag, figures dir, dada2 script paths, cost log). Then one production rerun of `assign_taxo` with `mini_db = FALSE`.
1. **First pass (week-scale)** — Q1.1 → Q1.2 → Q1.3 → Q1.4 → M1.1 (TAR/TDR addition) → Q3.1 → Q3.2 → M2.1. End state: a complete mock-community Results section answering Q1 + Q3, with a cost panel and the agreed F1+MCC+TAR/TDR metric set. Manuscript can be partly drafted. *(Done except the M2.1 re-validation.)*
2. **Database axis (week-scale)** — Q2.1 → Q2.2 → Q2.3 → Q2.5. Adds Q2 to the draft. Remember to filter out `mini_*` rows (decision 4). *(Done except the Q2.5 re-validation.)*
3. **Robustness via CV (week-scale)** — D1a fully. Re-runs Q1/Q2 figures on CV data; the agreement between mock and CV is itself a result. *(Pipeline run; figures open.)*
4. **In silico** — D1b-A (InSilicoSeq) first because it reuses the existing DADA2 pipeline. Then D1b-B (miaSim) which only varies community structure.
5. **Biological** — D1c (Taudière 2018 endophytes). The fastq fetch is the first blocker; resolve it before standing up `store_dada2_endophyte`.
6. **Statistical layer** — Q1.6, Q1.7, Q2.5 (the cost-vs-DB-size analysis), M1.2 (Edgar-style MC/OC/EPQ for CV only).
7. **Polishing** — M2.2 (CO₂eq via greenAlgoR), M2.3 (DADA2 cost split), Q1.5b (IdTaxa-exclusion note in the methods table), C1–C4, Structure Phases 2–4.

Stop at step 2 if the goal is a short methods note; go through step 5 for the full four-dataset manuscript per Bokulich 2020.

---

## 5. Done log

Kept verbatim from the pre-2026-09-10 ROADMAP. Items marked ⚠️ are ticked but depend on S1.1.

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

### Q1 — Algorithm & parameters

Already producible from `res_comp_tax`. What's missing is the **systematic figure set** and a stats layer.

- [x] **Q1.1** — Re-point `benchmark.qmd` at `tar_read(d_all_taxo, store = "store_assign_taxo")` (drop the `_all_taxo` grep, drop the `rename_ranks_pq(...)` block). Required before anything else in the qmd works on the parallel store.
- [x] **Q1.2** — Lock the **per-rank × method** figure: ACC, F1, MCC across `c("dada2", "sintax", "lca", "blastn")` faceted by rank K→S. `plot_tc_metrics_mock()` already builds this; pin a final version and save as `figures/fig_q1_methods.{pdf,png}`.
- [x] **Q1.3** — **min_bootstrap sweep figure** for dada2 / sintax: F1 or MCC as a function of `min_bootstrap ∈ {0.4, 0.5, 0.6}`, per rank. Aggregate from `res_comp_tax` (`bootstrap` column already extracted there).
- [x] **Q1.4** — **Vote-algorithm sweep figure** for blastn: F1 or MCC across `c("rel_majority", "abs_majority", "unanimity")` × DB × rank. The data is already in `res_comp_tax`; needs a dedicated faceted panel.
- [x] **Q1.5** — ~~Add IdTaxa~~ **Decided out** (decision 1). Action: keep the `idtaxa` block commented out in `script_assign_taxo_parallel.R::values_map`, and add a one-paragraph justification in the Methods section ("IdTaxa requires a training step and is therefore excluded; see Murali et al. 2018 for an IdTaxa-focused comparison").
- [x] **Q1.5b** — Update the methods table in `Taxonomic assignation.md` to mark IdTaxa as "evaluated elsewhere, not in this study" so reviewers see the decision was deliberate.
- [x] **Q1.6** — Pairwise method comparison **table**: per-rank ACC/F1 means (±sd) across DBs via `knitr::kable`. **Win matrix**: `slice_max` per (rank, db) → `count(tax_level, algo)` → pivoted table. Both added to `benchmark.qmd` § "Q1.6".
- [x] **Q1.7** — Statistical test: `lme4::lmer(values ~ algo + (1|db))` vs null per rank. LRT via `anova(m_null, m_full)`; χ², df, p-value table with significance stars. Added to `benchmark.qmd` § "Q1.7". Requires `lme4` on search path.

### Q2 — Database & simplification

The DB axis spans six labels (Unite, Unite_Fungi, EUK_ITS_v2, EUK_ITS_v2_Fungi, EUK_ITS_v2_Fungi_cut, EUK_SSU_v2_Fungi_cut) (decision 7; `v1_9_3` names at the time). The qmd currently has DB on the y-axis of the NA plots but no systematic accuracy view.

- [x] **Q2.1** — **Per-DB accuracy figure**: all 4 methods faceted (`fig_q2_db.{pdf,png}`). Update `method_for_q2` to the best method from Q1.2 before finalising. Color = simplification type (full/Fungi/cut/Fungi+cut) via `db_meta` tibble.
- [x] **Q2.2** — **Simplification-effect table** + figure: `db_pair_delta()` helper computes Δ F1/MCC for three pairs; `knitr::kable` summary + `fig_q2_simplification.{pdf,png}`.
- [x] **Q2.3** — `EUK_ITS_v1_9_3_Fungi` and `EUK_ITS_v1_9_3_Fungi_cut` added to `values_map` in both `script_assign_taxo_parallel.R` and `values_map_for_qmd.R` (files confirmed present in both `dada2_format/` and `sintax_format/`). **Re-run `tar_make` to compute the new targets.** Note: `EUK_SSU_99_*` variants exist only in `sintax_format/` — needs a separate values_map block for non-dada2 methods if desired.
- [x] **Q2.4** — ~~Decide on mini DBs~~ **Decided out** (decision 4). `derive_mini()` stays in `make_databases.R` as a smoke-test tool only. Make sure no Results figure or table references `mini_*` rows. Filter them out of `values_map` selections when finalising figures (`dplyr::filter(!startsWith(db, "mini_"))`).
- [x] ⚠️ **Q2.5** — Discussion-side: tie DB-size to runtime via `benchmark_costs` (does a smaller DB pay for itself in accuracy/cost?). Figure code added to `benchmark.qmd` § "Q2.5 — DB size vs compute cost and accuracy" (`fig_q2_5_size_cost.{pdf,png}`).

### Q3 — Consensus voting

The qmd already builds three consensus columns. The manuscript table and `comparpq::resolve_taxo_conflict` document **five strategies × `strict` flag × `nb_agree_threshold`**.

- [x] **Q3.1** — Apply all five strategies × `strict ∈ {FALSE, TRUE}`: `unanimity`, `consensus`, `abs_majority`, `rel_majority`, `preference`. Loop in benchmark.qmd (extra_consensus_configs) adds 9 new columns covering all missing combinations.
- [x] **Q3.2** — **Consensus vs single-method figure**: F1/MCC per rank with consensus strategies plotted alongside the best single methods (conservatism order: preference→rel_majority→abs_majority→consensus→unanimity). `figures/fig_q3_consensus_vs_single.{pdf,png}`.
- [x] **Q3.3** — `nb_agree_threshold` sweep on `rel_majority` with values **1, 2, 3** (decision 6). `rel_majority_nb2_consensus` and `rel_majority_nb3_consensus` columns added; `figures/fig_q3_nb_threshold.{pdf,png}`.
- [x] **Q3.4** — `all_consensus_suffixes` vector drives the `ranks_df` loop; `tc_metrics_mock` now runs on all 12 consensus columns. `is_consensus`, `consensus_strict`, `nb_agree` columns added to `res_comp_tax` for downstream filtering.
- [x] **Q3.5** — **Method×DB matrix for `preference`**: `geom_tile` heatmap of F1 (Gen_sp) with method on y and DB on x; the preferred cell (sintax × EUK_ITS_v1_9_3) highlighted with a red border; preference-consensus F1 annotated. `figures/fig_q3_heatmap.{pdf,png}`.

### D1a — Cross-validation

- [x] Wire `cross_val()` into its own targets script (`script_cross_val.R` → `store_cross_val`). Inputs: each DB in `values_map$db`; each method. Outputs: a tibble per (method, db) with good/wrong/NA proportions per rank (averaged over folds via `cv_to_tidy()`). Note: dada2 single-bootstrap branch in `cross_val()` does not apply the bootstrap filter (known upstream limitation).
- [x] Add the **leaked** variant (set `remove_tested_sequences = FALSE` in `cross_val`). Both `remove_tested = TRUE` (standard) and `FALSE` (leaked) rows are in `cv_values_map`; the `leaked_suffix` column drives target naming.
- [x] Aggregation target: `cv_results` via `tarchetypes::tar_combine(bind_rows(!!!.x))`. Each per-(method,db,variant) tibble has columns `tax_level`, `mean`, `sd`, `metric`, `method`, `db`, `remove_tested`, `min_bootstrap`. Run with: `tar_make(script = "script_cross_val.R", store = "store_cross_val")`.
- [x] Smoke-test size control: `max_seq` parameter added to `cross_val()` (subsamples the DB before folding). Wired through `run_cv()` in `script_cross_val.R` via `cv_max_seq` in `config.R`. Three publication knobs in `config.R`: `cv_fold_number = 10L`, `cv_fold_tested = 2L` (raise to `cv_fold_number` for publication), `cv_max_seq = 100L` (set to `NULL` for full DB). Fix also applied: duplicate sequences within a fold are deduplicated before `create_fake_pq_from_refseq()` to avoid MiscMetabar's refseq validation error on EUK SSU databases.

### M1 — Performance metrics

The set used in `tc_metrics_mock` already covers TP/FP/FN/TN/MCC/ACC/F1. The manuscript notes also mention Bokulich's TAR/TDR, Edgar's MC/OC/EPQ, and Bokulich-2018 over-/under-classification rates.

- [x] **M1.1** — Headline metric set fixed (decision 5): **F1 + MCC + TAR/TDR per rank** for mock-community results. TP/FP/FN/TN go to SM only. TAR (=PPV) and TDR (=TPR) added as `bind_rows` aliases after `res_comp_tax` creation in `benchmark.qmd`. `vote_algorithm` and `bootstrap_num` columns also added there for cleaner downstream filtering.
- [x] **M1.3** — Confirm that the **fake taxa** (`add_shuffle_seq_pq` + `add_external_seq_pq`) are correctly counted as the TN denominator across methods — there is a comment in `ieauieau_tmp.R` (`fake_taxa_cond`) suggesting an in-progress reimplementation. Resolve which version is in `comparpq` now. **Resolved:** `ieauieau_tmp.R` does not exist on disk; the only implementation is in `comparpq/R/compare_taxo.R::tc_metrics_mock_vec()`, which correctly sets `fake_taxa_cond <- taxa_names(physeq) %in% fake_taxa_names` (matching `^fake_|^external_`) and counts TN as NA assignments among those taxa.

### M2 — Environmental metrics

`benchmark_costs` target already aggregates `wall_time_s`, `peak_resident_mb`, `mean_cpu_pct` per (method, db, bootstrap, vote). Missing: usage in the analysis.

- [x] ⚠️ **M2.1** — Plot **accuracy vs cost**: x = `wall_time_s` (from `tar_read(benchmark_costs, ...)`), y = F1 (at Gen_sp), color = method, labels = DB. `fig_m2_cost_species.{pdf,png}` (main) + `fig_m2_cost_allranks.{pdf,png}` (SM, all ranks in one row).
- [x] ⚠️ **M2.2** — CO₂eq via `greenAlgoR` (the pqverse package listed in the workspace CLAUDE.md). Feed `benchmark_costs` through it; one row per assignment target.
- [x] ⚠️ **M2.3** — Also log `benchmark_costs` for the **dada2 store** so the DADA2 preprocessing cost is attributable separately from the assignment cost (currently only `data_final/autometric_log_assign_taxo.txt` is post-processed; `autometric_log_dada2.txt` is dormant).

### Cross-cutting / housekeeping

- [x] **C2** — Update `analysis/benchmark.qmd` so it loads `dbpq` (currently sources `comparpq` files only) — needed for `dbpq::format2dada2` if any rerun is triggered from the qmd.
- [x] **C3** — `script_dada2.R` still has inline copies of the `config.R` constants. Migrate it to `source(here("config.R"))` next time you re-execute the DADA2 store.
