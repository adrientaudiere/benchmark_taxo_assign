# ROADMAP — benchmar_assign_taxo

**The running todo, one line per item.** Split from the former ROADMAP on
2026-09-22 (item C22-P19): every ticked item, every dated narrative and every
past decision now lives in [`HISTORY.md`](HISTORY.md), which is the archive and
must be read for the measurements behind any item below. When an item here is
finished, tick it and move its narrative to `HISTORY.md`.

**Ids.**

| Prefix | Source |
|---|---|
| `0.x` | implementation plan of the 2026-09-17 design (`HISTORY.md` §0) |
| `S…` | structural review of 2026-09-10 ([`critique_fable/`](critique_fable/README.md)) |
| `C22-P…`, `C22-S…`, `C22-R…`, `C22-T…` | critique of 2026-09-22 ([`critique_2026-09-22.md`](critique_2026-09-22.md)), same numbering as that file: `P1–P13` conceptual, `P14–P20` and `S1–S7`, `R1–R6` code, `T9` a threat with no proposition, `P21` from its addendum, `N1` the second developer note of the addendum |
| `Q…`, `D1…`, `M1…`, `M2…`, `C…` | manuscript objectives |

Tags: `[Priority / facility]`, Priority = Critical → High → Medium → Low,
facility = easy (< 1 day) → moderate (a few days) → hard.

**Where the facts live.** The design is `docs/objectives_design.md` (decisions
1–27, "OD N"); the scoring is `docs/hleap_2021_metrics.md`; the positioning
against the literature is [`literature_2026-09-22.md`](literature_2026-09-22.md);
the glossary is [`CONTEXT.md`](CONTEXT.md); how to work here is
[`AGENTS.md`](AGENTS.md). A **new decision** is written in
`docs/objectives_design.md` §4 (as OD 28, 29, …) and summarised in
`HISTORY.md` §3d — not stated for the first time in this file.

---

## 0. Status (2026-09-22)

The 2026-09-17 grid is implemented: items 0.1–0.10 are done (`HISTORY.md` §0),
the per-unit truth (0.6), the Hleap scoring (0.7) and the chapters that read it
(0.8) landed on 2026-09-22. **The production `assign_taxo` run on the
nine-database grid landed on 2026-09-22** (12:26–15:58, 3 h 32, "587 completed,
53 skipped", 0 error, 0 warning): `store_assign_taxo` no longer holds the former
seven-database grid, and all 9 databases are present, 52 targets per `full`
database and 77 per `Fungi` / `Fungi+rep` (ITSx runs on the non-full databases
only). 19 `Unite_s_all_20250219_Fungi_rep` OTU targets are carried over from a
partial run the same morning (11:58), hash-valid. Still open in S6.5: `cross_val` on the new grid and
chapters 01–06 to re-render, so every figure must still be regenerated (S6.5,
S8.6). `store_dada2` was not rebuilt (unchanged since 2026-09-14). Numbers
quoted anywhere from `*_mini` stores are smoke-test numbers. The cost figures
of this run were corrected on 2026-09-22 (S1.6: the aggregation was at fault,
not the logging, so the existing logs were recoverable) and `benchmark_costs`
was rebuilt from them; whether that makes chapter 05 publishable without a
fresh run is the open half of S1.6.

Grid: 9 databases `[EUK_ITS_v2.1 | Unite_all_20250219 | Unite_s_all_20250219] ×
[all | fungi | fungi + rep]`, 3 inputs (ASV, OTU, ITSx ASV), 96 computations and
504 assignment rows per dataset, 21 parameter rows per database × input.

---

## 1. Before the first production run (blocking)

- [ ] **S6.5** `[Critical / moderate]` Production `dada2` → `assign_taxo` → `cross_val` on the nine-database grid, then re-render chapters 01–06 and re-tick the done-log items. → `HISTORY.md` §1 phase 6
- [ ] **S8.6** `[Critical / easy]` Production run budget: the dada2 × `EUK_ITS_v2.1` CV target needs ≥ 45 GB, `cv_n_workers_fast = 2` is the measured safe configuration; watch the first dada2 targets rather than trust the arithmetic. **The 45 GB predates the S1.6 correction of 2026-09-22 and should be re-derived before the CV launch**: on the corrected assignment-side logs `compute_dada2__EUK_ITS_v2.1` peaks at **58.1 GB**, and it is the CV cell of the same (method, database) pair that OOM-crashed six times. Re-read the budget off the corrected `benchmark_costs` peaks, not off the pre-correction figure. → `HISTORY.md` §1 phase 8
- [ ] **S1.6** `[High / moderate]` **Aggregation fixed and verified 2026-09-22, two halves left.** The samples were never mislabelled — `phase` is right in every row, only the *file name* lies — so the existing logs were recoverable and `benchmark_costs` has been rebuilt from them (504 rows, 0 NA). Correction measured on the 842 files: 20 of 330 phases had a wrong `peak_resident_mb` (9 **understated**, 11 overstated), worst `otu_compute_dada2__EUK_ITS_v2.1_Fungi_rep` at 1502 MB for a true 23 830 MB (15.9×). **The claim that wall time is unaffected was wrong**: 50 phases had a wrong `wall_time_s`, the same phase reading 3.0 s for a true 664 s. Open: **(a)** the leak guard in `with_autometric()`, deferred on purpose because it invalidates 604 targets (all 96 computations, a 3 h 32 re-run) while the aggregation fix invalidated only `benchmark_costs` — apply it in the same commit as the next production run, by adding `autometric::log_stop()` immediately before `autometric::log_start()` (a stop with nothing running is a no-op; without it a leaked logger makes the *next* target write no file at all). **(b)** decide whether the verification below retires the "clean run measured" condition, i.e. whether chapter 05 is publishable from the 2026-09-22 logs. → `HISTORY.md` §1 phase 1
- [ ] **S6.9** `[Critical / easy]` Open part only: the production `cross_val` run on the rebuilt `_cut` database (the code fix landed 2026-09-15). → `HISTORY.md` §1 phase 6

## 2. Manuscript framing (critique 2026-09-22)

The conceptual half of [`critique_2026-09-22.md`](critique_2026-09-22.md), read
together with its **addendum** (the developer's framing of 2026-09-22) and with
[`literature_2026-09-22.md`](literature_2026-09-22.md). Nothing here is applied:
each item is a proposal to triage.

- [ ] **C22-P1** `[Critical / easy]` Write the Introduction around 3–4 falsifiable predictions (P1a curation bias, P1b breadth not size, P1c abstention trade-off, P1d trimming no-op) instead of the "best trade-off" shopping question; the trade-off becomes the practical conclusion. The grid already tests all four at zero extra compute. → §1.1 C1
- [ ] **C22-P2** `[High / easy]` Pre-specify **one primary endpoint** and label everything else secondary or exploratory in the manuscript; the critique proposes F1 and PPV at Genus and Species on the two mocks. **Supersedes decision 5 and narrows decision 31** — left open on 2026-09-22 (developer), to settle when the first real figures exist; the accepted text goes to `docs/objectives_design.md` §4. → §1.1 C2
- [ ] **C22-P3** `[High / easy]` Held-out confirmation: any *recommendation* is selected on Pauvert and confirmed on Tedersoo, so a parameter set is never recommended from the mock that picked it. 0 compute; narrows decision 31 (which keeps the sweeps untuned on purpose). → §1.1 C3
- [ ] **C22-P4** `[Medium / easy]` Sensitivity run with the truthless units left out of the matrix (the `truth_depth` machinery already does this for ties); report whether the Q2 ranking survives. 42 % of the Tedersoo ASVs carry no truth. → §1.1 C4
- [ ] **C22-P5** `[High / easy]` Promote two findings buried in prose to Results: the **dada2 bootstrap collapse under taxonomic breadth** (`HISTORY.md` §2, "Databases building") and the **`_cut` no-op** (`docs/experimental_design.md` §5.5). Demote **Q5 only** — the critique's "demote Q3" is withdrawn by its addendum, Q3 is a primary axis. → §1.2 C5, addendum
- [ ] **C22-P6** `[Critical / moderate]` **Required, not optional** (addendum, because Q3 is a headline claim): analyse Q2 as paired within-source contrasts and write "UNITE vs EUKARYOME", not "databases" (n = 2 providers, the rest is a nested filtering ladder); report the Q3 gain against the **effective** number of independent voters, not the column count. → §1.3 T1
- [ ] **C22-P7** `[Medium / moderate]` One-species-out novelty test: remove every record of one mock species from a reference and score its ASVs. Turns a closed-book benchmark into an open one, answers Orsholm et al. 2026's central question, and doubles as the Q3 stress test. ≤ 1 grid run, or a scoring variant of the CV. → §1.3 T2, `literature_2026-09-22.md` §4
- [ ] **C22-P8** `[High / easy]` Headline F1/PPV/TPR **on real units only**; the controls report through their own rates (`NA_fake`, `ext_fungi`, `ext_correct`, `ctrl_assigned`). A reporting choice, not a rewrite — the per-unit cells are already separate. Otherwise the absolute metrics depend on `prop_fake = 0.5` and on the control count, which are design constants. → §1.3 T3
- [ ] **C22-P9** `[Medium / easy]` State the Pauvert/Tedersoo confound (ITS region, primer pair, sequencing run and community all differ at once): restrict the claim to "two mocks differing in ITS region and protocol"; only the in silico set or a third mock can separate region from community. → §1.3 T4
- [ ] **C22-P10** `[Medium / easy]` Sensitivity column with strict submitted-name-only scoring: synonym tolerance interacts with the UNITE/EUKARYOME nomenclature difference that Q2 measures. One-line variant of `tc_metrics_unit()`, the truth table already carries `Genus_accepted` / `Species_accepted`. → §1.3 T5
- [ ] **C22-P11** `[High / moderate]` Bootstrap CIs on the **paired** differences (10⁴ resamples of the ~150 scored units, Δ recomputed per resample) and McNemar on the discordant units. Cheap — the confusion cells are per unit — and it turns every Q2 claim into a supported one instead of eye-judged point estimates. → §1.3 T6
- [ ] **C22-P12** `[High / moderate]` D1a: give the cross-validation one explicit estimand in the Methods ("consistency of each method with each database's own labels at fixed reference size"), add the negative controls to the folds (own TODO in `R/cross_val.R`, finally a CV-side TNR), and report the per-database novelty rate beside every CV figure. → §1.3 T7
- [ ] **C22-P13** `[Medium / easy]` Keep M2 secondary: wall time and peak RAM as measured (after S1.6 and a clean rerun), CO₂eq only with a sensitivity note — if ±2× on the carbon intensity changes no ranking, say exactly that. → §1.3 T8
- [ ] **C22-T9** `[Low / easy]` One sentence in D1c saying what biological datasets can and cannot do: NA rates on real communities show that the mocks are not unrepresentative in their NA behaviour, they cannot support accuracy claims. → §1.3 T9
- [ ] **C22-N1** `[High / easy]` Declaration of interests: the harness and part of the benchmarked software share an author. Declare it, pin the package commits per run (the `session_info` targets already record them), and state that the engines (vsearch, BLAST+, dada2) are independent third-party tools — the wrapper can only flatter its own parameter plumbing, which `tests/test_assign_compute.R` audits against the direct calls. → addendum, `literature_2026-09-22.md` §4.4

## 3. Open scientific work

### After the first results

- [ ] **0.11a** `[Medium / easy]` Cross-validation databases (OD 27), chosen from the mock results; the CV grid keeps its 7 databases until then.
- [ ] **0.11b** `[Medium / moderate]` Final zoom (OD 16, 24, 25): `fungi + rep cut` on the two best databases (one UNITE, one EUKARYOME) × best methods, Pauvert only; compare `ext_fungi` against `fungi + rep` to see whether the k-mer methods are biased by the length difference.
- [ ] **0.11c** `[Low / easy]` Discussion notes: mono-nucleotide shuffle (M1.6), set-membership bias of the former metrics made moot by the per-unit truth, parameters reported on the whole mock rather than tuned on a third (hleap §5 point 8).
- [ ] **C22-P21** `[High / easy]` Q3 voter diversity, at zero compute on the existing columns: run `resolve_taxo_conflict()` on subsets (composition-only vs mixed families, one provider vs both, with and without the nested UNITE pair) to test the hypothesis the framing implies — **voting helps through error independence, so it gains from method-family and provider diversity, not from column count**. The Q3 figures cannot support a voting claim without it (C22-P6). → addendum, `literature_2026-09-22.md` §4.3

### Q4 — ITSx input

- [ ] `[Medium / easy]` Figures: F1 / MCC per rank with and without ITSx, per method × database; check whether the blastn gap closes and whether `min_cover` can return to 95 on ITSx-extracted ASVs. → `HISTORY.md` §2 Q4

### Q5 — ASVs vs OTUs

- [ ] `[Medium / easy]` Figures: F1 / MCC per rank for the three units (`ASV`, `OTU`, `OTU_derived`), per method × database. Frame as confirmation of Tosadori & Bosch 2026 with per-sequence truth (C22-P5). → `HISTORY.md` §2 Q5

### D1a — Cross-validation

- [ ] `[High / moderate]` CV figures parallel to Q1 and Q2 (F1 per rank, per method × database), after S6.5. → `HISTORY.md` §2 D1a
- [ ] `[Critical / moderate]` Close the CV reference-design item: the Bokulich design is implemented (`cv_reduce_reference = FALSE`, 3 folds), the memory budget is measured; what remains is the production run and the Methods text (C22-P12). → `HISTORY.md` §2 D1a
- [ ] `[Low / moderate]` Negative controls in the CV folds, for a CV-side true-negative rate (TODO block of `R/cross_val.R`; same item as C22-P12's second half).

### D1b — In silico (InSilicoSeq)

- [ ] `[Medium / moderate]` Stand up the `dada2_insilico` / `assign_taxo_insilico` projects in `_targets.yaml` (step 3 of `analysis/in_silico_simulation.qmd`); the fastqs of the three seeds exist since 2026-09-16.
- [ ] `[Medium / easy]` Score them against the per-record truth parsed from the UNITE headers with `comparpq::tc_metrics_unit()`, matched by the rule of `R/mock_truth.R::sanger_matches()`. Candidate use: separate ITS region from community (C22-P9).

### D1c — Biological communities

- [ ] `[Medium / moderate]` Decide which datasets enter the manuscript (GloSED, GSSP air, Korhonen wood, a leaf dataset, Tedersoo soil) and write the loader that turns their tables into phyloseq objects **with a refseq slot** — the methods need sequences. Read only through `NA_real`, `NA_fake` and `ext_fungi` (decision 26). → `HISTORY.md` §2 D1c
- [ ] `[Low / easy]` comparpq robustness: `add_shuffle_seq_pq()` silently requires a `tax_table` and a `sample_data` with ≥ 2 variables; two guards with actionable messages.

### D1d — Tedersoo Illumina mock, and a second mock

- [ ] `[High / easy]` Run `assign_taxo_tedersoo_illumina` (denoised since 2026-09-16: 238 ASVs; the truth table exists since 0.6). Blocked by S6.5 only in ordering.
- [ ] `[Medium / moderate]` Read retention on `ERR16773639` (15.3 % of pairs kept, against 53.7 % for `ERR16773638`): the 5′ reverse primer is found in 15 % of its R2 reads. Decide whether to keep both runs. → `HISTORY.md` §2 D1c
- [ ] `[Low / moderate]` SRR30413326 (decision 12) stays optional: fastqs not fetched, truth table on disk.

### M1 / M2 — metrics

- [ ] **M1.2** `[Medium / moderate]` Edgar-style misclassification / over-classification rate / EPQ as alternative columns, for cross-validation only (where "novel" vs "known" is defined).
- [ ] **M1.3** `[Medium / moderate]` Finalise the distance-to-true-community metric (`analysis/sandbox/distance_to_true_community.qmd`).
- [ ] **M1.6** `[Low / easy]` Discussion note on the mono-nucleotide shuffle (Hleap used `esl-shuffle`, di-nucleotide preserving), which likely inflates TNR for the k-mer methods. Same as 0.11c.
- [ ] **M2.1 / M2.2 / M2.3 / Q2.5** `[High / moderate]` Regenerate every cost figure after S1.6 and S6.5; keep them secondary (C22-P13).

### Side questions

- [ ] `[Low / easy]` Effect of `--maxaccepts 16` in sintax.
- [ ] `[Low / easy]` Effect of `--dbmask none` in sintax.

## 4. Code debt (critique 2026-09-22)

The code half of [`critique_2026-09-22.md`](critique_2026-09-22.md) §2. Its
verdict: ~7 500 lines of code for this grid is **not bloated**, and 2 000 lines
of tests is a virtue; what is not justified is where the size concentrates.

- [ ] **C22-S1 / C22-P14** `[High / moderate]` `cross_val()` is 625 lines: the metric block is written twice (vector and scalar `min_bootstrap` branches, ~120 lines each, diverging), `compute_by_tax_level` is duplicated in both and accumulates with `rbind` in a loop, and `res_assign_NA_classif` / `res_assign_NA_database` are computed then **dropped** although `docs/objectives_design.md` §1.6 promises them and `tests/test_cv_to_tidy.R` already feeds a `NA_classif` metric. Extract one `cv_score()` returning a long tibble, delete ~230 lines, return the two metrics or amend the design doc; also the dead `dada2_2steps` branch and the float `!=` invariant check. → §2.2 S1
- [ ] **C22-S2 / C22-P15** `[Medium / easy]` **Dead code done 2026-09-22**: `derive_clustered()` and `derive_no_parens()` deleted from `make_databases.R` (−40 lines; both defined, never called — verified by grep, not taken on the critique's word), withdrawn as §3 and §5 of `proposals_for_dbpq.md` under the rule *do not upstream dead code*, with the numbering kept stable because other documents cite those sections. `tests/run_all.R`: 12 files, 317 expectations, 0 failed. **Open**: replace the sed/awk header surgery (`general_header_rules`, `derive_no_pattern`, `derive_mini`) with the Biostrings + stringr idiom the same file already uses in `derive_fungi_rep()` — a French-locale awk once read `99.5` as `99`, and shell-string FASTA surgery has no type checking. → §2.2 S2
- [ ] **C22-S3 / C22-P16** `[High / moderate]` Three parsers of the same name grammar (`enrich_metrics()`, `R/score_assignments.R::assignment_db()`, the `na-prop` chunk of chapter 01) plus two definitions of "default settings" (`filter_default_settings()` / `default_names_for()`, the old item 0.4d). Keep the target names as identifiers, stop parsing them: `values_map` already carries `method`, `db`, `min_bootstrap`, `lca_cutoff`, `vote_algorithm`, `min_id` — join on `full_name`. ~90 lines and a class of silent corruption. → §2.2 S3
- [ ] **C22-S4 / C22-P17** `[Medium / moderate]` `analysis/01_load_and_clean.qmd` is 610 lines and four programs: move `clean_taxo()`, `add_gen_sp()` and `add_consensus_columns()` to `R/` next to `score_assignments.R` (where they become testable), keep the chapter to orchestration and invariant checks, and move the ~130 lines of exploratory figures to `sandbox/`. → §2.2 S4
- [ ] **C22-S5 / C22-P18** `[Low / easy]` `R/examples_cross_val.R` (134 lines, sourced by nothing) → `Archives/` or deleted; `Rplots.pdf`, `run_assign_taxo.log` and `v2.31.0.tar.gz` deleted or gitignored. → §2.2 S5
- [ ] **C22-S6 / C22-P19** `[Medium / moderate]` One source of truth for design facts. **Half done 2026-09-22**: the todo/log split is this file + `HISTORY.md`. Open: the same grid, computation counts and metric set are still stated in `README.md`, `AGENTS.md`, `CONTEXT.md`, `docs/objectives_design.md` and `docs/experimental_design.md` — make `docs/objectives_design.md` the only place a design fact is *defined*, and leave the others pointing at it. → §2.3 S6, S7
- [ ] **C22-P20** `[High / moderate]` Test balance: `cross_val()` (625 lines, all the CV metrics) has no test, while `tests/test_values_map.R` is 351 lines of byte-string pins that lock the name grammar without testing it. Add a synthetic 20-record fasta + a 4-query fold asserting `good + wrong + NA = 1` per rank and known counts, and property-test `full_name_for()` / `enrich_metrics()` as a round-trip. Note that several tests skip without vsearch / BLAST+ / conda, so a green `tests/run_all.R` proves less than it looks. → §2.4 R6
- [ ] **C22-R1** `[Low / easy]` `R/values_map.R:215-224`: callers must fake a vector (`rep("sintax", length(db))`) to satisfy an `ifelse()` inside `db_path_for()`; vectorise properly or use `dplyr::case_when()`.
- [ ] **C22-R2** `[Medium / easy]` `R/assign_compute.R:74-76`: `derive_assignment()` silently returns `physeq` unchanged when `computed` is NULL; warn at the derivation site, naming the (method, database, input).
- [ ] **C22-R3** `[Medium / easy]` `analysis/00_setup.R:119-122`: `db_shapes` slices a 9-shape vector with `[seq_along(db_list)]`, so a tenth database gets an `NA` shape and ggplot silently drops its points — the very failure the constant was added to prevent (S6.8). Recycle the vector instead of slicing it.
- [ ] **C22-R4** `[Low / moderate]` `pipelines/assign_taxo.R` and `pipelines/cross_val.R` duplicate the crew-controller block and the `ref_file_targets` machinery (~60 lines); one `R/targets_helpers.R`, or `tarchetypes::tar_map` / `tar_map_rep` instead of the hand-rolled `rlang::syms()` plumbing.
- [ ] **C22-R5** `[Low / easy]` `config.R` (318 lines) has started to behave like a module; the next dataset's columns go to a `data/data_raw/metadata/datasets.csv` read by `config.R`, not to another `tribble` column.

## 5. Structure & maintenance (2026-09-10, remaining)

Evidence in [`critique_fable/`](critique_fable/README.md); everything applied is
in `HISTORY.md` §1.

- [ ] **S6.6** `[Medium / easy]` Upstream reports: dbpq `format2sintax()` on UNITE general-release headers; MiscMetabar `assign_*()` parsing predictions by position (a missing rank shifts every column) and dropping empty taxa by default (`clean_pq = TRUE`, which silently skipped the controls); comparpq `tc_bar()` failing on an all-NA rank column; MiscMetabar `assign_blastn()` returning the input phyloseq unchanged with a message only when no hit passes `min_id` / `min_cover` / `bit_score` / `e_value`, so the caller gets a silent no-op instead of a taxonomy column — the failure that made every blastn result between 2026-05 and 2026-09-11 a measurement of the cover filter (S1.3, closed here on 2026-09-22). Proposal: warn naming the filter and the hits lost at each step, or an `on_empty = c("warn", "error")` argument; the local guard in `combine_taxo_assignments()` (2026-09-10) only catches it one level down.
- [ ] **S7.6** `[Medium / easy]` Report the dbpq defects listed in `proposals_for_dbpq.md` §0 (UNITE separator, kingdom-anchored `filter_db()`, `timeout = Inf` warning, nested EUKARYOME archives, name qualifiers passed through).

## 6. Housekeeping

- [ ] **C1** `[Low / moderate]` Run `make_databases.R::derive_all_variants()` end-to-end on a clean machine, into a test directory, and check the derivations match.
- [ ] **C4** `[Low / easy]` Consider proposing `combine_taxo_assignments()` and `cross_val()` to comparpq (noted in `proposals_for_dbpq.md` as not-for-dbpq). Out of scope for this benchmark.

## 7. Side analyses to focus on

Signals read off the production run of 2026-09-22 that no item above asks for,
all at **zero extra compute** — they are columns the grid already holds. Rates
below are the share of units given a **genus**, measured on the 195 mock ASVs,
the 97 shuffled `fake_` controls and the 100 non-fungal `external_` controls
separately (`analysis/sandbox/check_controls_by_simplification.R` for every
number below, `analysis/sandbox/check_blastn_min_cover.R` for the blastn ones).

**Read every control rate per simplification, never pooled.** Naming an
`external_` control is an error only on a `Fungi` database; on `full` it is
correct, and on `Fungi+rep` it is correct at Kingdom only (decisions 6, 15, 22).
Pooling the three also mixes the ITSx columns, which exist for the six non-full
databases only (`config.R::itsx_db_list`), with the raw ones, which exist for
all nine — a pooled mean therefore compares different database sets and
manufactures gaps that vanish when matched (found the hard way on 2026-09-22:
an apparent 3.5× ITSx effect on the `external_` controls was entirely this).

- [ ] **SA1** `[High / easy]` **The abstention trade-off, measured (P1c of C22-P1).** On the three `Fungi` databases at Genus, where naming a non-fungal control is a false positive by construction: blastn names 0.4 % of the `external_` controls and 0 % of the `fake_` ones for 0.753 of the mock ASVs, against lca 83.7 % / 8.2 % for 0.858, dada2 76.7 % / 8.8 % for 0.966 and sintax 35.2 % / 32.5 % for 0.966. blastn buys near-total abstention for about 21 points of genus recall, and lca — which has no identity threshold — names five of six non-fungal controls. This is the headline prediction of the Introduction proposed in C22-P1, already sitting in the data; it belongs in Results beside Q1, not in the control appendix. Pair it with C22-P8 (headline metrics on real units, controls through their own rates).
- [ ] **SA2** `[Medium / easy]` **ITSx makes lca markedly worse on the shuffled controls.** Matched on the six non-full databases at Genus (18 columns per cell), `lca` names 7.4 % of the `fake_` controls and `itsx_lca` 33.8 % — 4.6 times more, and its mock recall drops too (0.858 → 0.803, the only method ITSx costs anything) — while `itsx_dada2` improves on `dada2` (14.0 % → 7.2 %), sintax barely moves (31.2 % → 28.7 %) and blastn is unchanged at 0 %. Check first whether the denominator is intact (what an `itsx_*` column holds for a unit from which ITSx extracted nothing: dropped, kept whole, or NA) before reading any mechanism into it; a shuffled sequence with no recognisable ITS is exactly where that convention bites. Feeds Q4, whose current line only asks about F1 / MCC and the blastn gap.
- [ ] **SA3** `[Medium / easy]` **ITSx is a no-op for blastn on every unit class.** Matched by simplification, `blastn` and `itsx_blastn` give the same rates to three decimals (mock 0.753, `external_` 0.004 on `Fungi` and 0.175 on `Fungi+rep`, `fake_` 0). Q4 asks whether ITS extraction closes the blastn gap and whether `min_cover` can return to 95 on extracted ASVs: on the first half the production answer looks like "it changes nothing", which is a publishable negative result and a cheap one — confirm it on F1 rather than on assignment rate before writing it down.
