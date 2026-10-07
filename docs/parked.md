# Parked items

Moved out of `ROADMAP.md` on 2026-09-30. Open items that are **not planned**:
kept so nothing is forgotten, revisited at the end of the project or when
their trigger happens. Nothing here blocks the manuscript. Ids are kept as they
were in the ROADMAP, since other files cite them. An item taken back into work
goes back to `ROADMAP.md` as one line.

Tags: `[Priority / facility]`, Priority = Critical → High → Medium → Low,
facility = easy (< 1 day) → moderate (a few days) → hard.

## 1. Cross-validation — deferred to the end of the project (OD 28)

**Decision of 2026-09-22 (developer): the cross-validation is kept as an idea
and revisited at the end, run only if it is then judged to deserve its compute
cost.** Nothing of its design is retracted and no code is touched; the items
below are parked, not dropped. The reasons are in `docs/objectives_design.md`
OD 28: the full-reference CV has never completed and its dada2 × `EUK_ITS_v2.1`
cells need ~58 GB each on a 62 GB machine, while the question it answers —
consistency of each method with each database's own labels at fixed reference
size — is secondary to mock accuracy. `run_cv_night.R` stages a resumable run
for whenever it is revisited. The CV run launched 2026-09-22 at 18:30 was
stopped by the developer after 52 s, with 0 targets built. Its framing item is
C22-P12 (`docs/manuscript_framing.md`); its history is `HISTORY.md` §2 D1a.

- [ ] `[Deferred / moderate]` CV figures parallel to Q1 and Q2 (F1 per rank, per method × database). → `HISTORY.md` §2 D1a
- [ ] `[Deferred / moderate]` Close the CV reference-design item: the Bokulich design is implemented (`cv_reduce_reference = FALSE`, 3 folds), the memory budget is measured; what remains is the production run and the Methods text (C22-P12). → `HISTORY.md` §2 D1a
- [ ] `[Deferred / moderate]` Negative controls in the CV folds, for a CV-side true-negative rate (TODO block of `R/cross_val.R`; same item as C22-P12's second half).
- [ ] **0.11a** `[Medium / easy]` Cross-validation databases (OD 27), chosen from the mock results; the CV grid keeps its 7 databases until then.
- [ ] **S6.9** `[Deferred with D1a / easy]` Open part only: the production `cross_val` run on the rebuilt `_cut` database (the code fix landed 2026-09-15). The code fix itself is done and is not waiting on anything. → `HISTORY.md` §1 phase 6
- [ ] **M1.2** `[Medium / moderate]` Edgar-style misclassification / over-classification rate / EPQ as alternative columns, for cross-validation only (where "novel" vs "known" is defined).
- Memory rule if it is run: dada2 × full `EUK_ITS_v2.1` (~58 GB) plus two fast workers gives 58.1 + 6.2 + 5.9 = 70.2 GB with `cv_n_workers_fast = 2`, which is why `run_cv_night.R` stages those cells to run alone (S8.6, `HISTORY.md` §6).

## 2. Upstream reports (other pqverse packages)

- [ ] **S6.6** `[Medium / easy]` Upstream reports: dbpq `format2sintax()` on UNITE general-release headers; MiscMetabar `assign_*()` parsing predictions by position (a missing rank shifts every column) and dropping empty taxa by default (`clean_pq = TRUE`, which silently skipped the controls); comparpq `tc_bar()` failing on an all-NA rank column; MiscMetabar `assign_blastn()` returning the input phyloseq unchanged with a message only when no hit passes `min_id` / `min_cover` / `bit_score` / `e_value`, so the caller gets a silent no-op instead of a taxonomy column — the failure that made every blastn result between 2026-05 and 2026-09-11 a measurement of the cover filter (S1.3, closed on 2026-09-22). Proposal: warn naming the filter and the hits lost at each step, or an `on_empty = c("warn", "error")` argument; the local guard in `combine_taxo_assignments()` (2026-09-10) only catches it one level down.
- [ ] **S7.6** `[Medium / easy]` Report the dbpq defects listed in `reviews/proposals_for_dbpq.md` §0 (UNITE separator, kingdom-anchored `filter_db()`, `timeout = Inf` warning, nested EUKARYOME archives, name qualifiers passed through).
- [ ] `[Low / easy]` comparpq robustness: `add_shuffle_seq_pq()` silently requires a `tax_table` and a `sample_data` with ≥ 2 variables; two guards with actionable messages.
- [ ] **C4** `[Low / easy]` Consider proposing `combine_taxo_assignments()` and `cross_val()` to comparpq (noted in `reviews/proposals_for_dbpq.md` as not-for-dbpq). Out of scope for this benchmark.

## 3. Side questions and optional datasets

- [ ] `[Low / easy]` Effect of `--maxaccepts 16` in sintax.
- [ ] `[Low / easy]` Effect of `--dbmask none` in sintax.
- [ ] `[Low / moderate]` Read filtering of Rolling et al. 2022 (`maxEE = 8`, `truncQ = 8`) on datasets other than Pauvert (Tedersoo first), before recommending it: on Pauvert it changes the unit set only, so the default settings are kept there (developer, 2026-10-06). Rerun `analysis/sandbox/test_read_filtering.R` adapted to the dataset. Tosadori & Bosch 2026 already found no significant difference (`maxEE = 2` / `truncQ = 2` against 8 / 8) on simulated data and on the Pauvert mock, which lowers the value of this item (`docs/manuscript_framing.md` §2). → `ROADMAP.md` §1, `data/snapshots/read_filtering/`
- [ ] `[Low / moderate]` SRR30413326 (decision 12) stays optional: fastqs not fetched, truth table on disk.
- [ ] `[Low / moderate]` Sciensano fungal mocks (SRA PRJNA814881, public 2024-05-11): 37 mock communities of fungi common in indoor air, MiSeq, each sequenced on ITS1 and on ITS2 (e.g. `DMC_4_ITS1`, `DMC_26_ITS2`), built to compare databases and software. Would suit Q4 (ITS1 vs ITS2 on the same communities), but no article publishes the composition of each mock (developer, 2026-10-06), so there is no truth to score against. Usable only if the composition tables are obtained from the authors (Sciensano).
- [ ] **M1.3** `[Low / moderate]` Distance to the true community (Jaccard / Bray-Curtis per sample against the mock). Parked 2026-09-30 (developer): only useful with molecular abundances, which are outside the project's scope; the per-unit F1 / PPV covers the presence part. The stale sandbox `analysis/sandbox/distance_to_true_community.qmd` was deleted; it is in the git history.
- [ ] **C1** `[Low / moderate]` Run `make_databases.R::derive_all_variants()` end-to-end on a clean machine, into a test directory, and check the derivations match.
