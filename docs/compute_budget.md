# Compute budget: what the pipeline computes, and how to spend less

Estimates of 2026-09-11 for this machine (6 cores, 62 GB RAM). No full run
with the current databases exists yet, so the minutes below are
extrapolated from:

- the sequential run of 2025-02 on complete databases (UNITE, EUKARYOME
  v1.9.3), read from `data/data_final/autometric_log_assign_taxo.txt`
  (wall time and resident memory per target);
- the cross-validation store of 2026-05 (`store_cross_val` meta: 10 folds
  tested, every target on the 10 000-record `mini_Unite`, bug S6.1);
- the record counts of the current databases.

Treat them as orders of magnitude (±50 %). Production decisions are in the
ROADMAP (Phase 8).

## 1. What is computed

```
make_databases.R   once per release ............................ ~15 min (done)
      │
dada2 project      3 samples ─► d_asv (200 ASVs) ................ ~10-15 min
      │
assign_taxo project
      d_asv ─┬────────────► d_asv_for_assignation        200 ASVs + 100 fake_ + 100 external_
             └─► ITSx ────► d_asv_itsx_for_assignation   same controls, ITS1 only
      │
      │   2 inputs × 8 databases × 12 method rows  =  192 assignment targets
      │
      │   12 method rows = dada2  × min_bootstrap {0.4, 0.5, 0.6}   ← 3 identical dada2 runs, threshold applied after
      │                  + sintax × min_bootstrap {0.4, 0.5, 0.6}   ← 3 identical vsearch runs, threshold applied after
      │                  + lca    × min_bootstrap {0.4, 0.5, 0.6}   ← lca ignores min_bootstrap: 3 identical results
      │                  + blastn × vote {rel_, abs_majority, unanimity} ← 3 identical makeblastdb + blastn, vote applied after
      │
      └─► d_all_taxo, benchmark_costs ........................... ~5 h wall (4-8 h), ~8 h of target time
      │
cross_val project
      8 databases × 4 methods × 2 variants (standard, leaked) = 64 targets
      each target = cv_fold_tested folds (10) × [assign N/10 queries against 0.9 N references]
      N = min(records, cv_max_seq) ............................. NULL: weeks · 10 000: ~7-8 h wall
      │
analysis/ chapters (render) ..................................... not measured
```

Every assignment and every CV fold runs single-threaded today: the pipelines
pass neither `multithread` (dada2) nor `nproc` (sintax, lca, blastn), although
`config.R` sets `n_threads = 4`. The 2025 logs show dada2 at 17 % CPU, i.e.
one core out of six.

## 2. Where the assignment time goes

One input (raw or ITSx), single-threaded, minutes summed over the 3 rows of
each method:

| Database | Records | dada2 | lca | blastn | sintax | Total |
|---|---:|---:|---:|---:|---:|---:|
| `EUK_ITS_v2.1` | 1 877 003 | 54 | 28 | 9 | 4 | **95** |
| `EUK_ITS_v2.1_Fungi` | 1 076 870 | 27 | 12 | 3 | 1.5 | **43** |
| `Unite_s_all_20250219` | 266 589 | 27 | 3 | 1.5 | 1 | **32** |
| `Unite_s_all_20250219_Fungi` | 168 030 | 13 | 2 | 1 | 0.5 | **17** |
| `Unite_all_20250219` | 156 820 | 13 | 2 | 1 | 0.5 | **17** |
| `Unite_all_20250219_Fungi` | 102 137 | 9 | 1.5 | 1 | 0.5 | **12** |
| `EUK_SSU_v2.1_Fungi_cut` | 116 746 | 3 | 6 | 1 | 0.3 | **10** |
| `EUK_ITS_v2.1_Fungi_cut` | 51 602 | 4.5 | 4.5 | 1 | 0.3 | **10** |
| **Total** | | **150 (64 %)** | **59 (25 %)** | **18 (8 %)** | **9 (4 %)** | **~236** |

- dada2 dominates, and its 48 targets share **one** worker (`dada2_ctrl`):
  that lane sets the wall time (~5 h for both inputs).
- The two large EUKARYOME ITS databases account for 58 % of the time.
- Memory: dada2 on EUK ITS v1.9.3 peaked at 44 GB; v2.1 is larger. With swap
  it will finish, but slowly if it really swaps.

Cross-validation (2026-05 store, N = 10 000, 10 folds, per database, both
variants): lca 134 min, dada2 ~20 min, blastn ~3 min, sintax ~1 min, i.e.
~158 min of target time per database, ~21 h for 8 databases, ~7-8 h wall
with 3 workers. The cost grows at least linearly with N (every record is a
query once) and up to quadratically (references grow too): with
`cv_max_seq = NULL` (3.8 M records in total) the run takes weeks at best.

## 3. Levers

**Status 2026-09-11: L1-L5 are implemented** (ROADMAP S8.1-S8.3):
`pipelines/assign_taxo.R` runs 64 computations instead of 192 assignments,
dada2 uses `assign_threads_dada2` threads. One consequence of B19 goes the
other way: sintax and lca now also query the 200 negative controls, which
roughly doubles their (small) query cost. Decided and implemented on 2026-09-14: L6 (`config.R::cv_threads`; 10 folds, `cv_max_seq = 5000L`), L7 (ITSx input on the 5 Fungi-filtered databases only: 156 assignment targets, 52 computations) and L8 (complete CV grid kept). On 2026-09-15 the SSU database was removed (132 assignment targets, 44 computations, 56 CV targets), `EUK_ITS_v2.1_Fungi_cut` became a trimmed copy of all its Fungi records instead of a 51 602-record subset (ROADMAP S6.9), and the CV queries are trimmed with cutadapt (one extra cutadapt call on `cv_oversample × cv_max_seq` records per CV target). The tables below describe the grid before these decisions.

"Same results" means the stored columns and target names do not change.

| Lever | What changes | Saves | Same results? | Effort |
|---|---|---|---|---|
| **L1** lca once per database × input | the 3 lca rows are copies of one computation | 2/3 of lca (~40 min per input) | yes (lca has no bootstrap) | easy |
| **L2** dada2 once, 3 thresholds after | one `assignTaxonomy(outputBootstraps = TRUE)`, then `tax[boot < 100 * b] <- NA` for b in 0.4/0.5/0.6 | 2/3 of dada2 (~100 min per input) | yes (this is what `minBoot` does) | moderate |
| **L3** sintax once, 3 thresholds after | one `assign_sintax(behavior = "return_matrix")`, thresholds on `taxo_bootstrap` | 2/3 of sintax (~6 min per input) | yes | moderate |
| **L4** blastn once, 3 votes after | one raw blast table per database × input (`blast_pq()`), 3 votes with `resolve_vector_ranks()`; one `makeblastdb` per database | 2/3 of blastn (~12 min per input) | yes | moderate |
| **L5** use the threads | `multithread` for dada2, `nproc` for sintax / lca (and blastn through `blast_pq()`), and `n_workers × n_threads ≤ 6` cores | dada2 lane ÷ 2-3; CV lca ÷ 2-3 | yes | easy |
| **L6** CV threads and folds | `nproc` in `run_cv()`, `cv_fold_tested = 5` | CV ÷ 2 (folds) × ÷ 2-3 (threads) | fewer folds = less precise means | easy |
| **L7** narrower ITSx axis | ITSx input only for a subset (e.g. the `_Fungi` databases, or one parameter row per method) | up to half of the assignment time | changes the design (decision 19) | easy |
| **L8** narrower CV grid | CV on a subset of databases, or without the `leaked` variant | proportional | changes the design | easy |

L1-L4 keep the 192 target names: one *compute* target per method × database
× input (32 per input instead of 96), then cheap *derive* targets that apply
the thresholds or votes and carry the current `full_name`s, so the store,
`d_all_taxo` and the chapters stay as they are. Their equivalence can be
checked on the mini store column by column.

## 4. Scenarios

| Scenario | assign_taxo (wall) | cross_val (wall) | Total |
|---|---|---|---|
| Today (`cv_max_seq = 10000`, 10 folds) | ~5 h | ~7-8 h | ~12-14 h |
| L1-L4 | ~1.7 h | ~7-8 h | ~9-10 h |
| L1-L5 | ~0.75-1 h | ~3-4 h (L5 on CV lca) | ~4-5 h |
| L1-L6 | ~0.75-1 h | ~1.5-2 h | ~2.5-3 h |
| L1-L6 + `cv_max_seq = 5000` | ~0.75-1 h | ~0.5-1 h | ~1.5-2 h |

Plus ~15 min of DADA2 and the chapter rendering. The first production run
of `dada2` also has to rebuild from `filtered` onwards (Phase 0, S0.3).
