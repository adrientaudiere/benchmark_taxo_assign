# Code review of `benchmar_assign_taxo`, 2026-10-05

Written overnight on 2026-10-05 at the developer's request, as the
"Code review, rerun from scratch" item of `ROADMAP.md` §4. It also covers the
two items that line names: the stray files (C22-S5) and the design facts
repeated across the docs (C22-S6). The tone is deliberately harsh. Every
finding comes with the evidence and a fix precise enough to apply without
re-investigating.

**Nothing in the project was changed.** This file is the only output. Checked
at the end: `git status` shows this file as the only new entry, and no file
under `store_*/meta/`, `data/data_final/`, `analysis/_cache/` or `figures/`
was modified. The checks were run read-only against the current stores,
through scratch scripts kept outside the project. They are reproduced in the
appendix. Questions for the developer are grouped in §11; no
finding below depends on an answer unless it says so.

**What was read.** All of `config.R`, `make.R`, `make_databases.R`,
`.Rprofile`, `_targets.yaml`, `.gitignore`, `R/*.R` (19 files),
`pipelines/*.R` (5), `analysis/00_setup.R`, chapters 01–05, 07, the code of
08, 06 and the in silico chapter, `tests/` (15 files), the night scripts,
`patches/`, `AGENTS.md`, `ROADMAP.md`, and the 2026-09-22 critique. The
sandbox and the other documents were checked where a finding pointed to them.
Where the project calls them, the MiscMetabar, comparpq, greenAlgoR and dada2
sources were read too.

**What was run** (read-only): the test suite; refseq comparisons on the 9
ITSx inputs of the 5 stores; the bootstrap matrices of every dada2
computation; `benchmark_costs` per method; a replay of the external-control
draw; record lengths of the reference files; label checks on the Sanger
tables and reference headers; `tar_manifest()` against `tar_meta()` on 3
projects; git status of the project and of the pqverse checkouts; two fresh
dada2 computations (another seed, then the target's own seed) compared with
their stored twin; a tested process-tree memory sampler; and the test
snippets proposed in §7.

## Summary

The code is better than most research pipelines. It has file-tracked
references, isolated smoke-test stores, a fast green test suite (380
expectations, 38 s) and a per-unit truth, and it found and patched a real
dada2 bug. **The weak points are not where the effort went.** The
measurements and claims around the code have drifted away from what the code
does, and three of them reach the Results:

- the Q4 region comparison and its description (§1.1);
- the cost and CO₂ figures (§1.3, §1.4);
- the description of the controls in the Methods (§1.2).

A fourth, the seed variance of dada2 and sintax (§1.9), bounds how small an
effect the Results can claim.

| # | Finding | Severity | Cost |
|---|---|---|---|
| 1.1 | **Q4 compares regions on units ITSx never trimmed.** The ITS1 / ITS2 / full inputs keep the undetected units on their full amplicon (Tedersoo: 22 / 35 / 47 units). Chapters 01 and 07 and the ROADMAP say the opposite. Re-scored on the 150 units trimmed in every region (out of the 183 of Q4.4), the Tedersoo conclusion holds (Genus F1 full 0.935 > ITS1 0.932 > ITS2 0.923); only dada2's per-method order flips, within 0.003. The in silico "full" input is the ASV itself (0 of 71 cut by ITSx), which is already the full ITS: its UNITE templates have no SSU or LSU flank for ITSx to anchor on. It needs a sentence, not a fix. | Medium (Tedersoo), Low (in silico "full") | re-render 01, 07 |
| 1.3 / 1.4 | **The cost layer does not see vsearch or BLAST+** (child processes): their "peak RAM" is a flat ~0.9 GB, the idle R parent, at ~0 % CPU. **The CO₂ charges each job the machine's 62 GB.** Fix both before the M2.4 timing run (a tested sampler is given). | High (M2) | code now, measured by M2.4 |
| 1.6 | **`EUK_ITS_v2.1_Fungi_cut`, the fungal part of the planned zoom database, predates OD 32**: 11 records > 5 000 bp, including the 32 568 bp clone OD 32 removed. The zoom database itself (`fungi + rep cut`) has no derivation code yet. Root cause: in `make_databases.R` a file that exists counts as up to date, and writes are not atomic. | High (0.11b) | rebuild 2 files + write the zoom derivation |
| 1.7 | **The current scores come from uncommitted code**: OD 36 lives in a dirty comparpq checkout, and the project has no commit since 2026-09-22. `pqverse_versions()` cannot tell a dirty checkout from a clean one, and no chapter records versions. | High | commit + 20 lines |
| 1.9 | **Each dada2 column is one draw of its bootstrap.** Another seed moved dada2's Genus PPV on real Pauvert units from 0.937 to 0.924, the same order as several effects the chapters interpret. PPV is a co-primary endpoint. Even at a fixed seed, dada2 is not bit-reproducible (tie-break on `std::random_device`), but that part is small: 1–2 Species calls. | Medium | measure the seed spread now (~30 min) |
| 1.2 | **The external controls miss the 8 singleton phyla** (`sample(which(...), 1)` on a length-one vector): 65 of 73 phyla. The test passes by luck of its seed. | Medium | code now; regenerate = full rerun (Q2) |
| 3.1 | **The smoke-test project writes its cost logs into the production folder**: one case found, none affecting a current phase yet. | Medium | 1 line |
| 6 | **AGENTS.md, README and CONTEXT describe the pre-OD 34 grid** (96 computations, 504 rows, `d_vs` input, lca {0.8, 0.9, 1}, 10 folds, `n_workers = 3`). | Medium | docs; `ETAT.md` |

Everything else is in §2–§9. §10 orders all fixes by what they rerun.
**Every fix that invalidates assignments should go in before the M2.4
timing run, which reruns every computation anyway.**

## How to read this report

Each finding gives:

- **Severity**: Critical (a published number is wrong), High (a conclusion or a
  figure is biased, or the next run would be wasted), Medium (wrong but
  bounded, or a trap for the next change), Low (latent, style, hygiene).
- **Where**: `file:line`.
- **Evidence**: what was read or measured, and the measured output.
- **Fix**: code, or the exact change.
- **Cost**: what the fix invalidates. In a `{targets}` project this decides
  whether a fix is cheap. Three classes are used throughout:
  - **[0]**: no rerun (docs, dead code, a helper nothing has called yet);
  - **[R]**: re-render one or more chapters, no reassignment;
  - **[A]**: invalidates assignment targets (hours of night runs).

---

## 1. Confirmed defects (measured on the current stores)

### 1.1 Q4: the "region" inputs contain units ITSx never trimmed. **Medium** for Tedersoo, **Low** for the in silico "full" region [R]

**Where.** `pipelines/assign_taxo.R:225-235` (`itsx_pq(keep_undetected =
TRUE, duplicated_seqs = "merge")`), `R/itsx_inputs.R:34`. The wrong reading
appears in `analysis/07_q4_its_region.qmd:186-191`,
`analysis/01_load_and_clean.qmd:547-553` (and in `score_common_units()`,
`:556-579`), and in `ROADMAP.md:49` and `:53-54`.

**What the code does.** `itsx_pq(keep_undetected = TRUE)` keeps every unit
whose region ITSx did not find, **with its full, untrimmed sequence**
(MiscMetabar `R/itsx.R:398-426`: `keep` stays `TRUE`, and only `detected`
units get a new sequence). Such a unit is then assigned in the "ITS1" or
"ITS2" input on the whole amplicon. Nothing in the tax_table records it:
`ITSx_origin` is the organism call, not a detection flag.

**What the chapters say.** Chapter 07 (line 187): "an ASV without the region
is not assigned on it. **So the regions are scored on different subsets of
units**". Chapter 01 (line 548) says the same. ROADMAP line 49 reads "ITSx cuts
ITS1 from 213 and ITS2 from 199 of 237 real ASVs". These sentences are wrong:
213 and 199 are the unit counts **after merging** the units that became
identical once trimmed, and they include untrimmed units.

**Evidence** (the refseq of `d_<region>_asv_real` compared with the refseq of
`d_asv` under the same taxon name; "unchanged" = ITSx did not change the
sequence):

| Dataset | Region | `d_asv` | units | merged away | unchanged (assigned on the full amplicon) |
|---|---|---:|---:|---:|---:|
| Pauvert | ITS1 | 196 | 195 | 1 | 0 |
| Tedersoo | ITS1 | 238 | 213 | 25 | **22** |
| Tedersoo | ITS2 | 238 | 199 | 39 | **35** |
| Tedersoo | full | 238 | 237 | 1 | **47** |
| in silico | ITS1 | 71 | 71 | 0 | 1 |
| in silico | ITS2 | 71 | 68 | 3 | 1 |
| in silico | full | 71 | 71 | 0 | **71 (all)** |
| wood | ITS2 | 918 | 875 | 43 | 41 |
| air | ITS2 | 2 590 | 2 572 | 18 | 13 |

On Tedersoo (154 bp of 18S and 39 bp of 28S around the ITS), "unchanged" does
mean "region not detected". On wood and air, published ITS2 OTUs may simply
have no flank left to trim, so "unchanged" is not evidence of a failure there.

**Q4.4 (the "paired" reading) is impure too, but its conclusion survives.**
Of the 183 "ASVs common to the three regions" (`res_paired`), only **150** are
trimmed in all three regions. **20** are untrimmed in all three, so they enter
the ITS1, ITS2 and full scores with the same sequence and the same answers:
they dilute differences without creating any. 13 more are trimmed in some
regions and not in others. Re-scored with chapter 01's own functions
(`clean_taxo()`, `add_gen_sp()`, `build_ranks_df()`, `score_assignments()`,
default settings, the 36 columns). The 183-unit figures reproduce the cached
`res_paired` exactly:

| Genus F1, mean of 36 columns | ITS1 | ITS2 | full |
|---|---:|---:|---:|
| 183 common units (as in chapter 07) | 0.928 | 0.920 | 0.932 |
| 150 trimmed in every region | 0.932 | 0.923 | 0.935 |

The order full > ITS1 > ITS2 holds, so the conclusion on ROADMAP line 54
stands. Per method, the 150 shift every value by +0.002 to +0.005, and one
order changes: **dada2's best region goes from full (0.948) to ITS2 (0.949)**.
Its three regions lie within 0.003 on either set, which is smaller than the
dada2 seed spread of §1.9. At Gen_sp the means barely move (ITS1 0.476,
ITS2 0.454, full 0.466 on the 150). So for Tedersoo this is a
**method-description error** (Medium): the text must say what the units are,
and the paired reading should use the 150. The all-units comparisons Q4.1–Q4.3
were not re-scored. They hold 22, 35 and 47 untrimmed units per region, a
larger share than Q4.4 does.

**In silico "full" is the ASV itself, which is already the full ITS.**
(Corrected 2026-10-06: the first version of this report said the reads begin
inside ITS1. They do not.) ITSx writes a `full` sequence only when it finds
both flanks, the end of the SSU and the start of the LSU. The in silico
templates are records of `Unite_s_all_20250219`, which UNITE has already cut
to ITS1-5.8S-ITS2. ITSx 1.1.3, rerun with the pipeline's options
(`check_insilico_full.R` in the appendix), gives:

| Sequences | SSU found | ITS1 | 5.8S | ITS2 | LSU found | in `full.fasta` |
|---|---:|---:|---:|---:|---:|---:|
| 100 UNITE templates | 0 | 99 | 99 | 98 | 1 | 0 |
| 71 in silico ASVs | 0 | 70 | 70 | 70 | 1 | 0 |

Each ASV is its template less 1 base at each end (68 of 71 lie inside a
template). Without an SSU, ITSx lets ITS1 run from base 1 to the 5.8S
(positions file: `SSU: Not found  ITS1: 1-119  5.8S: 120-276`), and without
an LSU, ITS2 runs to the last base. So the in silico ITS1 and ITS2 inputs are
true ITSx regions. The "full" input is the untrimmed ASV, but that ASV spans
ITS1-5.8S-ITS2 and nothing else. The comparison holds; only the label is
wrong. Chapter 07 should say that this input was not cut by ITSx, because its
templates are already the full ITS.

**Fix (re-render only, no reassignment).**

1. Record the trimming status once, in chapter 01 `prepare_region()`, by
   comparing sequences (no store change needed):

   ```r
   # In prepare_region(), after reading d_asv_r. `orig` comes from the DADA2
   # store, read once per dataset in prepare_dataset().
   orig <- as.character(tar_read(d_asv, store = ds$dada2_store)@refseq)
   names(orig) <- taxa_names(tar_read(d_asv, store = ds$dada2_store))
   real_seq <- tar_read_raw(paste0("d_", region, "_asv_real"), store = ds$store)
   now <- stats::setNames(as.character(real_seq@refseq), taxa_names(real_seq))
   out$trimmed <- names(now)[now != orig[names(now)]]
   ```

2. In `score_common_units()`, intersect the trimmed sets, not the unit
   names:

   ```r
   common <- Reduce(intersect, lapply(preps, `[[`, "trimmed"))
   ```

   Report the 150 next to the 183. Optionally add a second paired reading
   that keeps the untrimmed units, so the difference shows.
3. In chapter 07, add a table "units trimmed / untrimmed per region"
   beside `n-real`, and rewrite lines 186-191. Add one sentence on the
   in silico "full" input: these are the ASVs as they are, already the full
   ITS of their UNITE templates, with no flank for ITSx to cut.
4. Correct ROADMAP lines 49 and 53-54 and the matching `HISTORY.md` §2 Q4
   narrative once the paired numbers are recomputed.
5. Longer term, upstream: `itsx_pq()` could write a logical `ITSx_detected`
   column next to `ITSx_origin`. This would invalidate the ITSx targets [A],
   so do it only with the next forced rerun.

Whether `keep_undetected = TRUE` stays for the **main** grid (Q1–Q3, where
keeping every unit is defensible) is a developer decision (§11, Q1). For Q4 it
must not stay silent.

### 1.2 The external controls skip the 8 rarest phyla (`sample()` on a length-one vector). **Medium** [A]

**Where.** `make_databases.R:493-497`, `derive_fake_ref()`:

```r
one_per_phylum <- vapply(
  unique(phyla),
  function(p) sample(which(phyla == p), 1),
  integer(1)
)
```

When a phylum has a single record at position `k`, `which()` returns `k` and
`sample(k, 1)` draws from `1:k`, which is any of the first `k` records. Three
lines below, the same file avoids this trap for `remaining` (comment at
`:500-501`), but not here.

**Evidence** (replaying the draw with `targets_seed` on
`sources/Unite_s_all_20250219_wo_Fungi.fasta`, restricted to the 13 kingdoms
of `external_control_kingdoms()`): the pool holds 90 525 records in 73 phyla,
**8 of them singletons** (Brachiopoda, Cryptista_phy_Incertae_sedis,
Developea, Gromida, Hemichordata, Nematomorpha, Streptophycophyta,
Streptophyta). All 8 "one per phylum" picks for them land on a record of
another phylum. The current `fake_ref_asv_100.fasta` holds 100 records in
**65 distinct phyla**, and **all 8 singleton phyla are missing**. The replay
reproduces the file, so the file was built by this code.

**Consequence.** The function's comment (`make_databases.R:465-470`) and the
name of its test both promise "one record per phylum, then random fill". The
controls in every store cover 65 of 73 phyla, and the 8 rarest are excluded. The
effect on the rates (`ext_fungi`, `ext_correct`) is probably small. The
documented design is not what ran, though, and the Methods must not say "one
per phylum".

**A second defect in the same function.** `n_topup` is computed from
`length(one_per_phylum)` *before* `unique()` (`:502-504`). When a singleton
pick collides with another phylum's pick, the file holds fewer than `n`
records. The production file has 100, so no collision happened there.

**Why the test did not catch it.** `tests/test_make_databases_helpers.R:247-258`
uses `tiny_fasta()` (lines 14-25), whose only singleton phylum (Basidiomycota)
sits at position 3. `sample(3, 1)` picks it with probability 1/3, and with
`seed = 1` it happens to. With another seed, the test either loses a phylum or
returns 2 records instead of 3.

**Fix.**

```r
one_per_phylum <- vapply(
  unique(phyla),
  function(p) {
    idx <- which(phyla == p)
    idx[sample.int(length(idx), 1L)]
  },
  integer(1)
)
remaining <- setdiff(seq_along(seqs), one_per_phylum)
n_topup <- min(length(remaining), max(0, n - length(unique(one_per_phylum))))
```

Make the test seed-proof: loop over `seed = 1:20` and assert, for each seed,
that the output holds every phylum and `min(n, length(seqs))` records.

**Cost.** Fixing the code alone changes nothing (`derive_all_variants()`
skips an existing file). Regenerating the file (`force = TRUE`) changes
`fake_ref_file`, then every `itsx_controls_*` target, every input and every
assignment of every store, plus the `_Fungi_rep` draws, which exclude
candidates near the controls. That is a full rerun of every dataset. Fix the
code now, and decide between "document 65/73 phyla" and "regenerate before
the final runs" (§11, Q2).

### 1.3 The cost layer does not measure vsearch or BLAST+. **High** for M2 [0 code, A to remeasure]

**Where.** `R/autometric_helpers.R:18-42`, `with_autometric()`, which calls
`autometric::log_start(path, seconds)`. Its `pids` argument defaults to
`c(local = Sys.getpid())` (checked with `args(autometric::log_start)`).

**What it implies.** sintax and lca run `vsearch`, and blastn runs `blastn`
and `makeblastdb`, as **child processes** through `system2()`. autometric
samples only the R process, which sits idle waiting for the child. dada2 runs
inside R (C++), so it alone is measured.

**Evidence** (`benchmark_costs`, one row per computation):

| Store | method | median wall (min) | median peak RSS (MB) | max peak RSS (MB) | median CPU % |
|---|---|---:|---:|---:|---:|
| wood | blastn | 0.38 | 934 | 978 | 0.13 |
| wood | lca | 0.33 | 908 | 910 | 0.02 |
| wood | sintax | 1.03 | 909 | 909 | 0.008 |
| wood | dada2 | 12.3 | 18 945 | 55 007 | 47.3 |
| Tedersoo | blastn | 0.52 | 936 | 987 | 0.02 |
| Tedersoo | lca | 0.36 | 936 | 937 | 0.015 |
| Tedersoo | sintax | 0.35 | 936 | 950 | 0.017 |
| Tedersoo | dada2 | 6.7 | 18 645 | 55 837 | 44.3 |

The non-dada2 peaks are flat at ~0.9 GB **whatever the database** (from
`_Fungi` to the 1.9 M records of `EUK_ITS_v2.1`), and their CPU is ~0 %. That
is the footprint of an R session holding the pqverse checkouts.
`docs/diagnostics_cv_and_databases.md` §7 measured 5.9–6.2 GB for sintax and
lca by other means. The "1.2–1.3 GB for sintax, lca and blastn" of
`AGENTS.md:195` and `docs/experimental_design.md` §3 is therefore not a
measurement of those methods. **Wall time is valid**, because the parent waits
for the child.

**Fix, before the M2.4 final timing run** (ROADMAP §4). Two measurements are
missing, and each has a cheap fix:

- **CPU of the children, exact and free.** R already accounts the CPU of
  every child it has waited for: `proc.time()` has `user.child` and
  `sys.child`. Take the difference around `expr`.
- **Peak RSS of the process tree.** It must be sampled from **another**
  process, because the R process is blocked inside `system2()`. A small
  sampler script, launched with `processx`, polls the crew worker and all its
  descendants with `ps`.

```r
# R/sample_tree.R (run as a script): Rscript R/sample_tree.R <pid> <csv> <seconds>
# Appends "epoch,rss_mb" every <seconds> for <pid> and all its descendants,
# the sampler itself excluded (it is one of those descendants), until <pid>
# exits.
args <- commandArgs(trailingOnly = TRUE)
root <- ps::ps_handle(as.integer(args[[1]]))
path <- args[[2]]
self <- Sys.getpid()
cat("epoch,rss_mb\n", file = path)
while (ps::ps_is_running(root)) {
  kids <- tryCatch(ps::ps_children(root, recursive = TRUE), error = function(e) list())
  kids <- Filter(function(h) ps::ps_pid(h) != self, kids)
  procs <- c(list(root), kids)
  rss <- vapply(procs, function(h) {
    tryCatch(ps::ps_memory_info(h)[["rss"]], error = function(e) 0)
  }, numeric(1))
  cat(sprintf("%.0f,%.1f\n", as.numeric(Sys.time()), sum(rss) / 1e6), file = path, append = TRUE)
  Sys.sleep(as.numeric(args[[3]]))
}
```

```r
# In with_autometric(), around force(expr):
tree_csv <- sub("\\.txt$", "__tree.csv", path)
sampler <- processx::process$new(
  file.path(R.home("bin"), "Rscript"),
  c(here::here("R/sample_tree.R"), Sys.getpid(), tree_csv, seconds)
)
on.exit(sampler$kill(), add = TRUE)
cpu0 <- proc.time()
value <- force(expr)
cpu <- proc.time() - cpu0
write.csv(
  data.frame(
    phase = phase,
    cpu_self_s = cpu[["user.self"]] + cpu[["sys.self"]],
    cpu_child_s = cpu[["user.child"]] + cpu[["sys.child"]],
    wall_s = cpu[["elapsed"]]
  ),
  sub("\\.txt$", "__cpu.csv", path), row.names = FALSE
)
value
```

Tested in a scratch session, the sampler excluding itself. The tree reads
68 MB while the parent idles (the parent alone: 69 MB, which is all
autometric records today). It peaks at 730 MB while a child `Rscript` fills a
400 MB vector.
`summarise_autometric_costs()` then reads `max(rss_mb)` from the
`__tree.csv` and the CPU-seconds from the `__cpu.csv`. `elapsed` replaces
the first-to-last-sample span, so sub-second computations no longer get NA
costs. A lighter alternative, if MiscMetabar can prefix its external calls:
`/usr/bin/time -f "%M %e %U %S"` around `vsearch` and `blastn` (max RSS in kB,
elapsed, user and system seconds).

Apply the S1.6a guard (`autometric::log_stop()` before `log_start()`,
`R/autometric_helpers.R:24-37`) in the same commit. Both change the code of
a function every compute target depends on, so put them in the commit that
precedes the M2.4 run, as that comment already plans.

### 1.4 CO₂: every computation is charged the machine's whole RAM. **High** for M2 [R]

**Where.** `analysis/05_m2_costs.qmd:155-170`, `add_co2()`. It calls
`greenAlgoR::ga_footprint()` without `memory_ram`, and `ga_footprint()`
(greenAlgoR `R/ga_footprint.R:262-264`) then uses
`benchmarkme::get_ram()`, the **total RAM of the machine** (62 GB here).

**Magnitude.** Memory power is 62 GB × 0.3725 W/GB = 23.1 W for every
computation, against 12 W × `n_cores` of CPU (12 W for every 1-thread
sintax, lca or blastn job, 36 W for dada2). So for the fast methods, about
two thirds of the reported CO₂ comes from RAM they never used. The ranking
between methods is compressed toward wall time alone, and the absolute
values are inflated by ~3×.

**Fix.** Pass the measured peak in GB, once §1.3 makes it real, or a
documented constant if it stays unmeasured:

```r
co2eq_gCO2 = greenAlgoR::ga_footprint(
  runtime_h     = wall_time_s / 3600,
  location_code = "WORLD",
  n_cores       = if (method == "dada2") assign_threads_dada2 else assign_threads_fast,
  memory_ram    = peak_tree_rss_mb / 1000,
  add_ref_values = FALSE
)$carbon_footprint_total_gCO2
```

Say in the Methods which memory figure enters the formula. Until §1.3 is
fixed, `peak_resident_mb` is wrong for three methods out of four, so using
it would only swap one error for another.

### 1.5 Q2.5: the wall time per database is a blastn-weighted mean, split in two points per database. **Medium** [R]

**Where.** `analysis/03_q2_databases.qmd:211-218`:

```r
costs_for(v$dataset, res) |>
  filter(phase %in% res$method_db, !startsWith(db, "mini_")) |>
  mutate(db_size_mb = file.size(here(db_path)) / 1e6) |>
  group_by(db, db_size_mb) |>
  summarise(mean_wall_h = mean(wall_time_s, na.rm = TRUE) / 3600, ...)
```

**Two defects.**

1. `benchmark_costs` has one row per **benchmark row** (21 per database and
   input), and the rows of one computation share its cost
   (`pipelines/assign_taxo.R:342-349`). The mean is therefore weighted 12/21
   blastn, 3/21 each for dada2, sintax and lca, so "mean wall time" mostly
   tracks blastn.
2. `db_path` is the dada2-format file for dada2 rows and the sintax-format
   file for the others, so `group_by(db, db_size_mb)` gives **two points per
   database**. One is dada2 alone at the dada2-file size. The other is the
   blastn-weighted mean of the three other methods at the sintax-file size.

**Fix.**

```r
db_size_tbl <- costs_for(v$dataset, res) |>
  filter(phase %in% res$method_db) |>
  distinct(compute_name, method, db, wall_time_s) |>
  left_join(
    tibble(db = db_list,
           db_size_mb = file.size(here(db_path_for("sintax", db_list))) / 1e6),
    by = "db"
  )
```

Then plot one point per method × database (`colour = method`). Or sum the
four computations per database, if the question is "what does this database
cost for the four methods".

### 1.6 `EUK_ITS_v2.1_Fungi_cut` still holds the records OD 32 removed, and the zoom database is not built at all. **High** for item 0.11b [0 now]

**Where.** `make_databases.R:53-68` (`skip_if_exists()`), `:784-811` (the
length cap applies to the whole releases and to `_Fungi`), `:812-827` (the
`_cut` variants).

**Evidence.**

| File (`sintax_format/`) | built | records > 5 000 bp | longest |
|---|---|---:|---:|
| `EUK_ITS_v2.1_Fungi.fasta` | 2026-09-24 (after OD 32) | 0 | 4 984 |
| `EUK_ITS_v2.1_Fungi_cut.fasta` | **2026-09-15** (before OD 32) | **11** | **32 568** |

The 32 568 bp record is the Laccaria genomic clone that OD 32 removed because
it "made sintax name every shuffled control on `EUK_ITS_v2.1_Fungi`" (SA5).
The `_cut` file was trimmed on 2026-09-15 from the uncapped `_Fungi` file.
`derive_all_variants()` never rebuilds a file that exists, and
`cap_record_length()` is never applied to the `_cut` variants. Trimming does
not save the long records, because `cut_discard_untrimmed = FALSE` keeps the
records without primer sites at full length.

**What it means for the zoom (0.11b).** The planned zoom database is
`fungi + rep cut`, that is, the `_Fungi_cut` records plus untrimmed
representatives (`docs/objectives_design.md` §1.3, decision 16), on one UNITE
and one EUKARYOME database. Today:

- the only `_Fungi_cut` file, the EUKARYOME one, which would be its fungal
  part, is stale;
- **nothing builds the zoom database itself**: `make_databases.R` has no
  `_Fungi_rep_cut` path;
- `config.R::zoom_dbs` lists no UNITE `_Fungi_cut`.

Item 0.11b is therefore bigger than "moderate": it needs a derivation
function, a cap, and a `zoom_dbs` row per source (§11, Q8).

**The root cause is general.** In `make_databases.R`, a file that exists
counts as up to date. Three consequences:

- a change of rule (a new cap, a header fix, a new `rep_share`) never reaches
  the files already derived unless someone remembers `force = TRUE`, file by
  file;
- an interrupted derivation leaves a **truncated** output, which the next run
  skips silently. Every `derive_*` writes straight to its final path:
  `dbpq::format2sintax(output_path = output)`, `filter_db(output = ...)`, the
  `awk ... > output` lines;
- `force = TRUE` deletes the old file **before** rebuilding (`:54-56`), so a
  failed rebuild loses the good copy.

**Fix.**

1. Now: rebuild the two `_Fungi_cut` files (`derive_cutadapted(...,
   force = TRUE)` for `dada2` and `sintax`), then cap them like the others.
   In `derive_all_variants()`, the `cut_sources` loop should end with:

   ```r
   if (s %in% names(max_record_length)) {
     cap_record_length(refseq_path(paste0(s, "_Fungi_cut"), fmt), max_record_length[[s]])
   }
   ```

2. Make the skip rule a staleness rule, and the writes atomic:

   ```r
   # Rebuild when the output is missing or older than one of its inputs.
   is_stale <- function(output, inputs) {
     !file.exists(output) || any(file.mtime(inputs) > file.mtime(output))
   }
   # Write to output.part, then rename: an interrupted run leaves no output.
   with_atomic_output <- function(output, write_fun) {
     part <- paste0(output, ".part")
     unlink(part)
     write_fun(part)
     if (!file.rename(part, output)) stop("Could not move ", part, " to ", output)
     normalizePath(output)
   }
   ```

   `config.R` constants are inputs too: compare against
   `file.mtime(here("config.R"))`, or better, against a small stamp file per
   rule (the cap, the header rules, `rep_share`).
3. The structural fix is a `{targets}` project for the databases
   (`format = "file"` targets whose commands read the constants). The
   invalidation logic this file reimplements by hand would then come for free.
   About a day of work; worth it only if more releases are coming.

**Dead derived files** (no grid, no zoom, no dataset uses them):
`EUK_SSU_v2.1*` (SSU dropped on 2026-09-15) and `*_Fungi_cut_full_ITS*`
(`bio_datasets$cut_suffix` is `NA` everywhere). That is about 3.6 GB in
`dada2_format/` and `sintax_format/`.

### 1.7 Every current score depends on uncommitted code, and nothing records it. **High** (reproducibility) [0]

**Evidence.**

- **comparpq** (`pqverse_pkg/comparpq`, HEAD `8535d66`, 2026-10-01) has 4
  modified files. The diff of `R/compare_taxo_unit.R` (+26 lines) **is
  OD 36**: the `foreign` column, `n_foreign`, `foreign_named`, and
  `scored_real <- is_real & !deeper_than_truth & !is_foreign`. Every number in
  `res_comp_tax` since 2026-10-02 was produced by this uncommitted code,
  including the Tedersoo Genus F1 that went from 0.52–0.59 to 0.89–0.93.
- **This project's** last commit is `3492332` (2026-09-22, "A new step").
  Since then, all of OD 34, 35 and 36, the patched dada2, three datasets,
  chapters 07–08 and 4 test files are uncommitted. That is 40 modified files
  and these untracked ones: `R/assignment_helpers.R`, `R/bio_loaders.R`,
  `R/dada2_patched.R`, `R/insilico_input.R`, `R/itsx_inputs.R`,
  `pipelines/load_bio.R`, `analysis/07_q4_its_region.qmd`,
  `analysis/08_q5_units.qmd`, `patches/`, `docs/manuscript_framing.md`,
  `docs/parked.md`, `tests/test_{assignment_helpers,bio_loaders,insilico_input,itsx_inputs}.R`,
  and three metadata CSVs.
- **`pqverse_versions()`** (`R/load_pqverse.R:29-48`) records
  `git rev-parse HEAD` only. A dirty checkout reports the same commit as a
  clean one. The claim "commits are pinned by the `session_info` targets"
  (`AGENTS.md:54`, ROADMAP C22-N1) is false whenever a checkout is dirty, and
  comparpq is dirty today.
- **The scoring happens in Quarto, not in a target.** No `session_info`
  covers chapter 01 at all, and no chapter prints or caches a version
  (`rg "session|pqverse_versions" analysis/` finds nothing).

**Fix.**

1. Commit comparpq and this project now (the developer commits).
2. Record the dirty state:

   ```r
   # In pqverse_versions(), inside the lapply:
   status <- tryCatch(
     system2("git", c("-C", shQuote(path), "status", "--porcelain"), stdout = TRUE),
     error = function(e) NA_character_
   )
   diff_hash <- tryCatch({
     d <- system2("git", c("-C", shQuote(path), "diff", "HEAD"), stdout = TRUE)
     if (length(d) == 0) NA_character_ else rlang::hash(d)
   }, error = function(e) NA_character_)
   # ... then add to the tibble:
   #   dirty = length(status) > 0 && !all(is.na(status)), diff_hash = diff_hash
   ```

3. In chapter 01, before `cache_write()`:

   ```r
   attr(res_comp_tax, "versions") <- list(
     pqverse = pqverse_versions(),
     project = system2("git", c("rev-parse", "HEAD"), stdout = TRUE),
     project_dirty = length(system2("git", c("status", "--porcelain", "--", "R", "analysis", "config.R"), stdout = TRUE)) > 0,
     rendered = Sys.time()
   )
   ```

   Print that attribute at the end of every chapter (one line in
   `00_setup.R`), so each HTML says which code produced it.
4. Optionally, `load_pqverse(strict = TRUE)` stops on a dirty checkout, and
   the final renders use it.

**Cost.** Changing `pqverse_versions()` invalidates only the `session_info`
targets, which are cheap.

### 1.8 The stores keep hundreds of orphan objects, and a sandbox script reads one. **Medium** [0]

**Evidence** (`tar_manifest()` against `tar_meta()`, all projects up to date):

| Project | targets in the pipeline | meta rows | orphan rows | orphan objects on disk |
|---|---:|---:|---:|---:|
| `assign_taxo` | 484 | 1 198 | 714 | 464 (101 MB) |
| `assign_taxo_korhonen_wood` | 253 | 738 | 485 | 385 (268 MB) |
| `assign_taxo_tedersoo_illumina` | 1 398 | 1 493 | 95 | 0 |

`analysis/sandbox/check_controls_by_simplification.R:26` runs
`tar_read(d_all_taxo, store = "store_assign_taxo")`. That object is an orphan
built on **2026-09-22 15:58**: before OD 34, before the patched dada2, and
before the OD 32 and OD 35 caps. The script is cited as evidence in
`HISTORY.md:828`, and it is in the modified set, so it is still in use.
`compare_lca_stampa_vs_common.R:89`, `check_threads_determinism.R:20`,
`timing_query_x_refsize.R:144` and `timing_by_query_class.R:32` name
`d_asv_for_assignation`, which is an orphan too.

**Cause.** `make.R` prunes through `run_project()`, but production now runs
through `run_bio_night.R`, which never prunes (§3.2).

**Fix.**

- In `run_bio_night.R`, after a successful `tar_make()`
  (`status == "ok"` and `length(still) == 0`), call `targets::tar_prune()`.
  A script that reads a removed target then fails loudly instead of reading
  three-week-old data.
- Port the five sandbox scripts to the OD 34 names (`d_all_taxo_ITS1_asv`,
  `d_ITS1_asv_for_assignation`), or move them to `Archives/`.
- Expect one large diff of the tracked `store_*/meta/meta` files the first
  time.

### 1.9 dada2 depends on its seed, and is not exactly reproducible even with it. **Medium** for the co-primary PPV endpoints [0 to measure, A to patch]

Two effects, which must not be confused.

**(a) Bootstrap variance: the seed matters.** dada2 draws its bootstrap k-mers
with `Rcpp::runif()` (dada2 1.40.0 `src/taxonomy.cpp:283-284`), so they follow
`set.seed()`. `targets` gives each target its own seed
(`tar_seed_create(name, global_seed = 22)`), so every dada2 column of every
store is **one draw** from the bootstrap distribution. A fresh run of
`ITS1_asv_compute_dada2__Unite_all_20250219_Fungi` (Pauvert,
`d_ITS1_asv_for_assignation`, 392 units) with **another seed** (22) gave the
same best-hit Genus everywhere, but other bootstraps on 302 of 392 units.
Scored against the Pauvert truth (real units only, hleap §6 rules):

| `min_bootstrap` | rank | TP stored / other seed | FP stored / other seed | PPV stored / other seed |
|---|---|---|---|---|
| 0.4 | Genus | 134 / 133 | 12 / 12 | 0.918 / 0.917 |
| 0.5 | Genus | 133 / 133 | 9 / 11 | **0.937 / 0.924** |
| 0.6 | Genus | 133 / 133 | 9 / 10 | 0.937 / 0.930 |
| 0.5 | Species | 72 / 72 | 35 / 36 | 0.673 / 0.667 |
| 0.6 | Species | 71 / 70 | 33 / 32 | 0.683 / 0.686 |

F1 hardly moves, because a changed call swings between a wrong name and NA,
which moves a unit between FP and FN. **PPV moves by up to 0.013** between two
seeds on one combination. PPV at Genus and Species is a co-primary endpoint
(OD 37), and several effects the chapters read are of the same order: the Q4
regions are "within 0.006–0.021 per method", and the Q1 bootstrap word
threshold is 0.03. sintax has the same kind of variance, and a larger one:
`config.R:345-349` records that two unseeded sintax runs against
`EUK_ITS_v2.1_Fungi` "gave a different genus to 11 of 71 in silico ASVs".
`sintax_seed` (SA6) made that reproducible, but did not make it smaller. A
single seed hides this variance for both k-mer methods.

**(b) Non-reproducibility at a fixed seed: small.** Ties between equally
scored genera are broken by an engine seeded from `std::random_device`
(`src/taxonomy.cpp:80-81`, in `get_best_genus()`, the function the overflow
patch already edits), which no R seed reaches. Evidence:

- On 30 `data_fungi_mini` ASVs, two `assignTaxonomy()` calls after the same
  `set.seed(1)` return different bootstrap matrices.
- Rerunning the same Pauvert target **with its own target seed** and its
  3 threads, the bootstraps differ on 44 of 392 units, by at most 9 points.
  The Genus calls are identical at 0.4, 0.5 and 0.6, 1–2 Species calls of
  real units change, and PPV is unchanged except at Species 0.6
  (0.683 → 0.689).

So the stores are reproducible to within one or two Species calls. The
`session_info` cannot promise bit-identity for dada2, and the Methods should
say so in one sentence.

**Fix.**

1. **Measure (a) now [0, about 30 min of compute].** Run 3 extra seeds × 3
   databases on Pauvert, for dada2 and sintax, outside the stores (a sandbox
   script calling `compute_assignment()` with
   `set.seed(tar_seed_create(name, s))` for `s` in 23–25). Report the seed
   spread of Genus and Species PPV next to every dada2 or sintax comparison.
   If the spread is comparable to an effect, the effect is not a result.
2. **Fix (b) only if bit-identity matters [A, with M2.4].** A second patch
   that seeds the tie-break engine from R. `get_best_genus()` runs inside
   `RcppParallel`, so draw the seeds on the main thread next to `unifs`
   (`:283-286`) and pass them down:

   ```cpp
   // In C_assign_taxonomy2, beside the bootstrap uniforms:
   Rcpp::NumericVector tie_seeds = Rcpp::runif(nseq * (NBOOT + 1));
   // ... pass the relevant element down to get_best_genus(), then, in it:
   std::mt19937 gen((uint32_t) (tie_seed * 4294967295.0));
   ```

   Offer it upstream with the overflow PR (same function).

### 1.10 A tracked input differs from the one that was used. **Medium** [0]

`data/data_raw/insilico/iss_input_fungi.fasta` sits under a gitignored folder
(`.gitignore:220`) but is tracked (commit `56b1439`). The working copy, the
second generation used by `store_dada2_insilico`, differs from HEAD by
+760/−822 lines. A clone would hand InSilicoSeq the archived first
generation. Commit it. The same class covers `Rplots.pdf` and
`run_assign_taxo.log`: both are tracked although they are stray outputs (§5).

---

## 2. Latent defects, verified absent from the current stores. **Low**

Each of these is a real defect in the code, but a check on every current
store found **no affected case**. They are listed so that nobody
re-investigates them, and so they are fixed before a new dataset or release
triggers them.

| # | Where | Defect | Verified on the stores | Fix |
|---|---|---|---|---|
| L1 | `R/assign_compute.R:121-126` | `derive_dada2_row()` keeps the ranks whose bootstrap passes and **packs them to the left**, as `dada2::assignTaxonomy()` itself does (`taxes[[i]][boots[i, ] >= minBoot]`, checked in the 1.40.0 source). With a non-monotone bootstrap (a homonym at a deep rank under two lineages, which the EUKARYOME header fix makes more likely by stripping the `(Fungi)` qualifiers), a Genus would land in the Class column. | 0 non-monotone bootstrap rows in all `*_compute_dada2__*` targets of the 5 stores | Keep the dada2 semantics but **stop at the first failing rank**: `k <- cumprod(computed$boot[i, ] >= 100 * min_bootstrap) == 1`. This differs from dada2 only on the cases dada2 gets wrong; say so in the Methods. [A] |
| L2 | `R/assign_compute.R:139` | `taxo[computed$taxo_bootstrap < min_bootstrap] <- NA` also compares the character column `taxa_names` with a number (`"Taxa_1" < 0.4` is a string comparison, `FALSE` because `"T" > "0"`). It works by accident. | structure checked (`taxo_bootstrap$taxa_names` is `chr`) | `rank_cols <- setdiff(names(taxo), "taxa_names"); taxo[rank_cols][computed$taxo_bootstrap[rank_cols] < min_bootstrap] <- NA` [A] |
| L3 | `R/assignment_helpers.R:72-77` | `gen_sp_for()` decides **per column** whether Species already carries the genus (all non-NA values hold `_`), so one value without `_` makes the whole column `Genus_Genus_species`. And `out[grepl("NA", out)] <- NA` erases any name that contains the letters `NA`. | 0 mixed columns and 0 `NA`-substring names over 189 Species columns × 4 inputs | Decide per value: `out <- ifelse(is.na(sp), NA, ifelse(grepl("_", sp), sp, paste(gen, sp, sep = "_"))); out[is.na(gen) & !grepl("_", sp)] <- NA` [R] |
| L4 | `R/mock_truth.R:392-399` | `mock_truth_table()` stops the truth at the first rank where the GBIF lineage is `NA`, **even when Genus and Species are known**. A strain whose order or family GBIF leaves empty (common for fungi incertae sedis) would lose its Genus and Species truth. | only `Sordariomycetes_sp` (Pauvert) has an NA above Genus, and that stop is intended | Treat an NA above Genus as "not scored at that rank" (needs a comparpq option), or fill it from the Table S1 lineage, and write the rule in hleap §6. [R] |
| L5 | `R/mock_truth.R:280-303`, `make_databases.R:560-581` | vsearch truncates labels at the first space, and the results are matched back on the full name. | Sanger IDs: 0 with a space (Pauvert 189, Tedersoo 103). **2026-10-06: hit on the Rolling mock (OD 42)**, whose strain names had spaces: every *Candida* strain came back as "Candida" and no ASV got a truth; worked around by underscores in the truth csv, the code still truncates. `sintax_format/EUK_ITS_v2.1.fasta`: **3** headers with a space (Peronospora, e.g. `s:atriplicis.hastatae minor`), and they are Stramenopila candidates for the `_Fungi_rep` draw, so `records_near_external()` could miss them | Add `"--notrunclabels"` to both `system2()` calls. [0 for the truth; A for `_Fungi_rep` only if rebuilt] |
| L6 | `R/bio_loaders.R:124, 191` | `tax_rows[match(keep, ...), ]` and `taxo[match(keep, taxo$OTU), ]` give NA rows, hence NA taxon names, for an OTU missing from the taxonomy table. Today the four OTUs missing from Korhonen's taxonomy are present in no sample, so they never reach `keep`. | wood and air stores build | `stopifnot(!anyNA(match(keep, taxo$OTU)))` after each `match()`. [0] |
| L7 | `R/insilico_input.R:55`, `make_databases.R:482, 690` | `set.seed()` inside a function changes the caller's RNG stream. | — | `withr::local_seed(seed)` [0] |
| L8 | `R/insilico_input.R:61-62` | When `n` cuts a round mid-turn, the phyla that sort first alphabetically get the extra records: a deterministic bias toward A–C phyla. | — | Shuffle the phylum order within each turn: `order(pool$turn, sample(nrow(pool)))`. [0; the file exists] |

---

## 3. Reproducibility and pipeline structure

### 3.1 The smoke-test project writes its costs into the production folder. **Medium** [0]

**Where.** `pipelines/assign_taxo.R:77-80`:

```r
autometric_dir_assign <- here(
  "data/data_final/autometric",
  if (is.null(bio_cfg)) "assign_taxo" else paste0("assign_taxo_", bio_dataset)
)
```

For `TAR_PROJECT = "assign_taxo_mini"`, `config.R:216-225` gives
`bio_dataset = "mini"`, which is not a row of `bio_datasets`, so `bio_cfg` is
`NULL` and the folder is **`assign_taxo`, the production one**. The phase
names are the same in both projects (`compute_name` ignores `mini_db`). So a
smoke test leaves the newest log of each phase in the production folder, and
`summarise_autometric_costs()` keeps the newest log. This is the same failure
the in silico run caused on 2026-09-22, which the comment at `:70-76`
describes.

**Evidence.** Each log file was matched to the target of the same name whose
recorded build ended within 2 min after the file's stamp. **168 of the 1 004
files in `data/data_final/autometric/assign_taxo/` match the last build of
their phase in `store_assign_taxo_mini`**, and none matches a production
build. That is a lower bound: the meta keeps only the last build of each
target. An example is `compute_sintax__EUK_ITS_v2.1_Fungi__20260924_145940.txt`,
built by the mini store at 2026-09-24 14:59:41. Exactly one of the 168 is
still the newest log of its phase, and that phase is an orphan. The OD 34
production computations (2026-09-30 to 10-01) are newer than the last mini
run (2026-09-27), so **no current cost is wrong**. A smoke test run between
the M2.4 timing run and the last rebuild of `benchmark_costs` would make it
wrong.

**Fix** (no invalidation: the value stays `.../assign_taxo` for production and
`.../assign_taxo_<dataset>` for the datasets):

```r
autometric_dir_assign <- here(
  "data/data_final/autometric",
  Sys.getenv("TAR_PROJECT", unset = "assign_taxo")
)
```

Move the 168 identified mini logs out of the production folder, using the
rule above: same phase name, and a mini build `time` within
`[stamp, stamp + seconds + 120 s]`. Do not select them by date range,
because production runs happened on the same days.

### 3.2 `make.R` is no longer the entry point, and nothing rebuilds the dataset stores. **Medium** [0]

`make.R` runs `dada2`, `assign_taxo` and `cross_val`, the last of which is
parked (OD 28). It runs **none** of the eight dataset projects. Production
actually runs through `run_bio_night.R`, which runs `assign_taxo_<dataset>`
only. It does not run `dada2_<dataset>` or `load_bio.R` first, and it does not
prune (§1.8). `run_insilico_night.R` duplicates
`Rscript run_bio_night.R insilico`, and its header is stale ("48 ASVs", "same
concurrency as … 2026-09-22"). A fresh clone has no single command that
rebuilds what the paper reads.

**Fix.** One function, one table of projects, used by both `make.R` and the
night script:

```r
# R/run_project.R
production_projects <- function() {
  ds <- bio_datasets$dataset
  c("dada2", "assign_taxo", as.vector(rbind(paste0("dada2_", ds), paste0("assign_taxo_", ds))))
}
run_projects <- function(projects = production_projects(), prune = TRUE) {
  for (p in projects) {
    withr::with_envvar(c(TAR_PROJECT = p), {
      targets::tar_make(reporter = "timestamp")
      if (prune && length(targets::tar_outdated(reporter = "silent")) == 0) {
        targets::tar_prune()
      }
    })
  }
}
```

`make.R` becomes `derive_all_variants(); run_projects()`. The night script
calls `run_projects(projects = <its arguments>)`. Delete
`run_insilico_night.R`. Keep `cross_val` out of `make.R` while OD 28 is
parked.

### 3.3 The Pauvert mock is a special case in eight places. **Medium** [0 / A]

The production dataset is "no `bio_cfg`". That gives `is.null(bio_cfg)`
branches at `config.R:232, 245, 281, 286` and `pipelines/assign_taxo.R:79`, a
second copy of the dataset → regions/units mapping in
`analysis/00_setup.R:57-93` (`dataset == "pauvert"`), a third in
`run_bio_night.R:33-37`, and **two DADA2 pipelines**. `pipelines/dada2.R`
(271 lines) is an 85 % copy of `pipelines/dada2_bio.R` (388 lines). A `diff`
of their `tar_plan()` shows only paths, fastq patterns and the misnamed
`d_idtaxa` branch (a `clusterize` clustering). The `"mini"` bug of §3.1 comes
from exactly this special case.

**Fix, in two steps.**

- **[0] now.** `00_setup.R` must not re-derive regions and units. Have
  `config.R` expose one function, `dataset_config(dataset)`, that returns
  `itsx_regions`, `itsx_units`, the stores and the autometric folder for any
  dataset, `"pauvert"` included. `config.R` itself and `00_setup.R` both call
  it.
- **[A] only if the harness is meant to be reused (§11, Q5).** Make
  `"pauvert"` a row of `bio_datasets` (`raw_dir = "data/data_raw/rawseq"`,
  `sam_data = "sam_data.csv"`, the ITS1F/ITS2 primers, `name_strip = "_R.*"`,
  `itsx_region = "ITS1"`, `otu_input = TRUE`). Point the `dada2` project of
  `_targets.yaml` at `dada2_bio.R`, keep the store names, and delete
  `pipelines/dada2.R`. The `dada2` store would rebuild (about 4 min). If
  `d_asv` comes out byte-identical, nothing downstream reruns. Check that with
  `tar_outdated()` on a copy of the store first.

### 3.4 Guards that run in the wrong process, or not at all. **Low** [0]

- `assert_dada2_patched()` runs in the **main** process
  (`pipelines/assign_taxo.R:30`, `pipelines/cross_val.R:33`), but dada2 runs
  in the **crew worker**. The worker gets the patched build only because crew
  1.3.2 starts it with `r_arguments = c("--no-save", "--no-restore")` (no
  `--vanilla`) in the project root, where `.Rprofile` applies. That is
  correct today, and checked nowhere. Fix without invalidation: a new,
  independent target on the worker,

  ```r
  tar_target(
    dada2_build_on_worker,
    {
      p <- utils::packageDescription("dada2")$Patched
      if (is.null(p) || !grepl("taxonomy_overflow", p)) stop("worker loads ", find.package("dada2"))
      find.package("dada2")
    },
    resources = if (single_worker) tar_resources() else tar_resources(crew = tar_resources_crew(controller = "dada2_ctrl"))
  )
  ```

- `pipelines/dada2.R` and `pipelines/dada2_bio.R` do not call
  `assert_dada2_patched()` at all. Their seed taxonomy (`tax_tab`) runs
  `assignTaxonomy()` on `Unite_s_all_20250219_Fungi`. It is only a
  placeholder tax_table, but the guard costs one line.
- `pipelines/load_bio.R` has no `session_info` target, against
  `AGENTS.md:116` ("Every store also has a `session_info` target").
- `patches/install_dada2_patched.R:17-29` downloads
  `bioconductor.org/packages/3.23/bioc/src/contrib/dada2_1.40.0.tar.gz` with
  no checksum. Bioconductor moves a superseded patch release to
  `src/contrib/Archive/dada2/`, so this URL breaks at the next dada2 1.40.x.
  Record `tools::md5sum(tarball)` in the script and stop on a mismatch, fall
  back to the `Archive/` URL, or keep the tarball in `patches/`.
- `pipelines/load_bio.R:62-78` writes `selected_units_<dataset>.csv` into
  `data/data_raw/metadata/`, so a pipeline writes into raw data. Write it to
  `data/data_final/`.

### 3.5 Design constants that rely on a package default. **Low** [A at the next forced rerun]

`docs/objectives_design.md:341` lists `maxrejects = 32` among the lca settings
of OD 33, but `config.R:324-326` does not set it, and `compute_assignment()`
does not pass it. The 32 is MiscMetabar's default
(`assign_vsearch_lca(maxrejects = 32)`, `R/vsearch.R:1693`). A changed default
in the checkout would silently change the design. Add
`lca_maxrejects <- 32L` to `config.R` and pass it. The same reasoning applies
to `blast_pq()` and `assign_sintax()` defaults the design depends on (word
size, e-value): list the ones the Methods quote and pin them.

---

## 4. Analysis layer

### 4.1 Win matrix and "best" phrases break ties by row order. **Medium** [R]

`analysis/02_q1_methods.qmd:193-199` (Q1.6 win matrix) uses
`slice_max(values, n = 1, with_ties = FALSE)`. So do `best_phrase()`
(`00_setup.R:375`), `04_q3_consensus.qmd:42` and
`05_m2_costs.qmd:61, 65, 69`. Ties are common. On the default-settings ASV
rows of `res_comp_tax` (36 dataset × region × database cells per rank), the
best F1 is shared by two or more methods in **16 of 36 cells at Class, 18 at
Order, 18 at Family, 8 at Phylum**, 4 at Genus and 1 at Gen_sp. At the middle
ranks, half of the "wins" go to whichever method comes first in the rows.

**Fix.** `with_ties = TRUE`, then count a shared win as `1 / n_tied`, or show
ties as their own column. In `best_phrase()`, name every tied winner
(`"dada2 = sintax"`).

### 4.2 Warnings are hidden in every chapter. **Medium** [R]

`analysis/_quarto.yml:22` sets `warning: false` for the whole project. That
hides lme4's singular fits (Q1.7, 9 databases as random effect), ggplot's
"Removed n rows", `left_join()` many-to-many warnings and `summarise()`
regrouping. These are exactly the signals that a filter or a join went wrong.
Turn warnings on in chapter 01 at least (`#| warning: true` in its
`setup` chunk, or drop the project-wide key). Silence the known ones locally
with `suppressWarnings()` around the call, with a comment.

### 4.3 Small defects in the chapters. **Low** [R]

- **Chapter 06 renders stale results.** `_quarto.yml:15` renders `0*.qmd`, so
  `06_cross_validation.qmd` reads `store_cross_val` (7-database grid,
  2026-09-15/16, parked by OD 28) into the published HTML as if it were
  current. Exclude it (`- "!06_cross_validation.qmd"`) or open it with a
  banner.
- **A check that checks nothing.** `01_load_and_clean.qmd:648-655`
  (`check-cells`) states that the four cells sum to a constant per group, then
  only prints a table. Make it fail:

  ```r
  cells <- res_comp_tax |>
    filter(metrics %in% c("TP", "FP", "FN", "TN")) |>
    group_by(dataset, region, unit, external_scoring, tax_level, method_db) |>
    summarise(n_units = sum(values), .groups = "drop")
  stopifnot(all(
    count(distinct(cells, dataset, region, unit, external_scoring, tax_level, n_units),
          dataset, region, unit, external_scoring, tax_level)$n == 1
  ))
  ```

  On the current `analysis/_cache/res_comp_tax.rds` it passes: 168 groups,
  0 violations.

- **"Default settings" are still defined twice**: `filter_default_settings()`
  (`00_setup.R:201-209`) and `default_names_for()` (`01:105-114`). This was
  2026-09-22 S3. Replace both with one function over `name_map`:

  ```r
  default_full_names <- function(vm) {
    vm$full_name[which(
      (vm$method %in% c("dada2", "sintax") & vm$min_bootstrap == 0.5) |
        (vm$method == "lca" & vm$lca_cutoff == 0.51) |
        (vm$method == "blastn" & vm$vote_algorithm == "rel_majority" & vm$min_id == 95)
    )]
  }
  filter_default_settings <- function(res) {
    dplyr::filter(res, method_db %in% default_full_names(name_map))
  }
  ```

- **`assignment_db()`** (`R/score_assignments.R:26-30`) still parses the
  database out of the column name with two regexes, although `name_map`
  carries `db`. Join instead, and let consensus columns get `NA` from the
  failed join.
- **`default_view()`** (`00_setup.R:280-282`) defaults `keep_dataset` to
  `res$dataset[[1]]`, so the default depends on row order. Make it a required
  argument.
- **`costs_for()`** (`00_setup.R:402-413`) stops only when **no** name
  matches. A half-stale store passes silently. Stop when
  `!all(default_full_names(...) %in% costs$phase)`.
- **Q3.2's caption** (`04_q3_consensus.qmd:142`) says "mean ± SE across DBs",
  but a consensus has one value per view, so it has no error bar.
- **Q1.7** (`02_q1_methods.qmd:221-256`). `tryCatch()` turns an lme4 error
  into a row reading "ERROR". The LRT treats the 9 databases as independent
  replicates, while they are 3 sources × 3 nested filters (2026-09-22 T1).
  Either drop the p-values, or replace them with the paired unit bootstrap of
  2026-09-22 P11.

---

## 5. Dead code, duplication and stray files (C22-S5)

**Dead code** [0, unless marked]:

| Where | What | Action |
|---|---|---|
| `pipelines/assign_taxo.R:171-201`, `config.R:238-249`, `bio_datasets$taxo_mock` | `file_taxo_mock` / `taxo_mock` targets: the set-membership truth replaced by `R/mock_truth.R`. No chapter reads them. The comment at `:243-244` ("Tedersoo's own truth … still to be fetched") is stale | delete the targets, the constant and the column, plus `data/data_raw/metadata/taxo_mock*.csv` if nothing else reads them |
| `pipelines/dada2.R:218-231`, `pipelines/dada2_bio.R:327-338` | `d_vs`, `d_idtaxa` (a `clusterize` clustering), `d_*_mumu`: no consumer since OD 34 (assign_taxo clusters the ITSx ASVs itself). `dada2_bio.R:6-8` ("`d_vs` stays, since Q5 rests on the vsearch OTUs") and `00_setup.R:55` ("`d_asv` and `d_vs`") are stale | keep `track_df` rows only if the Methods use them, else delete |
| `config.R:9-11` | `fw_primer_AM`, `rev_primer_AM` (Glomeromycota primers, unused) | delete |
| `config.R:188-193`, `make_databases.R:829-862` | `bio_datasets$cut_suffix` is `NA` in every row, so the 2b-bis loop never runs | delete the column and the loop (they are in git if the zoom ever needs them) |
| `R/values_map.R:138-139, 144` | columns `cutadapted_db`, `db_filter`, `do_clean_pq`: read by nothing (`clean_pq = FALSE` is hard-coded in `compute_assignment()`) | delete (`tar_eval()` does not read them, so nothing is invalidated) |
| `R/cross_val.R:13-32, 314-325` | `cross_val_param()`, the French TODO, the `dada2_2steps` branch that `stop()`s | see §8 |
| `R/examples_cross_val.R` (134 lines) | sourced by nothing; 2026-09-22 S5, still open | `Archives/` |
| `run_insilico_night.R` | duplicate of `run_bio_night.R insilico` | delete (§3.2) |
| `R/assign_compute.R:82-84` | `derive_assignment()` returns `physeq` unchanged when `computed` is `NULL`; only `combine_taxo_assignments()` catches it later, without naming the method × database × input (2026-09-22 R2, still open) | `warning("No result for ", suffix)` [A: changes the function] |

**The sandbox has become a second code base.** `analysis/sandbox/` holds 37
files and 3 995 lines, as much as `R/` (2 909) and `pipelines/` (1 326)
together. HISTORY cites 19 of them as evidence, and five read orphan targets
(§1.8). Each script that HISTORY cites should start with one line naming what
it reads and where its result is cited. Move the others to `Archives/`.

**Stray files at the root:**

| File | State | Action |
|---|---|---|
| `Rplots.pdf` | **tracked**, stray device output | `git rm` |
| `run_assign_taxo.log` | **tracked** (committed before `run_*.log` was ignored) | `git rm --cached` |
| `run_bio_weekend_nohup.out` (345 KB) | untracked, **not ignored** | add `*.out` to `.gitignore`, move to `logs/` |
| `length_cap_test_dada2.log` | untracked, **not ignored** (the pattern is `run_*.log`) | ignore `*.log` |
| `benchmark_grid.html` (18 KB) | untracked, not ignored | move to `docs/` or delete |
| `v2.31.0.tar.gz` (609 KB, vsearch source) | ignored, still at the root | delete |
| 8 `run_*.log` (~1.5 MB) | ignored | move to `logs/`, which the night scripts write to |
| `reviews/critique_2026-09-22.md`, `reviews/literature_2026-09-22.md`, `reviews/proposals_for_dbpq.md`, `reviews/critique_fable/`, this file | review documents at the root | one `reviews/` folder |

Also: `.gitignore:52` says the figures are regenerated by
`analysis/benchmark.qmd`, which does not exist. Of the dead derived databases
(§1.6), about 3.6 GB are on disk.

---

## 6. Documentation drift (C22-S6)

`docs/objectives_design.md` is right: OD 33 (lca cutoffs {0.51, 0.8, 1}) and
OD 34 (ITSx inputs only) are both written there. The four documents an agent
or a newcomer reads first still describe the **pre-OD 34 grid**: three inputs
(raw ASV, `d_vs` OTU, ITSx ASV), 96 computations, 504 rows, lca
{0.8, 0.9, 1}. AGENTS.md is the file every agent session is told to read
first, and it describes a pipeline that no longer exists.

| Where | States | Is (code, 2026-10-05) |
|---|---|---|
| `AGENTS.md:25` | Tedersoo "denoised …, not assigned yet" | assigned (`store_assign_taxo_tedersoo_illumina`, 1 398 targets) |
| `AGENTS.md:41`, `README.md:22, 24, 34-35`, `CONTEXT.md:23, 127, 135`, `docs/experimental_design.md:5-6, 38, 72, 196` | inputs ASV + OTU (`d_vs`) + ITSx ASV; 96 computations, 504 rows per dataset; lca `lca_cutoff` ∈ {0.8, 0.9, 1} | OD 34: inputs `<region>_<unit>` only, 36 computations and 189 rows **per input** (Pauvert 2 inputs: 72 / 378; Tedersoo and in silico 6: 216 / 1 134; wood and air 1: 36 / 189); lca {0.51, 0.8, 1} (`R/values_map.R:20`) |
| `AGENTS.md:43` | chapter 01 scores `d_all_taxo` (raw and ITSx columns) and `d_all_taxo_otu` | `d_all_taxo_<region>_asv` / `_otu` |
| `AGENTS.md:66-68` | `R/` lists `itsx.R` | deleted; missing `assignment_helpers.R`, `bio_loaders.R`, `dada2_patched.R`, `insilico_input.R`, `itsx_inputs.R` |
| `AGENTS.md:99-106` | a 6-row table (7 projects); Tedersoo "first run planned the night of 2026-09-30" | 13 projects in `_targets.yaml` (in silico, air and wood missing from the table); all run |
| `AGENTS.md:112` | `cross_val`: "10 folds" | `cv_fold_number <- 3L` (`config.R:364`) |
| `AGENTS.md:167, 169-170, 173`, `CONTEXT.md:119` | `d_asv_shuffled` → `d_asv_for_assignation`; `d_vs_file` → `d_vs_for_assignation`; `itsx_db_list`; `preprocess`; `d_all_taxo` combined on `d_asv_for_assignation` | none of these targets or columns exists; `with_controls()` on `d_<region>_<unit>_real` (`pipelines/assign_taxo.R:221-286`) |
| `AGENTS.md:195`, `docs/experimental_design.md:263` | "1.2–1.3 GB for sintax, lca and blastn" | the R parent process, not vsearch / BLAST+ (§1.3) |
| `AGENTS.md:200-204`, `README.md:61` | "twelve files"; lists `test_itsx.R` | 15 files, 380 expectations; `test_itsx.R` deleted; 4 new files not listed |
| `AGENTS.md:209`, `docs/diagnostics_cv_and_databases.md:311` | `n_workers = 3` | `n_workers <- 2` and `single_worker <- TRUE` (`config.R:302, 310`), so `n_workers` is ignored |
| `AGENTS.md:128-130` | in silico chapter "exploratory … the targets projects still to create" | projects built and scored |
| `config.R:215` | `bio_dataset` is "NULL for the production projects" | it is a string (`"assign_taxo"`, `"mini"`, …); `bio_cfg` is the one that is `NULL` |
| `config.R:270-274` | "Only ITS1 is kept" | per-dataset regions since OD 34 |
| `config.R:243-244` | "Tedersoo's own truth (103 species) is still to be fetched" | `mock_sanger_tables` row since 0.6 |
| `R/autometric_helpers.R:33` | "invalidates 604 targets" | a different number on every store since OD 34; drop the number |
| `R/mock_truth.R:29-46` | a block of measurements in a code comment, with a broken sentence ("… and 47 without truth. The / The GBIF accepted name …", lines 40-41) | measurements belong in HISTORY or ETAT (below) |
| `analysis/01_load_and_clean.qmd:547-549`, `analysis/07_q4_its_region.qmd:186-191`, `ROADMAP.md:49, 53-54` | "an ASV without the region is not assigned on it" | it is, untrimmed (§1.1) |

**Why it keeps happening, and the fix.** Counts and grid sizes are written as
prose in five documents and some code comments, so every design change has
to be hunted down by hand, and it never is. Your own convention already has
the cure: **`ETAT.md`**, generated by a script, with the dated figures; the
other documents state decisions and point to it. Concretely:

1. `tools/etat.R` computes, per project, `nrow(build_values_map(...))` and
   the number of distinct `compute_name`, the number of targets
   (`tar_manifest()`), outdated and orphan targets, units per input, truthed
   and foreign units (from `res_comp_tax`), and the number of test files and
   expectations. It prepends today's block to `ETAT.md`.
2. In AGENTS, README and CONTEXT, replace every count with "see `ETAT.md`",
   and every grid description with a pointer to `build_values_map()` and OD
   34. Keep the *why*, drop the *how many*.
3. Run `tools/etat.R` at the end of `run_bio_night.R`.

Until that exists, fix the rows of the table above. AGENTS.md first, since it
steers every agent session.

---

## 7. Tests

**Run tonight**, from the root, as AGENTS says: `Rscript tests/run_all.R` ✓.
That is 15 files, 380 expectations, 0 failed, 0 errors, **0 skipped**
(vsearch, BLAST+, cutadapt and ITSx were all found), in 38 s. The suite is
fast and green. Its gaps are where the scientific risk is.

**What the suite does not test:**

1. **The two methods the paper ranks first are not checked against their
   direct call.** `tests/test_assign_compute.R:33, 62` proves that the
   derived **lca** and **blastn** rows equal the direct assignment. Nothing
   does the same for `derive_sintax_row()` or `derive_dada2_row()`, the
   reimplementations of `assign_sintax(min_bootstrap = )` and
   `assignTaxonomy(minBoot = )`. sintax is deterministic with `seed` and one
   thread, so its test is direct:

   ```r
   test_that("sintax rows derived from one search equal a direct call at each min_bootstrap", {
     skip_if_not(dbpq::is_vsearch_installed())
     computed <- compute_assignment(pq, "sintax", ref_fasta = ref, nproc = 1)
     for (b in c(0.4, 0.5, 0.6)) {
       suffix <- paste0("_sintax_", b)
       derived <- derive_assignment(pq, computed, "sintax", suffix, min_bootstrap = b)
       direct <- MiscMetabar::assign_sintax(
         pq, ref_fasta = ref, behavior = "add_to_phyloseq", min_bootstrap = b,
         clean_pq = FALSE, nproc = 1, seed = sintax_seed, suffix = suffix
       )
       expect_equal(new_columns(derived, suffix), new_columns(direct, suffix))
     }
   })
   ```

   dada2 cannot be compared reliably with a second call. Even at a fixed
   seed, its tie-break makes two calls differ (§1.9 b), so such a test would
   be flaky. Pin its **rule** instead, on a hand-made `computed` with a
   non-monotone bootstrap row. Then whoever changes the left-packing of §2 L1
   has to change the test on purpose:

   ```r
   test_that("derive_dada2_row keeps the ranks dada2 keeps", {
     tax <- matrix(c("Fungi", "Asco", "Sordario", "Hypo", "Nectria", "Fusarium", "oxysporum"), 1,
                   dimnames = list("ACGT", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")))
     boot <- matrix(c(100L, 100L, 40L, 40L, 40L, 60L, 30L), 1, dimnames = dimnames(tax))
     pq1 <- phyloseq::phyloseq(
       phyloseq::otu_table(matrix(1, 1, 1, dimnames = list("s", "u1")), taxa_are_rows = FALSE),
       phyloseq::tax_table(matrix("x", 1, 1, dimnames = list("u1", "seed"))),
       phyloseq::refseq(Biostrings::DNAStringSet(c(u1 = "ACGT")))
     )
     out <- as(derive_dada2_row(pq1, list(tax = tax, boot = boot), 0.5, "_d")@tax_table, "matrix")
     # dada2 packs the kept ranks to the left: Genus lands in the Class column.
     expect_equal(unname(out[1, c("Kingdom_d", "Phylum_d", "Class_d")]), c("Fungi", "Asco", "Fusarium"))
   })
   ```

   Both snippets pass against the current code (scratch run, with the setup
   lines of `tests/test_assign_compute.R`: 3 + 1 expectations). Add them to
   that file.

2. **The analysis helpers that define every figure have no test.**
   `enrich_metrics()`, `default_view()`, `filter_default_settings()`,
   `default_names_for()`, `view_rows()`, `score_common_units()` and the three
   `sweep_*` functions are untested. Chapter 01 checks two invariants at
   render time (`check-inputs`, `check-species`), and a third one only prints
   (§4.3). A `tests/test_analysis_helpers.R` built on a 12-row fake
   `res_comp_tax` (2 datasets × 2 regions × 3 units, two consensus readings)
   would pin the filters that hold the default-settings and
   one-reading-of-each-consensus semantics. Their failure mode is silent:
   "every row a consensus row, `filter_default_settings()` returns nothing".
3. **`derive_fake_ref()` passes by luck** (§1.2). Make its test seed-proof.
4. **No test for derivation staleness** (§1.6): there is no staleness rule to
   test yet.
5. Still open from 2026-09-22: `cross_val()` (710 lines) has no test, and
   `tests/test_values_map.R` (362 lines) pins byte strings rather than the
   name grammar (R6). A round-trip test would cover the grammar:
   `full_name_for()` → join on `name_map` → the original columns.

---

## 8. Cross-validation (parked, OD 28)

Nothing in the Results depends on it today. `store_cross_val` predates the
nine-database grid, and the developer parked the CV axis. So, briefly:

- **Unchanged since 2026-09-22 S1**: 710 lines, the metric block written
  twice (vector and scalar `min_bootstrap`), `NA_classif` / `NA_database`
  computed then dropped (`R/cross_val.R:411-412, 560-561`), and a
  floating-point `!=` on a sum of proportions (`:472-477, 600-605`).
- **New, in `compute_by_tax_level`** (`:422-464`, `:571-598`):
  - **FP and FN are swapped**. `FP` counts the wrong answers among the
    queries whose *truth* is the taxon, which is FN for that taxon. `FN`
    counts them among the queries *assigned* the taxon, which is FP. F1 is
    symmetric in FP + FN, so it survives; the labelled TP / FP / FN rows do
    not.
  - **Kingdom is dropped**: `colnames(assign_res$taxo[, -1])` removes a
    first column that is a rank, because `taxa_names` was already removed at
    `:364`.
  - **NA answers are never FN**: `sum(!good, na.rm = TRUE)` drops them.
  - Only `R/examples_cross_val.R`, which nothing sources, turns this branch
    on.
- **The vector branch NA-masks cumulatively** (`:406`):
  `assign_res$taxo[assign_res$taxo_bootstrap < min_bootstrap[i]] <- NA`
  modifies the table it reads at the next `i`. It is correct only when
  `min_bootstrap` is sorted ascending. `c(0.6, 0.4)` gives 0.4 the 0.6
  result.
- **`create_fake_pq_from_refseq()`** parses up to 7 ranks with
  `str_split_fixed(…, n_taxrank)`. A header with more fields puts the rest
  into Species.

**Recommendation.** Do not fix this code. When OD 28 is revisited, rewrite
`cross_val()` around the scoring the mocks already use. A CV query **is** a
unit with a per-unit truth, so `comparpq::tc_metrics_unit()` gives the same
F1 / PPV per rank as the mocks, with controls if they are added to the folds
(2026-09-22 P12). About 400 lines go away, and CV and mock figures become
comparable by construction. Until then, move `R/cross_val.R`,
`R/examples_cross_val.R`, `R/cv_to_tidy.R`, `pipelines/cross_val.R`,
`run_cv_night.R` and chapter 06 to `Archives/cross_val/`, together with the
store pointer, so that nobody reads stale CV figures as current (§4.3).

---

## 9. Status of the code items of 2026-09-22

| Item | Status 2026-10-05 | Evidence |
|---|---|---|
| S1 `cross_val()` size, duplicated blocks, dropped `NA_classif` | **open** | §8 |
| S2 dead `derive_*` functions; sed/awk FASTA surgery | dead functions **done** (deleted 2026-09-22); sed/awk **open** (`general_header_rules`, `derive_no_pattern`, `derive_mini`) | `make_databases.R:227-247, 328-393` |
| S3 metadata parsed from target names | **half**: `enrich_metrics()` joins `name_map`; `assignment_db()` still parses, and "default settings" is still defined twice | §4.3 |
| S4 chapter 01 is four programs | **half**: helpers moved to `R/assignment_helpers.R`; consensus construction, scoring orchestration and three sweeps still inline, and the chapter grew from 610 to 759 lines | `analysis/01_load_and_clean.qmd` |
| S5 stray files, `examples_cross_val.R` | **open, worse** (more logs, a tracked `Rplots.pdf` and log, a 4 000-line sandbox) | §5 |
| S6 design facts in several documents | **open, worse** since OD 34 | §6 |
| S7 ROADMAP as changelog | **done** (2026-09-30 split) | `ROADMAP.md` |
| R1 `rep("sintax", length(db))` trick for `ifelse()` | **open** | `R/values_map.R:213-222` |
| R2 `derive_assignment()` silent on `NULL` | **open** | `R/assign_compute.R:82-84` |
| R3 `db_shapes` sliced, not recycled | **open**, harmless while there are 9 databases | `analysis/00_setup.R:179-182` |
| R4 crew and `ref_file_targets` blocks duplicated across pipelines | **open** | `pipelines/assign_taxo.R:33-62, 166-169`; `pipelines/cross_val.R:36-58, 175-196` |
| R5 `config.R` as a module | **open, grown**: `bio_datasets` has 13 columns, 2 of them `NA` everywhere | §5 |
| R6 test balance | **open** | §7 |

---

## 10. Fix plan, grouped by what it reruns

**The M2.4 final timing run is the natural moment for every [A] fix.** It
changes the threads (`assign_threads_*`) and, through §1.3, the code of
`with_autometric()`. Both enter the command of every compute target, so M2.4
reruns every computation of every store it touches **anyway**. Any [A] fix
made before M2.4 costs nothing extra. Any [A] fix made after M2.4 costs a
second full rerun. Hence this order.

**[0] Now, no rerun:**

1. Commit this project and comparpq (§1.7). Add the `dirty` flag to
   `pqverse_versions()` and the version attribute to `res_comp_tax`.
2. Rebuild and cap `EUK_ITS_v2.1_Fungi_cut` before the zoom (§1.6). Add the
   staleness rule and the atomic writes to `make_databases.R`.
3. `autometric_dir_assign` from `TAR_PROJECT` (§3.1); delete the mini logs
   from the production folder. *Done 2026-10-07 for the folder; the old
   logs stay in `assign_taxo/`, superseded by the newer M2.4 runs.*
4. `run_projects()` with pruning (§3.2, §1.8); port or archive the five
   sandbox scripts that read orphans.
5. Fix the `derive_fake_ref()` code and make its test seed-proof (§1.2); the
   file stays until Q2.
6. Docs: AGENTS.md first, then README, CONTEXT and experimental_design
   (§6). Start `tools/etat.R` → `ETAT.md`.
7. Dead code and stray files (§5); archive the CV code (§8).
8. The sintax equivalence test and the dada2 rule test (§7).

**[R] Re-render chapters, no reassignment:**

9. Q4 on trimmed units only: the `trimmed` flag in chapter 01 and the paired
   reading on the 150 (expected: same order, full > ITS1 > ITS2), the
   one sentence on the in silico "full" input, then re-render 01 and 07, and
   correct ROADMAP and HISTORY (§1.1). *Done 2026-10-07; the 150 and the
   183 are both shown.*
10. Q2.5 one point per computation (§1.5); CO₂ with `memory_ram` (§1.4, once
    §1.3 gives real peaks); ties (§4.1); warnings on in 01 (§4.2); chapter
    06 out of the render list (§4.3); `check-cells` as `stopifnot()`.
11. Measure the seed spread of dada2 and sintax on a few combinations and
    report it beside their comparisons (§1.9, step 1; this costs compute,
    not reassignment of the stores).

**[A] Bundle into the M2.4 rerun:**

12. Process-tree RSS sampler + child CPU + the S1.6a `log_stop()` guard
    (§1.3). *Done 2026-10-07 (`R/sample_tree.R`, `add_tree_costs()`).*
13. ~~The seeded tie-break patch for dada2 (§1.9, step 2).~~ Not wanted
    (Q4 answer, 2026-10-06).
14. `derive_fake_ref(force = TRUE)` with the fixed code, then
    `derive_fungi_rep(force = TRUE)`, which excludes candidates near the
    controls (§1.2; Q2 answer, 2026-10-06: regenerate).
15. `--notrunclabels`, `lca_maxrejects`, the `taxa_names`-free bootstrap
    mask in `derive_sintax_row()`, and the `warning()` in
    `derive_assignment()` (§2, §3.5, §5). *Done 2026-10-07; the new mask
    gives identical tables on the 27 stored Pauvert sintax cases.*
16. Optionally `ITSx_detected` written by `itsx_pq()` (§1.1, step 5), and
    the dada2 worker guard target (§3.4). *Guard target done 2026-10-07
    (`dada2_build_on_worker`); `ITSx_detected` not done.*

---

## 11. Questions for the developer

**Q1. `keep_undetected = TRUE`.** For Q1–Q3 and Q5, keeping the units whose
region ITSx did not find, on their full sequence, is defensible: no unit is
lost. For Q4 it mixes regions. Is the plan of §1.1 right: keep them in the
main grid, and score Q4 on the units trimmed in every region? Or should the
region inputs drop the undetected units (an [A] change, best bundled with
M2.4)?
*Answer (2026-10-06): keep them in the main grid, and score Q4 on the units
trimmed in every region.*

**Q2. The external controls cover 65 of 73 phyla, the 8 rarest missing
(§1.2).** Do we document it, or regenerate the controls with the fixed code
in the M2.4 rerun? Regenerating changes every control-based rate (`NA_fake`
is unaffected; `ext_fungi` and `ext_correct` are affected) and the
`_Fungi_rep` draws.
*Answer (2026-10-06): regenerate them in the M2.4 rerun.*

**Q3. The in silico "full" input.** It is the ASV itself, which is already
the full ITS of its UNITE template (§1.1, corrected 2026-10-06). Keep it in
chapter 07 with one sentence saying so, or drop it anyway?
*Answer (2026-10-06): keep it, with one sentence.*

**Q4. Seed variance of the k-mer methods (§1.9).** A different seed moves
dada2's Genus PPV by about 0.01 on one combination. Shall I script the
seed-spread measurement (3 seeds × 3 databases × dada2 and sintax, Pauvert,
about 30 min)? And do you want the second dada2 patch for bit-identity at a
fixed seed? It is a small effect (1–2 Species calls), and one more patch to
defend in the Methods.
*Answer (2026-10-06): measure the seed spread, without the patch. Also
measure the impact of the dada2 read-filtering settings recommended for
fungal ITS (maximum expected errors = 8, reads truncated at the first quality
score ≤ 8; Rolling et al., 2022). Today `pipelines/dada2.R:95-110` (Pauvert)
and `pipelines/dada2_bio.R:164-183` (Tedersoo, in silico) pass only
`trimLeft = 1, trimRight = 1` to `filter_trim()`, so `filterAndTrim()` runs
on its defaults (`maxEE = Inf`, `truncQ = 2`). Wood and air are published
OTU tables (`pipelines/load_bio.R`), so the setting cannot apply to them.
Scope: see Q9 and Q10.*

**Q5. Is the code meant to be reused** (published as a harness for other
datasets), or only to back this paper? That decides whether the Pauvert-as-
a-dataset refactor (§3.3) and the database `{targets}` project (§1.6, step 3)
are worth a day each.
*Answer (2026-10-06): paper only. Skip both refactors; keep the staleness rule
and the atomic writes of §1.6.*

**Q6. May a future session prune the stores?** The `store_*/meta/meta` files
are tracked (decision 9), so the first prune makes a large diff, and the five
sandbox scripts of §1.8 will then fail loudly.
*Answer (2026-10-06): yes, after the M2.4 rerun, with `tar_prune()` in
`run_bio_night.R`; archive or fix the five sandbox scripts first.*

**Q7. Where do review documents live?** At the root, as now
(`reviews/critique_2026-09-22.md`, `reviews/literature_2026-09-22.md`, this file), or in a
`reviews/` folder (§5)? And do you want `ETAT.md` (§6)?
*Answer (2026-10-06): a `reviews/` folder, and yes to a script-generated
`ETAT.md`.*

**Q8. The zoom database (§1.6).** `fungi + rep cut` needs a UNITE `_Fungi_cut`
(none is built), a derivation that appends the untrimmed representatives of
`_Fungi_rep` to the `_Fungi_cut` records (none exists), and the OD 32 cap on
both. Should the representatives be the **same draw** as in `_Fungi_rep`?
Reusing `sources/<db>_reps.csv` is what makes `fungi + rep` against
`fungi + rep cut` a paired comparison.
*Answer (2026-10-06): the same draw. Since Q2 regenerates the controls in
M2.4, which redraws `_Fungi_rep`, the zoom database must be built after that
and read the new `sources/<db>_reps.csv`.*

**Q9. Datasets for the read-filtering measurement (Q4 answer).** The
Rolling et al. settings change the ASVs, so each dataset needs a new dada2
store, its ITSx inputs, a new per-unit truth and new assignments. Which of
the three fastq datasets (Pauvert, Tedersoo, in silico)?
*Answer (2026-10-06): Pauvert only, on the 2 best databases with every
method (the final-zoom rule).*

**Q10. What the read-filtering measurement decides.** Is it a sensitivity
analysis reported beside the main grid, or a choice to make before M2.4, so
that the main grid reruns with the better settings?
*Answer (2026-10-06): a choice made before M2.4. Consequences for the order
of work:*

- *the seed spread (Q4) comes first or alongside, because it is the noise
  floor: a filtering effect smaller than the seed spread is no reason to
  switch;*
- *both arms use the current databases and controls, so the comparison is
  paired; the regenerated controls (Q2) wait for M2.4;*
- *if M2.4 adopts the stricter filtering, it applies to Tedersoo and
  in silico too, on Pauvert evidence alone; the Methods must say so.*

*Outcome (2026-10-06): measured with `analysis/sandbox/test_read_filtering.R`
(no change on the units common to both arms; the Genus gain on all units comes
from the unit set). Decision: default filtering kept; the Rolling et al.
settings on other datasets are a side question in `docs/parked.md` §3.*

---

## Appendix: the checks behind the findings

Read-only scripts, run from the project root with `Rscript <file>.R` (no store is modified; `check_dada2_rerun_seeded.R` runs one dada2 computation, about 3 min, and `check_q44_trimmed.R` re-scores Tedersoo, about 5 min). Copy them to `analysis/sandbox/` if you want to keep them.

### `check_itsx_undetected.R` (§1.1, units left untrimmed by ITSx per store and region)

```r
suppressMessages({library(phyloseq); library(targets)})
stores <- c(pauvert = "store_assign_taxo", tedersoo = "store_assign_taxo_tedersoo_illumina",
            insilico = "store_assign_taxo_insilico", wood = "store_assign_taxo_korhonen_wood",
            air = "store_assign_taxo_gssp_air")
for (ds in names(stores)) {
  s <- stores[[ds]]
  d_asv <- tar_read(d_asv, store = s)
  orig <- setNames(as.character(d_asv@refseq), taxa_names(d_asv))
  regs <- sub("^d_(.*)_asv_real$", "\\1", grep("^d_.*_asv_real$", tar_meta(store = s)$name, value = TRUE))
  for (r in regs) {
    x <- tar_read_raw(paste0("d_", r, "_asv_real"), store = s)
    sq <- setNames(as.character(x@refseq), taxa_names(x))
    same <- sq == orig[names(sq)]
    cat(sprintf("%-9s %-5s d_asv=%4d units=%4d merged_away=%4d untrimmed(full seq kept)=%4d\n",
      ds, r, ntaxa(d_asv), ntaxa(x), ntaxa(d_asv) - ntaxa(x), sum(same)))
  }
}
```

### `check_insilico_full.R` (§1.1, why ITSx writes no full ITS for the in silico ASVs)

```r
suppressMessages({library(phyloseq); library(targets); library(Biostrings)})
asv <- tar_read(d_asv, store = "store_assign_taxo_insilico")@refseq
tmpl <- readDNAStringSet("data/data_raw/insilico/iss_input_fungi.fasta")
run_itsx <- function(seqs, label) {
  names(seqs) <- paste0("s", seq_along(seqs))
  fa <- tempfile(fileext = ".fasta"); out <- tempfile()
  writeXStringSet(seqs, fa)
  cmd <- paste0("source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxenv && ",
                "ITSx -i ", fa, " -o ", out, " -t all --cpu 4",
                " --preserve F --save_regions all --graphical F --silent T")
  system2("bash", c("-c", shQuote(cmd)), stdout = FALSE, stderr = FALSE)
  pos <- readLines(paste0(out, ".positions.txt"))
  found <- sapply(c("SSU", "ITS1", "5.8S", "ITS2", "LSU"), function(r) {
    sum(grepl(paste0(r, ":"), pos) & !grepl(paste0(r, ": Not found"), pos))
  })
  full <- paste0(out, ".full.fasta")
  n_full <- if (file.size(full) > 0) length(readDNAStringSet(full)) else 0
  cat(label, ":", length(seqs), "sequences; found", paste(names(found), found), "; full.fasta", n_full, "\n")
}
run_itsx(tmpl, "UNITE templates")
run_itsx(asv, "in silico ASVs")
inside <- sapply(as.character(asv), \(s) any(grepl(s, as.character(tmpl), fixed = TRUE)))
cat("ASVs inside a template:", sum(inside), "/", length(asv), "\n")
```

### `check_common_units.R` (§1.1, trimming status of the Q4.4 common units)

```r
suppressMessages({library(phyloseq); library(targets)})
s <- "store_assign_taxo_tedersoo_illumina"
d_asv <- tar_read(d_asv, store = s)
orig <- setNames(as.character(d_asv@refseq), taxa_names(d_asv))
trimmed <- lapply(c("ITS1", "ITS2", "full"), function(r) {
  x <- tar_read_raw(paste0("d_", r, "_asv_real"), store = s)
  sq <- setNames(as.character(x@refseq), taxa_names(x))
  names(sq)[sq != orig[names(sq)]]
})
units <- lapply(c("ITS1", "ITS2", "full"), function(r) taxa_names(tar_read_raw(paste0("d_", r, "_asv_real"), store = s)))
common <- Reduce(intersect, units)
common_trim <- Reduce(intersect, trimmed)
cat("common real units:", length(common), "\n")
cat("of which trimmed in all three regions:", length(intersect(common, common_trim)), "\n")
cat("untrimmed in ITS1 / ITS2 / full among common:",
    sapply(1:3, function(i) length(setdiff(common, trimmed[[i]]))), "\n")
cat("untrimmed in all three among common:", length(setdiff(common, Reduce(union, trimmed))), "\n")
x <- tar_read_raw("d_ITS1_asv_for_assignation", store = s)
cat("ITS1 for_assignation ntaxa:", ntaxa(x), " fake:", sum(grepl("^fake_", taxa_names(x))), " external:", sum(grepl("^external_", taxa_names(x))), "\n")
```

### `check_fake_ref.R` (§1.2, replay of the external-control draw)

```r
suppressMessages(source("make_databases.R"))
pool <- file.path(sources_dir, paste0(fake_ref_source, "_wo_Fungi.fasta"))
kingdoms <- external_control_kingdoms()
cat("kingdoms:", paste(kingdoms, collapse = ", "), "\n")
seqs <- Biostrings::readDNAStringSet(pool)
rk <- stringr::str_match(names(seqs), "k__([^;]+);")[, 2]
seqs <- seqs[rk %in% kingdoms]
phyla <- stringr::str_match(names(seqs), "p__\\s*(.*?)\\s*;c__")[, 2]
tab <- table(phyla, useNA = "ifany")
cat("pool records:", length(seqs), " phyla:", length(tab), " singleton phyla:", sum(tab == 1), " NA phyla:", sum(is.na(phyla)), "\n")
print(tab[tab == 1])
fr <- Biostrings::readDNAStringSet(here(fake_ref_fasta))
fp <- stringr::str_match(names(fr), "p__\\s*(.*?)\\s*;c__")[, 2]
cat("fake_ref records:", length(fr), " distinct phyla in file:", length(unique(fp)), "\n")
cat("pool phyla missing from fake_ref:", paste(setdiff(names(tab), fp), collapse = ", "), "\n")
# replay the draw with the current code
set.seed(targets_seed)
seqs2 <- Biostrings::readDNAStringSet(pool)
rk2 <- stringr::str_match(names(seqs2), "k__([^;]+);")[, 2]
seqs2 <- seqs2[rk2 %in% kingdoms]
ph2 <- stringr::str_match(names(seqs2), "p__\\s*(.*?)\\s*;c__")[, 2]
one <- vapply(unique(ph2), function(p) sample(which(ph2 == p), 1), integer(1))
cat("one_per_phylum picks whose phylum differs from the intended one:",
    sum(ph2[one] != unique(ph2)), "of", length(one), "\n")
cat("replayed draw identical to file names:", identical(sort(names(seqs2)[unique(c(one))]) %in% names(fr), rep(TRUE, length(unique(one)))), "\n")
```

### `check_orphans.R` (§1.8, orphan and outdated targets)

```r
for (p in c("assign_taxo", "assign_taxo_tedersoo_illumina", "assign_taxo_korhonen_wood")) {
  Sys.setenv(TAR_PROJECT = p)
  man <- suppressMessages(targets::tar_manifest(fields = "name", callr_function = NULL))
  meta <- targets::tar_meta(fields = c("name", "bytes", "type"))
  meta <- meta[meta$type != "function" & meta$type != "object" | is.na(meta$type) | TRUE, ]
  orphans <- setdiff(meta$name, man$name)
  orphan_targets <- intersect(orphans, meta$name[!is.na(meta$bytes) & meta$bytes > 0])
  store <- targets::tar_config_get("store")
  obj <- file.path(store, "objects", orphan_targets)
  cat(sprintf("%-32s manifest=%4d meta=%5d orphan_meta_rows=%5d orphan_objects_on_disk=%4d (%.1f MB)\n",
    p, nrow(man), nrow(meta), length(orphans), sum(file.exists(obj)), sum(file.size(obj[file.exists(obj)])) / 1e6))
  out <- targets::tar_outdated(callr_function = NULL, reporter = "silent")
  cat("  outdated now:", length(out), if (length(out)) paste(head(out, 5), collapse = ", ") else "", "\n")
}
```

### `check_mini_logs_exact.R` (§3.1, mini-run logs in the production cost folder)

```r
suppressMessages(library(dplyr))
dir <- "data/data_final/autometric/assign_taxo"
files <- list.files(dir, pattern = "\\.txt$")
logs <- tibble(
  file = files,
  phase = sub("__\\d{8}_\\d{6}\\.txt$", "", files),
  stamp = as.POSIXct(sub("^.*__(\\d{8}_\\d{6})\\.txt$", "\\1", files), format = "%Y%m%d_%H%M%S")
)
mini <- targets::tar_meta(store = "store_assign_taxo_mini", fields = c("name", "time", "seconds"))
prod <- targets::tar_meta(store = "store_assign_taxo", fields = c("name", "time", "seconds"))
# A log belongs to the run whose target of the same name finished within
# [stamp, stamp + seconds + 120 s].
match_run <- function(meta) {
  logs |>
    inner_join(meta, by = c("phase" = "name")) |>
    filter(time >= stamp, time <= stamp + seconds + 120) |>
    distinct(file)
}
m <- match_run(mini)
p <- match_run(prod)
cat("log files:", nrow(logs), "\n")
cat("matching the LAST mini build of their phase:", nrow(m), "\n")
cat("matching the LAST production build of their phase:", nrow(p), "\n")
cat("matching both:", length(intersect(m$file, p$file)), "\n")
newest <- logs |> group_by(phase) |> slice_max(stamp, n = 1, with_ties = FALSE) |> ungroup()
cat("newest-per-phase logs that are mini logs:", sum(newest$file %in% m$file), "\n")
```

### `check_dada2_rerun_seeded.R` (§1.9 (b), dada2 rerun with the target's own seed)

```r
suppressMessages(library(targets))
s <- "store_assign_taxo"
name <- "ITS1_asv_compute_dada2__Unite_all_20250219_Fungi"
db <- "Unite_all_20250219_Fungi"
stored <- tar_read_raw(name, store = s)
input <- tar_read_raw("d_ITS1_asv_for_assignation", store = s)
seed <- targets::tar_seed_create(name, global_seed = 22)
cat("target seed:", seed, "\n")
run <- function() {
  set.seed(seed)
  dada2::assignTaxonomy(
    input@refseq,
    refFasta = file.path("data/data_raw/refseq/dada2_format", paste0(db, ".fasta")),
    minBoot = 0, outputBootstraps = TRUE, multithread = 3
  )
}
fresh <- run()
saveRDS(fresh, file.path(tempdir(), "fresh_dada2.rds"))
cat("same target seed, same threads (3): bootstraps differ on", sum(rowSums(stored$boot != fresh$boot) > 0), "of", nrow(stored$boot), "units\n")
cat("max absolute bootstrap difference:", max(abs(stored$boot - fresh$boot)), "\n")
cat("best-hit Genus differs on", sum(stored$tax[, "Genus"] != fresh$tax[, "Genus"], na.rm = TRUE), "units\n")
```

### `check_q44_trimmed.R` (§1.1, Q4.4 re-scored on the units trimmed in every region)

```r
suppressMessages({
  source(here::here("analysis/00_setup.R"))
  source(here("R/mock_truth.R"))
  source(here("R/score_assignments.R"))
  source(here("R/assignment_helpers.R"))
})
ds <- "tedersoo_illumina"
store <- here("store_assign_taxo_tedersoo_illumina")
dstore <- here("store_dada2_tedersoo_illumina")
cfg <- mock_sanger_tables[mock_sanger_tables$dataset == ds, ]
d_asv <- tar_read(d_asv, store = dstore)
truth <- build_mock_truth(d_asv, cfg)
orig <- setNames(as.character(d_asv@refseq), taxa_names(d_asv))
vm <- values_map_for(ds)
default_names <- function(inp) {
  vm |>
    filter(input == inp) |>
    filter(
      (method %in% c("dada2", "sintax") & min_bootstrap == 0.5) |
        (method == "lca" & lca_cutoff == 0.51) |
        (method == "blastn" & vote_algorithm == "rel_majority" & min_id == 95)
    ) |>
    pull(full_name)
}
regions <- c("ITS1", "ITS2", "full")
preps <- lapply(regions, function(r) {
  d <- tar_read_raw(paste0("d_all_taxo_", r, "_asv"), store = store) |>
    clean_taxo() |>
    add_gen_sp()
  real <- tar_read_raw(paste0("d_", r, "_asv_real"), store = store)
  now <- setNames(as.character(real@refseq), taxa_names(real))
  list(region = r, d = d, trimmed = names(now)[now != orig[names(now)]])
})
real_names <- lapply(preps, function(p) {
  tn <- taxa_names(p$d)
  tn[!grepl(control_pattern, tn)]
})
common <- Reduce(intersect, real_names)
common_trim <- intersect(common, Reduce(intersect, lapply(preps, `[[`, "trimmed")))
cat("common:", length(common), " trimmed in all three:", length(common_trim), "\n")
score_set <- function(p, keep, label) {
  cols <- default_names(paste0(p$region, "_asv"))
  ranks_df <- build_ranks_df(p$d, cols, ranks = tax_order)
  kept <- prune_taxa(taxa_names(p$d) %in% keep | grepl(control_pattern, taxa_names(p$d)), p$d)
  score_assignments(
    kept, ranks_df, truth, db_meta,
    external_truth = external_control_truth(taxa_names(p$d)),
    truth_ranks = truth_ranks
  ) |>
    mutate(region = p$region, set = label)
}
res <- bind_rows(
  lapply(preps, score_set, keep = common, label = "183 common"),
  lapply(preps, score_set, keep = common_trim, label = "150 trimmed")
) |>
  mutate(tax_level = ifelse(tax_level == "Species", "Gen_sp", tax_level)) |>
  left_join(distinct(vm, full_name, method), by = c("method_db" = "full_name"))
saveRDS(res, file.path(tempdir(), "q44_res.rds"))
for (rk in c("Genus", "Gen_sp")) {
  cat("\n==", rk, "F1, mean over the 36 default columns ==\n")
  print(res |> filter(metrics == "F1_score", tax_level == rk) |>
    group_by(set, region) |> summarise(F1 = round(mean(values), 4), .groups = "drop") |>
    tidyr::pivot_wider(names_from = region, values_from = F1))
  cat("-- per method --\n")
  print(res |> filter(metrics == "F1_score", tax_level == rk) |>
    group_by(set, method, region) |> summarise(F1 = round(mean(values), 4), .groups = "drop") |>
    tidyr::pivot_wider(names_from = region, values_from = F1) |> arrange(method, set), n = 20)
}
cat("\n== check against the cached res_paired (183, chapter 01) ==\n")
print(cache_read("res_paired") |> filter_default_settings() |>
  filter(metrics == "F1_score", tax_level == "Genus") |>
  group_by(region) |> summarise(F1 = round(mean(values), 4)))
```
