# 01 — Bugs and silent traps

Ordered by impact on the manuscript. ROADMAP ticket ids in brackets.

## B1 — `benchmark_costs` aggregates a stale log; the crew workers never log a phase [S1.1] — Confirmed

**Symptom.** The M2 cost layer (`benchmark_costs`, and therefore figures
M2.1 / M2.2 / Q2.5) is built from measurements of the *sequential* pipeline
run of 2025-02-11, not from the parallel run of 2026-05-21 that produced the
current `d_all_taxo`.

**Evidence** (from `data/data_final/autometric_log_assign_taxo.txt`):

| Rows with a 2026 timestamp | Phase value |
|---|---|
| 20 825 | `__DEFAULT__` |
| 6 | `prepare: <target>` / `conclude: <target>` (targets' own hooks, main process) |
| 0 | any `full_name` such as `dada2__Unite___0.5` |

Every row whose phase equals a `full_name` is dated 2025-02-11 and carries
the old database names (`EUK_ITS_v1_9_3`, `mini_Unite`). Yet
`store_assign_taxo/meta/meta` shows `dada2__Unite___0.5` and `d_all_taxo`
were (re)built on day 20594 = 2026-05-21. The stored `benchmark_costs`
object is 1 229 bytes, consistent with an almost-empty inner join.

**Cause.** `log_start()` is called under `if (tar_active())` in the main
process only. Each assignment target runs in a `crew_controller_local`
worker, a separate R process where no logger is running, so
`autometric::log_phase_set(full_name)` has no effect. targets itself sets
`prepare:`/`conclude:` phases in the main process, which is the only phase
traffic the 2026 run produced.

**Second-order problem.** The log file is append-only across runs (rows from
2025-02-07 to 2026-05-21 coexist). `wall_time_s = max(time) - min(time)`
per phase therefore spans every run that ever used that phase name.

**Fix proposal.**

1. Start the logger *inside* the target, in the worker, with one file per
   target and per run:

   ```r
   with_autometric <- function(full_name, expr, dir = here("data/data_final/autometric")) {
     dir.create(dir, showWarnings = FALSE, recursive = TRUE)
     path <- file.path(dir, paste0(full_name, "__", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))
     autometric::log_start(path = path, seconds = 1)
     on.exit(autometric::log_stop(), add = TRUE)
     autometric::log_phase_set(full_name)
     force(expr)
   }
   ```

   and in the `tar_eval` body: `with_autometric(full_name, add_new_taxonomy_pq(...))`.
2. `benchmark_costs` reads the directory (`list.files(..., full.names = TRUE)`
   then `lapply(autometric::log_read)`), keeps the newest file per
   `full_name`, and aggregates. Make the directory a `format = "file"`
   target so a rerun invalidates the aggregation.
3. Alternative to evaluate: crew's built-in autometric integration
   (`crew_controller_local(options_metrics = crew_options_metrics(path = ..., seconds_interval = 1))`).
   It logs worker resource use without touching the targets, but whether
   targets tags the phase with the target name inside the worker must be
   verified before relying on it.
4. Apply the same treatment to `script_dada2.R` (`benchmark_costs_dada2`),
   which shares the design but runs single-process, so it is currently
   correct except for the append-across-runs problem.
5. After the fix, rerun and regenerate the M2 figures. Until then treat
   M2.1, M2.2, M2.3 and Q2.5 as **not done** despite their ticks in the
   ROADMAP.

## B2 — `mini_db <- TRUE` is the committed default, and stores cannot tell mini from full [S0.1] — Confirmed

`config.R` line 66 sets `mini_db <- TRUE`. `values_map` and `cv_values_map`
switch `db_path` to `mini_<db>.fasta` when it is TRUE, but `full_name` (the
target name) is built from `db`, not `db_name`. Consequences:

- `store_assign_taxo` and `store_cross_val` as they stand were built on
  2026-05-21 with the flag on; nothing in the store records that.
- Flipping the flag does not invalidate anything by target *name*; targets
  will still rerun because `db_path` changed in the command, but the two
  runs overwrite the same objects and only one survives.

Fix: (a) read the flag from an environment variable defaulting to FALSE
(`mini_db <- as.logical(Sys.getenv("BENCHMARK_MINI_DB", "FALSE"))`); (b)
route mini runs to their own stores (`store_assign_taxo_mini`,
`store_cross_val_mini`) via `_targets.yaml` projects (see
[03_layout_and_files.md](03_layout_and_files.md)); (c) print one
`message()` line at the top of each pipeline script stating the flag value
and the store path.

## B3 — `figures/` does not exist; every `ggsave()` in `benchmark.qmd` errors [S0.2] — Confirmed

`benchmark.qmd` has 26 `ggsave(here("figures/..."))` calls (Q1.2, Q1.3,
Q1.4, Q2.1, Q2.2, Q2.5, Q3.2, Q3.3, Q3.5, M2.1, M2.2). No `figures/`
directory exists and none is created in the notebook. With the installed
ggplot2 4.0.3 the default is `create.dir = FALSE`, so the first call aborts
the render. Fix: `dir.create(here("figures"), showWarnings = FALSE)` in the
setup chunk (or `create.dir = TRUE`), plus a decision on whether figures are
tracked in git (see open questions in the ROADMAP).

## B4 — `script_dada2.R` points at files and directories that no longer exist [S0.3] — Confirmed

1. `file_refseq_taxo <- here("data/data_raw/refseq/", refseq_file_name)`
   resolves to `data/data_raw/refseq/Unite_Fungi.fasta`, which is not on
   disk; the file lives in `data/data_raw/refseq/dada2_format/`. The store
   dates from 2024-12-19, so any rerun (D1b-A InSilicoSeq, D1c endophytes,
   both planned as copies of this script) fails at `assignTaxonomy`.
2. `lapply(list.files("~/Nextcloud/IdEst/Projets/MiscMetabar/R/", ...), source)`
   targets a directory that no longer exists (the package moved to
   `pqverse/pqverse_pkg/MiscMetabar/`). `list.files()` on a missing directory
   returns `character(0)`, so the line silently does nothing and the script
   runs on whatever `library("MiscMetabar")` finds in the user library, while
   the other two pipelines `load_all()` the development checkout. The two
   stores can therefore be produced by two different MiscMetabar versions
   without any warning.
3. Suspected: `filter_trim(output_fw = paste(getwd(), here("/data/data_intermediate/filterAndTrim_fwd"), sep = ""))`
   concatenates two absolute paths (`here()` already returns an absolute
   path). The output directories exist at the *expected* location, so either
   `filter_trim` normalises the path or the store predates this line. Replace
   with `here("data/data_intermediate/filterAndTrim_fwd")` and rerun once to
   check.

## B5 — Contradictory comments in `config.R` and CLAUDE.md about CV knobs [S0.4] — Confirmed

`config.R`: `cv_fold_tested <- cv_fold_number    # smoke-test default; set to cv_fold_number for the publication run`
— the value *is* already the publication value. `cv_max_seq <- NULL` while
CLAUDE.md says the default is 100. The header of `config.R` still says
`script_dada2.R` has inline copies; it sources `config.R` since C3. Anyone
following the docs to "smoke-test first" gets the opposite of what they
expect.

## B6 — Rank named `Phyla` in the seed taxonomy, `Phylum` everywhere else [S2.3] — Confirmed

`script_dada2.R` passes `taxLevels = c("Kingdom", "Phyla", ...)` to
`assignTaxonomy`. `benchmark.qmd` then needs the regex
`"Kingdom_|Phyla_|Phylum_|..."` to strip rank prefixes, and every downstream
join on rank names has to remember the exception. Rename to `Phylum` at the
source (invalidates `tax_tab` and `d_asv`, which is fine given B4 forces a
rerun anyway).

## B7 — `d_asv` is read across stores without a file dependency [S1.2] — Confirmed

`tar_target(d_asv, tar_read(d_asv, store = here::here("store_dada2")))` in
`script_assign_taxo_parallel.R`: targets sees a constant command, so a
rebuilt `store_dada2` never invalidates `store_assign_taxo`. Fix:

```r
tar_target(d_asv_file, here("store_dada2/objects/d_asv"), format = "file"),
tar_target(d_asv, readRDS(d_asv_file))
```

(`d_asv` is stored as rds; check `tar_meta(d_asv, store = "store_dada2")$format` once.)

## B8 — Stores carry stale targets from three pipeline generations [S0.7] — Confirmed

`store_assign_taxo/meta/meta` has 509 rows, including `*_all_taxo` targets
from the sequential script, `dada2__mini_Unite___0.4NANA` from an early
`full_name` formula, `idtaxa__*` targets, and every `v1_9_3` name.
`store_cross_val` similarly mixes `v1_9_3` and `v2` names. `tar_prune()` on
each store removes the orphans; a `tar_prune()` step in `make.R` after each
`tar_make()` keeps it that way. Because `store_*/meta/meta` is tracked in
git, the pruning also shrinks the diff noise that currently hits three of
the five commits in the history.

## B9 — Session-state files at the root [S0.6] — Confirmed

- `.RData` (31 MB) and `.Rhistory` at the root, with
  `RestoreWorkspace: Default` / `SaveWorkspace: Default` in the `.Rproj`.
  Every RStudio open restores 31 MB of unknown objects into the global
  environment, which is exactly the situation `targets` exists to avoid.
  Set both to `No` and delete the two files.
- `renv/` (539 MB) is still on disk although `.Rprofile` is empty and
  `.gitignore` says renv is no longer used. Delete it.

## B10 — Second mock dataset present but undocumented [S3.2] — Confirmed

`data/data_raw/metadata/taxo_mock_SRR30413326.csv`, the empty directory
`data/data_raw/rawseq/SRR30413326/`, and `data/data_raw/mock_hleap2021/Mocks/`
(five `*_realized.fa/.txt` pairs) exist on disk. No script, notebook, or doc
mentions them. Either they are a planned dataset (then add a ROADMAP item)
or leftovers (then delete). A `data/README.md` listing each subfolder with
its provenance would prevent the next reader from having to guess.

## B11 — All 18 blastn targets of the current store added no taxonomy column, silently [S1.3] — Confirmed (found 2026-09-10 while rendering the split notebook)

`tar_read("blastn__Unite___0.5...rel_majority...100", store = "store_assign_taxo")`
has the same 7 `tax_table` columns as `d_asv_for_assignation`; the target ran
4.6 s and targets recorded neither warning nor error. Same for the other 17
blastn rows, while every dada2 / sintax / lca target added its 7 columns.

**Cause.** `MiscMetabar::assign_blastn()` (`R/blast.R`, lines ~997–1001 and
~1074–1078) emits `message("None blast query match the score filters")` and
`return(physeq)` unchanged when no hit passes the filters. Under targets a
message is not a failure, so the pipeline "succeeds" and the missing columns
only surface in the analysis (`ranks_df` check of `01_load_and_clean.qmd`).
The current store was built with the `mini_*` databases (10 000 records),
which is the likely reason for zero hits; the 2025-02 production run did
produce blastn columns.

**Fix applied.** `combine_taxo_assignments()` now stops and names every
assignment that added zero columns (`allow_empty = TRUE` to opt out); tested
in `tests/test_combine_taxo_assignments.R`. So `d_all_taxo` can no longer be
built from a silently empty assignment. Still to do: check on the production
rerun that blastn does hit with the full databases, and consider a loud
warning upstream in `assign_blastn()` (MiscMetabar).
