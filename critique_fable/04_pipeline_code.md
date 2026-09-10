# 04 — Pipeline code: duplication, naming, tests

## `values_map` exists three times [S2.1]

| Copy | Lines | Differences |
|---|---|---|
| `script_assign_taxo_parallel.R` | 305–347 | adds `controller` column |
| `R/values_map_for_qmd.R` | 14–55 | identical minus `controller`; asserts `here::i_am("analysis/benchmark.qmd")` |
| `script_cross_val.R` (`cv_values_map`) | 458–497 | same DB vector, same `db_path` formula, different method grid |

The six-element DB vector, including the four commented-out SSU variants, is
copy-pasted in all three. Q2.3 in the ROADMAP already records one incident
where the two first copies had to be edited in lockstep. The comparpq
`preference_pattern` and the `db_meta` tibble in `benchmark.qmd` (line 765)
are two more places that list the same databases by hand.

Proposal, `R/values_map.R`:

```r
build_methods_grid <- function(min_bootstrap = c(0.4, 0.5, 0.6),
                               vote_algorithm = c("rel_majority", "abs_majority", "unanimity"),
                               nb_voting = 100) { ... }          # the expand_grid + full_join

db_path_for <- function(method, db_name) { ... }                 # dada2_format vs sintax_format

build_values_map <- function(dbs = db_list, mini_db = FALSE, methods = build_methods_grid()) { ... }

build_cv_values_map <- function(dbs = db_list, mini_db = FALSE, remove_tested = c(TRUE, FALSE)) { ... }

full_name_for <- function(method, db, min_bootstrap, vote_algorithm, nb_voting) { ... }
```

with `db_list` (and `db_meta`, the simplification-type table used by Q2.1)
moved into `config.R`. Then:

- the two pipeline scripts and the qmd `source(here("R/values_map.R"))`;
- `R/values_map_for_qmd.R` is deleted;
- `tests/test_values_map.R` pins the exact set of `full_name` strings the
  store is expected to contain, which is the cheapest guard against the
  "silently renamed targets" failure mode.

`full_name_for()` also fixes the `gsub("...NA...NA", "", ...)` trick: build
the suffix with `paste` only for non-NA parts instead of pasting NA and
deleting it afterwards.

## Naming

- `R/functions.R` contains one function, `create_fake_pq_from_refseq()`.
  Rename the file after the function. `[S2.3]`
- `Phyla` vs `Phylum`: see B6 in critique 01. `[S2.3]`
- `script_assign_taxo_parallel.R`: the `_parallel` suffix only made sense
  while the sequential script was alive next to it; that one is in
  `Archives/`. With `_targets.yaml` the file name stops mattering, so rename
  to `script_assign_taxo.R` (or `pipelines/assign_taxo.R`) at the same time.
- `d_asv_for_assignation_fake` is the *intermediate* (shuffled only) and
  `d_asv_for_assignation` the final (shuffled + external). Reversed reading
  on first sight; `d_asv_shuffled` and `d_asv_with_controls` say what they
  hold.
- Target `benchmark_costs` (assign store) and `benchmark_costs_dada2` (dada2
  store): fine, but both will change with B1; name the new ones after what
  they measure (`assignment_costs`, `denoising_costs`).

## Reproducibility details

- `targets_seed <- 22` is used by `tar_option_set()`; `derive_fake_ref()` in
  `make_databases.R` has its own literal `seed = 22`. Pass `targets_seed`.
- `cross_val()` writes every fold to the fixed path
  `tempdir()/test_refseq.fasta`. Safe under crew (one `tempdir()` per worker
  process) but not if someone runs two `cross_val()` calls with `mirai` or
  `future` inside one process. `tempfile(fileext = ".fasta")` costs nothing.
- `cross_val()` still has the French TODO block and a dead `dada2_2steps`
  branch that `stop()`s. Either delete the branch or move the TODO into the
  ROADMAP (done in this review: see "Open scientific work").
- `n_workers = 3` × `n_threads = 4` plus the dedicated `dada2_ctrl` worker
  = up to 16 threads. Fine on the current machine; make the arithmetic
  explicit in the config comment so the next machine is sized correctly.

## Tests [S4.1]

Today: one file, four expectations, for `combine_taxo_assignments()`. Not
covered although pure and fast:

| Function | Fixture needed | Why it matters |
|---|---|---|
| `cv_to_tidy()` | a fake `cross_val()` result list | the only bridge between CV output and `cv_results` |
| `create_fake_pq_from_refseq()` | 3-record sintax fasta string | used by every CV fold |
| `make_databases.R::derive_no_pattern()` / `derive_mini()` / `derive_fake_ref()` | 6-record two-line and multi-line fasta | shell one-liners are where silent truncation happens (proposals_for_dbpq.md already documents one such bug) |
| `build_values_map()` (after S2.1) | none | pins target names |
| `benchmark_costs` aggregation (after S1.1) | two small autometric log files | B1 would have been caught by a test that feeds a log with two runs of the same phase |

Add `tests/run_all.R`:

```r
here::i_am("tests/run_all.R")
testthat::test_dir(here::here("tests"), reporter = "summary")
```

and mention `Rscript tests/run_all.R` in README instead of the single file.

## `make.R`

- Its header says "run this file to (re)build everything" while CLAUDE.md
  and CONTEXT.md say it runs only the first two pipelines and that CV "is
  run independently". The file runs all four steps. Pick one and align the
  docs (critique 05 lists the sentences).
- Add `tar_prune()` after each `tar_make()` (B8) and a final
  `tar_meta()` export if the meta files stop being tracked.
