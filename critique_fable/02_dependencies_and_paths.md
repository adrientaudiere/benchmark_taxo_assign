# 02 — Dependencies and paths

## The problem: four ways to load the same three packages

| File | MiscMetabar | comparpq | dbpq | greenAlgoR |
|---|---|---|---|---|
| `script_dada2.R` | `library()` + `source()` of a directory that no longer exists | – | – | – |
| `script_assign_taxo_parallel.R` | `library()` **and** `devtools::load_all(<abs path>)` | `source()` of 2 files by absolute path | – | – |
| `script_cross_val.R` | `library()` **and** `devtools::load_all(<abs path>)` | `source()` of 2 files by absolute path | – | – |
| `analysis/benchmark.qmd` | `library()` + `load_all(<abs path>)` | `source()` of 3 files by absolute path | `load_all(<abs path>)` | `load_all(<abs path>)` |
| `In_silico_simulation.qmd` | `library()` + `load_all(~ path)` | `source()` of 2 files by `~` path | – | – |
| `make_databases.R` | – | – | `library("dbpq")` (installed) | – |

Fifteen absolute paths of the form
`/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/...` are spread
over four files (3 + 3 + 6 + 3), plus the dead
`~/Nextcloud/IdEst/Projets/MiscMetabar/R/` path in `script_dada2.R`. Sourcing comparpq file by file (`compare_taxo.R`,
`fake_creation.R`, `taxtab_modification.R`, `compare_taxo_plot.R`) couples
the benchmark to comparpq's internal file layout: renaming a file in
comparpq breaks the benchmark, and any comparpq helper defined in a file
that is not sourced is silently missing. Calling `library()` and then
`load_all()` on the same package attaches the installed version and then
masks it, so which code runs depends on load order.

## Proposal: one loader, one environment variable [S2.2]

`R/load_pqverse.R`:

```r
# Loads the sister pqverse packages from their development checkouts.
# Set PQVERSE_PKG_DIR in .Renviron to override the default location.
pqverse_pkg_dir <- function() {
  Sys.getenv("PQVERSE_PKG_DIR", "~/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg")
}

load_pqverse <- function(pkgs = c("MiscMetabar", "comparpq", "dbpq"),
                         quiet = TRUE) {
  root <- normalizePath(pqverse_pkg_dir(), mustWork = TRUE)
  for (pkg in pkgs) {
    pkgload::load_all(file.path(root, pkg), quiet = quiet)
  }
  invisible(root)
}
```

Then every entry point does exactly two lines:

```r
source(here::here("R/load_pqverse.R"))
load_pqverse()                          # add "greenAlgoR" in the qmd
```

and `library("MiscMetabar")` disappears from the scripts. `pkgload` is what
`devtools::load_all()` calls; using it directly avoids attaching devtools in
every crew worker.

Points to decide:

- Whether the *installed* MiscMetabar (CRAN release) or the dev checkout is
  the reference for the manuscript. Today it is the checkout for two
  pipelines and the installed one for `script_dada2.R` (see B4). Whichever
  is chosen, `sessionInfo()` or `pkgload::pkg_version()` should be captured
  into a target (`tar_target(session_info, sessioninfo::session_info())`)
  so the manuscript can quote versions.
- `here::here()` is used for project files; the pqverse root is the only
  path outside the project. Keeping it in one function is enough; no need
  for a config entry per package.

## Other path issues

- `R/values_map_for_qmd.R` calls `here::i_am("analysis/benchmark.qmd")`.
  A helper in `R/` should not assert which file is calling it; this breaks
  the moment the helper is sourced from a pipeline or a test. (Resolved by
  S2.1, which replaces the file.)
- `benchmark.qmd` calls `here::i_am()` twice and then `setwd(here::here())`
  inside a chunk. Quarto supports `execute: dir: project` in `_quarto.yml`
  (see [05_analysis_and_docs.md](05_analysis_and_docs.md)); with it, both
  the `setwd()` and the second `i_am()` go away.
- `script_dada2.R` and `script_assign_taxo_parallel.R` pass a *relative*
  path to `log_start(path = "data/data_final/...")` but `here()` everywhere
  else. Harmless while `tar_make()` runs from the root; make it `here()` for
  consistency.
- `make_databases.R::derive_fake_ref()` defaults its `input` and `output`
  to relative paths while the orchestrator passes `here()` paths. Same fix.
- `cutadapt_conda_prelude` hard-codes `~/miniforge3/...`. Fine as a
  default, but read it from an env var too (`BENCHMARK_CONDA_PRELUDE`) so a
  second machine only needs a `.Renviron`, not a code edit.

## What a second machine needs today (for the README)

Nothing in the repo lists it end to end. Collected while reading:

1. R >= 4.1 with targets, tarchetypes, crew, autometric, here, conflicted,
   tidyverse core, Biostrings, dada2, vegan, lme4, patchwork, ComplexUpset,
   miaSim, testthat, pkgload.
2. Checkouts of MiscMetabar, comparpq, dbpq, greenAlgoR under one directory
   (`PQVERSE_PKG_DIR`).
3. `vsearch` on PATH (the tracked `v2.31.0.tar.gz` is the vsearch source
   tarball; see [03_layout_and_files.md](03_layout_and_files.md)).
4. BLAST+ on PATH (for `assign_blastn`), not mentioned anywhere.
5. conda env `cutadaptenv` with cutadapt; InSilicoSeq for D1b-A.
6. ~13 GB of reference fasta under `data/data_raw/refseq/` (three source
   files, the rest derived by `make_databases.R`), raw fastq under
   `data/data_raw/rawseq/`.
