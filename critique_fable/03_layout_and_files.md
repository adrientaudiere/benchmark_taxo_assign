# 03 — Layout, stray files, git hygiene

## Current root (39 tracked files, 13 of them at the root)

```
.Rprofile (empty)              ROADMAP.md                     script_assign_taxo_parallel.R
.gitignore                     benchmar_assign_taxo.Rproj     script_cross_val.R
41576_2023_679_MOESM2_ESM.ods  config.R                       script_dada2.R
Archives/                      data/                          store_assign_taxo/
CLAUDE.md                      design_analysis.svg            store_cross_val/
CONTEXT.md                     make.R                         store_dada2/
In_silico_simulation.qmd       make_databases.R               tests/
R/                             proposals_for_dbpq.md          v2.31.0.tar.gz
README.md                      analysis/
```

Plus untracked / ignored on disk: `.RData` (31 MB), `.Rhistory`, `renv/`
(539 MB), `.claude/settings.local.json`.

## Stray files [S0.5, S0.6, S3.2]

| File | What it is | Proposal |
|---|---|---|
| `v2.31.0.tar.gz` (600 KB, tracked since commit 17ecc9e) | vsearch 2.31.0 source tarball | `git rm --cached`, add to `.gitignore`, document "vsearch >= 2.31 on PATH" in README. A build artefact does not belong in a data-analysis repo. |
| `41576_2023_679_MOESM2_ESM.ods` (tracked) | Springer Nature supplementary table (MOESM2 naming; journal id 41576 = Nature Reviews Genetics, 2023). Not referenced by any script or doc. | Move to `docs/references/<firstauthor>_2023_supp2.ods` with one line in `docs/references/README.md` saying what it is used for, or delete if it was only consulted once. |
| `design_analysis.svg` (110 KB, tracked) | Design diagram of the analysis, February 2025 | Move to `docs/`, link it from README. |
| `In_silico_simulation.qmd` | The D1b notebook | Move to `analysis/` next to `benchmark.qmd`. Two notebooks in two places with two `here::i_am()` conventions is one too many. |
| `.RData`, `.Rhistory` | RStudio session state | Delete; set `RestoreWorkspace: No`, `SaveWorkspace: No`, `AlwaysSaveHistory: No` in the `.Rproj`. |
| `renv/` | Legacy library | Delete (already ignored). |
| `.claude/settings.local.json` | Per-machine Claude permissions | Add `.claude/settings.local.json` to `.gitignore`. |
| `.Rprofile` (empty) | – | Either delete, or use it for `Sys.setenv(TAR_PROJECT = ...)` defaults once `_targets.yaml` exists. |

## Proposed tree

Changes are marked `+` (new) and `>` (moved/renamed). Everything else stays.

```
benchmar_assign_taxo/
├── README.md  ROADMAP.md  CLAUDE.md  CONTEXT.md
├── config.R
├── make.R
├── make_databases.R
├── + _targets.yaml            # three named projects, see below
├── > pipelines/               # optional; see "Should the scripts move?"
│   ├── > dada2.R
│   ├── > assign_taxo.R
│   └── > cross_val.R
├── R/
│   ├── + load_pqverse.R       # single dependency loader (critique 02)
│   ├── + values_map.R         # build_methods_grid(), build_values_map(), build_cv_values_map()
│   ├── > create_fake_pq_from_refseq.R   # was functions.R
│   ├── combine_taxo_assignments.R
│   ├── cross_val.R
│   └── cv_to_tidy.R
├── analysis/
│   ├── + _quarto.yml
│   ├── > benchmark.qmd        # split into chapters, see critique 05
│   ├── > in_silico_simulation.qmd
│   └── + sandbox/             # "En chantier" material out of the render path
├── + figures/                 # ggsave() output; tracked or not is a ROADMAP question
├── + docs/
│   ├── > design_analysis.svg
│   ├── > examples_cross_val.R # was R/examples_cross_val.R (not sourced by anything)
│   ├── proposals_for_dbpq.md  # (moving it is optional)
│   └── references/
│       └── > <named>.ods
├── data/                      # unchanged, plus + data/README.md
├── store_dada2/  store_assign_taxo/  store_cross_val/   # unchanged
│   (+ store_*_mini/ for smoke runs, via _targets.yaml)
├── tests/
│   ├── + run_all.R
│   └── test_*.R
├── Archives/                  # unchanged
└── critique_fable/            # this review
```

## `_targets.yaml`: name the pipelines once [S3.1]

targets supports several pipelines per repo through `_targets.yaml`. It
removes the `script = ..., store = ...` pair from every `tar_make()`,
`tar_read()`, `tar_visnetwork()`, `tar_poll()` call in `make.R`, CLAUDE.md,
the ROADMAP tables and the notebooks:

```yaml
dada2:
  script: script_dada2.R
  store: store_dada2
assign_taxo:
  script: script_assign_taxo_parallel.R
  store: store_assign_taxo
assign_taxo_mini:
  script: script_assign_taxo_parallel.R
  store: store_assign_taxo_mini
cross_val:
  script: script_cross_val.R
  store: store_cross_val
cross_val_mini:
  script: script_cross_val.R
  store: store_cross_val_mini
```

Usage: `Sys.setenv(TAR_PROJECT = "assign_taxo"); tar_make()` or
`tar_make(project = "assign_taxo")` is not a real argument, so `make.R`
becomes:

```r
run_project <- function(name) {
  withr::with_envvar(c(TAR_PROJECT = name), {
    targets::tar_make()
    targets::tar_prune()
  })
}
run_project("dada2")
run_project("assign_taxo")
run_project("cross_val")
```

The `*_mini` projects solve B2: the script reads
`mini_db <- grepl("_mini$", Sys.getenv("TAR_PROJECT"))` and the smoke-test
run can never overwrite the production store.

## Should the scripts move into `pipelines/`?

Pros: the root goes from 13 files to 8; the three scripts stop being
confused with `make.R` and `make_databases.R`. Cons: three `here::i_am()`
lines and the `_targets.yaml` paths change; CLAUDE.md and the ROADMAP
command tables need a search-and-replace. With `_targets.yaml` in place the
scripts are never typed by hand again, so the move is cheap. Recommended
but not required; do it together with S3.1 or not at all.

## Git hygiene

- `store_*/meta/meta` are tracked on purpose (the per-store `.gitignore`
  keeps `meta/` and drops `objects/`). That is a reasonable choice for
  provenance, but it means every `tar_make()` dirties the tree. If that
  becomes annoying, an alternative is to keep `tar_meta()` exported as CSV
  (`data/data_final/tar_meta_<store>.csv`) by `make.R` and ignore the raw
  meta files. Decision for the developer (see ROADMAP open questions).
- The folder is `benchmar_assign_taxo` (missing `k`), the `.Rproj` follows
  the folder, and README's title is `benchmark_taxo_assign`. Renaming the
  folder breaks Nextcloud sync history and the Claude memory path, so keep
  the folder and align the README title and the `.Rproj` display name
  instead. `[Low / easy]`
- Commit 17ecc9e added the tarball; if history size matters, a
  `git filter-repo` pass would drop it, but at 600 KB it is not worth the
  rewrite. Just stop tracking it going forward.
