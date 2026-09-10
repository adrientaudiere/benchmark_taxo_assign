# benchmark_taxo_assign

Benchmark of taxonomic-assignment methods × reference databases for ITS fungal
metabarcoding, supporting a manuscript-in-preparation in the
[pqverse](https://github.com/adrientaudiere/pqverse) ecosystem.

Four assignment methods (`dada2`, `sintax`, `lca`, `blastn`) are run against
nine reference databases (UNITE, Eukaryome ITS/SSU, with `_Fungi`-filtered,
cutadapt-trimmed, and 99%-clustered variants) using a `{targets}` pipeline.
Performance is evaluated against a mock community (Pauvert et al. 2019),
cross-validation, in silico simulations (InSilicoSeq + miaSim), and tree
endophyte data ([Taudière et al. 2018](https://doi.org/10.1016/j.funeco.2018.07.008)).
Five consensus-voting strategies and their parameters are compared via
`comparpq::resolve_taxo_conflict()`.

## Quick start

```r
# Build the DBs (only required once per machine, idempotent)
source("make_databases.R")
derive_all_variants()

# Run the three pipelines in order (named projects in _targets.yaml)
source("make.R")

# Or one pipeline at a time
Sys.setenv(TAR_PROJECT = "assign_taxo"); targets::tar_make()

# Smoke test on the mini_* databases (own stores, never overwrite production)
Sys.setenv(TAR_PROJECT = "assign_taxo_mini"); targets::tar_make()

# Tests
Rscript tests/test_combine_taxo_assignments.R
Rscript tests/test_values_map.R
```

The analysis is a Quarto project in `analysis/` (`quarto render analysis`, or
one chapter at a time starting with `01_load_and_clean.qmd`, which feeds the
others through `analysis/_cache/`). Figures are written to `figures/`
(not tracked).

## Documentation

- **`CLAUDE.md`** — architecture, pipeline layout, external dependencies.
- **`ROADMAP.md`** — running todo list: structural debt (section 1), open
  manuscript work (section 2), scope decisions of 2026-05-19, and the done
  log.
- **`critique_fable/`** — evidence files behind the structural items of the
  ROADMAP (2026-09-10 review).
- **`proposals_for_dbpq.md`** — six helpers still implemented locally that
  could be upstreamed to [`dbpq`](https://github.com/adrientaudiere/dbpq).
- **`docs/`** — design diagram and external reference documents.

## Layout

```
_targets.yaml                 # named targets projects (dada2, assign_taxo, cross_val, *_mini)
pipelines/dada2.R             # ASV pipeline -> store_dada2
pipelines/assign_taxo.R       # 4 methods x 6 DBs in parallel -> store_assign_taxo
pipelines/cross_val.R         # k-fold CV pipeline -> store_cross_val
make.R                        # DB derivation + the three production projects in order
make_databases.R              # idempotent DB derivation (delegates to dbpq::)
config.R                      # shared constants (primers, threads, db_list, seed, CV params)
R/                            # load_pqverse, values_map, autometric helpers, combine, cross_val
analysis/                     # Quarto project: 00_setup.R + chapters 01..06, sandbox/
figures/                      # save_fig() output of the chapters (not tracked)
tests/                        # Rscript-runnable testthat fixtures
docs/                         # design diagram, external references
critique_fable/               # structural review (evidence for ROADMAP section 1)
Archives/                     # superseded scripts (kept for reference)
```

## Dependencies

- R ≥ 4.1
- [`MiscMetabar`](https://github.com/adrientaudiere/MiscMetabar),
  [`dbpq`](https://github.com/adrientaudiere/dbpq),
  [`comparpq`](https://github.com/adrientaudiere/comparpq) — pqverse packages
- [`targets`](https://docs.ropensci.org/targets/),
  [`tarchetypes`](https://docs.ropensci.org/tarchetypes/),
  [`crew`](https://wlandau.github.io/crew/),
  [`autometric`](https://wlandau.github.io/autometric/)
- `cutadapt` installed in a conda env named `cutadaptenv`
- `vsearch` (>= 2.31) and BLAST+ (`blastn`, `makeblastdb`) on `PATH`; install
  them from your package manager or from source, they are not shipped in this
  repository
