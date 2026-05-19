# benchmark_taxo_assign

Benchmark of taxonomic-assignment methods × reference databases for ITS fungal
metabarcoding, supporting a manuscript-in-preparation in the
[pqverse](https://github.com/adrientaudiere/pqverse) ecosystem.

Four assignment methods (`dada2`, `sintax`, `lca`, `blastn`) are run against
seven reference databases (UNITE, Eukaryome ITS/SSU, with `_Fungi`-filtered,
cutadapt-trimmed, and 99%-clustered variants) using a `{targets}` pipeline.
Performance is evaluated against a mock community (Pauvert et al. 2019),
cross-validation, in silico simulations (InSilicoSeq + miaSim), and tree
endophyte data ([Cregger et al. 2018](https://doi.org/10.1016/j.funeco.2018.07.008)).
Five consensus-voting strategies and their parameters are compared via
`comparpq::resolve_taxo_conflict()`.

## Quick start

```r
# Build the DBs (only required once per machine, idempotent)
source("make_databases.R")
derive_all_variants()

# Run both pipelines (parallel; uses crew_controller_local from config.R)
source("make.R")

# Test the combine helper
Rscript tests/test_combine_taxo_assignments.R
```

The headline analysis lives in `analysis/benchmark.qmd`.

## Documentation

- **`CLAUDE.md`** — architecture, pipeline layout, external dependencies.
- **`benchmark_plan.md`** — running todo list mapping each manuscript question
  to concrete tasks, with scope decisions taken on 2026-05-19.
- **`proposals_for_dbpq.md`** — six helpers still implemented locally that
  could be upstreamed to [`dbpq`](https://github.com/adrientaudiere/dbpq).

## Layout

```
script_dada2.R              # ASV pipeline -> store_dada2
script_assign_taxo_parallel.R # 4 methods x 7 DBs in parallel -> store_assign_taxo
make_databases.R            # idempotent DB derivation (delegates to dbpq::)
config.R                    # shared constants (primers, threads, paths, seed)
R/                          # combine helper, cross-validation, helpers
analysis/                   # benchmark.qmd
tests/                      # Rscript-runnable testthat fixtures
Archives/                   # superseded scripts (kept for reference)
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
- `vsearch` on `PATH`
