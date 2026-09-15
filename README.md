# benchmar_assign_taxo

Benchmark of taxonomic-assignment methods × reference databases for ITS fungal
metabarcoding, supporting a manuscript-in-preparation in the
[pqverse](https://github.com/adrientaudiere/pqverse) ecosystem.

Four assignment methods (`dada2`, `sintax`, `lca`, `blastn`) are run against
seven reference databases (UNITE 19.02.2025 with and without singletons,
EUKARYOME ITS v2.1, with `_Fungi`-filtered variants and one variant trimmed to
the ITS1F–ITS2 amplicon; see `docs/reference_databases.md`) using a
`{targets}` pipeline. A post-clustering question (97 % OTUs whose taxonomy
comes from their member ASVs) is derived from the same assignments.
Performance is evaluated against a mock community (Pauvert et al. 2019),
cross-validation, in silico simulations (InSilicoSeq + miaSim), and tree
endophyte data ([Taudière et al. 2018](https://doi.org/10.1016/j.funeco.2018.07.008)).
Five consensus-voting strategies and their parameters are compared via
`comparpq::resolve_taxo_conflict()`.

## Quick start

```r
# Reference databases (idempotent; see docs/reference_databases.md)
source("make_databases.R")
download_reference_sources()   # general FASTA releases listed in config.R
derive_all_variants()          # dada2 / sintax formats, Fungi and cut variants

# Everything in order: databases, then the three production projects (hours)
source("make.R")

# One project at a time (sourcing R/run_project.R runs nothing by itself)
source("R/run_project.R")
run_project("assign_taxo")        # tar_make() + tar_prune() on that project

# Smoke test on the mini_* databases (own stores, never overwrite production)
run_project("assign_taxo_mini")
run_project("cross_val_mini")     # lower the CV knobs of config.R first

# Tests (six files, or one at a time with Rscript tests/test_<name>.R)
Rscript tests/run_all.R
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
- **`docs/reference_databases.md`** — where each reference database comes
  from (UNITE and EUKARYOME releases, URLs, DOIs), how it is converted, and
  what to change when a new release comes out.
- **`docs/`** — also the design diagram and external reference documents.

## Layout

```
_targets.yaml                 # named targets projects (dada2, assign_taxo, cross_val, *_mini)
pipelines/dada2.R             # ASV pipeline -> store_dada2
pipelines/assign_taxo.R       # 4 methods x 8 DBs in parallel -> store_assign_taxo
pipelines/cross_val.R         # k-fold CV pipeline -> store_cross_val
make.R                        # sourcing it runs DB derivation + the three production projects
R/run_project.R               # run_project(name): tar_make() + tar_prune() on one project
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
- `7z` (p7zip) on `PATH`, to extract the EUKARYOME general FASTA archives
- ITSx 1.1.3 in a conda env named `itsxenv`
  (`conda create -n itsxenv -c conda-forge -c bioconda itsx`), for the ITS1
  extraction of the ASVs (ROADMAP Q4)
- `vsearch` (>= 2.31) and BLAST+ (`blastn`, `makeblastdb`) on `PATH`; install
  them from your package manager or from source, they are not shipped in this
  repository
