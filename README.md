# benchmar_assign_taxo

Benchmark of taxonomic-assignment methods × reference databases for ITS fungal
metabarcoding, supporting a manuscript-in-preparation in the
[pqverse](https://github.com/adrientaudiere/pqverse) ecosystem.

The question: **which method × reference database (including database
modifications) × parameter combinations give the best trade-off between
assignment quality, compute time and compute power** for the taxonomic
assignment of Fungi. The benchmark tests that one step and nothing else: no
community ecology.

## Design (decided 2026-09-17; inputs since OD 34)

The full design, objective by objective, is in
[`docs/objectives_design.md`](docs/objectives_design.md); the scoring in
[`docs/hleap_2021_metrics.md`](docs/hleap_2021_metrics.md). In short:

| Axis | Values |
|---|---|
| Methods | `dada2`, `sintax`, `lca`, `blastn` |
| Parameters | dada2 / sintax bootstrap ∈ {0.4, 0.5, 0.6}; lca `lca_cutoff` ∈ {0.51, 0.8, 1}; blastn vote (relative, absolute majority, unanimity) × `min_id` ∈ {90, 92, 95, 97} — derived from one computation per method × database × input |
| Databases | EUKARYOME ITS v2.1, UNITE 19.02.2025 without singletons, UNITE with singletons × {all, Fungi only, Fungi + 5 % non-fungal representatives} = 9 databases |
| Inputs | ITSx inputs only (OD 34): each ITS region of the dataset's amplicon (`config.R::bio_datasets`) × ASVs and their 97 % OTUs; grid sizes per project in [`ETAT.md`](ETAT.md) |
| Negative controls | shuffled sequences (`fake_`) and 100 non-fungal sequences (`external_`) added to every input |
| Datasets | mock communities of Pauvert et al. 2019 (ITS1), Tedersoo et al. 2026 (full ITS) and Rolling et al. 2022 (ITS1), scored per unit against the sequence of each strain as in Hleap et al. 2021; in silico reads (InSilicoSeq); biological datasets (air, wood) read through their share of unassigned sequences and controls; cross-validation parked (OD 28) |
| Endpoints | F1 and PPV at Genus (OD 43); Species in its own section, by exact and near match (chapter 09) |
| Cost | wall time, peak memory and CPU per computation, child processes included (`autometric`, `R/sample_tree.R`), CO₂eq (`greenAlgoR`) |

Trimming the references to the amplicon (cutadapt) is left out of the grid and
will be tried at the end on the best combinations, on the Pauvert mock only.
Consensus voting across columns is compared with
`comparpq::resolve_taxo_conflict()`.

**Status.** What has run and what is left: [`ROADMAP.md`](ROADMAP.md); the
dated counts (targets, outdated targets, units per input, tests):
[`ETAT.md`](ETAT.md).

## Quick start

```r
# Reference databases (idempotent; see docs/reference_databases.md)
source("make_databases.R")
download_reference_sources()   # general FASTA releases listed in config.R
derive_all_variants()          # dada2 / sintax formats, Fungi, Fungi + rep and cut variants

# Everything in order: databases, then every dataset (a day or more)
source("make.R")
# At night, per dataset (its DADA2 store, then its grid; log in logs/):
#   Rscript run_bio_night.R pauvert korhonen_wood ... species

# One project at a time (sourcing R/run_project.R runs nothing by itself)
source("R/run_project.R")
run_project("assign_taxo")        # tar_make() + tar_prune() on that project

# Smoke test on the mini_* databases (own stores, never overwrite production)
run_project("assign_taxo_mini")

# Tests (every file, or one at a time with Rscript tests/test_<name>.R)
Rscript tests/run_all.R
```

The analysis is a Quarto project in `analysis/` (`quarto render analysis`, or
one chapter at a time starting with `01_load_and_clean.qmd`, which feeds the
others through `analysis/_cache/`). Figures are written to `figures/`
(not tracked). `BENCHMARK_ANALYSIS_MINI=TRUE` renders the chapters on the
smoke-test stores, with their own cache and figures.

## Documentation

- **`docs/objectives_design.md`** — the decided design: for each objective,
  the inputs, controls, databases, methods, parameters and metrics, with the
  developer's decisions 1–27.
- **`docs/hleap_2021_metrics.md`** — how Hleap et al. 2021 score an
  assignment (paper and code), the per-unit truth of the two mocks, and the
  scoring rules adopted here.
- **`docs/experimental_design.md`** — the ITS-region constraints of each dataset
  and the measured compute cost (updated to the nine-database grid on
  2026-09-22, but its cost figures were read on the former 44-computation grid
  and describe less than half the current work).
- **`AGENTS.md`** — architecture, pipeline layout, external dependencies, and
  the framing of the manuscript. `CLAUDE.md` points to it.
- **`CONTEXT.md`** — glossary of the project's terms and key decisions.
- **`ROADMAP.md`** — the todo, one line per open item: what blocks the first
  production run, the manuscript framing, the open scientific work, the code
  debt and the remaining structural work.
- **`HISTORY.md`** — the archive: every ticked item with its dated narrative and
  its measurements, and every past decision (sections 3 to 3d).
- **`ETAT.md`** — the dated figures (targets per project, outdated and
  orphan targets, computations, units per input, tests), newest first.
  Written by `tools/etat.R`, never by hand, at the end of every
  `run_bio_night.R`; the other documents point to it for counts.
- **`reviews/`** — the review documents: `code_review_2026-10-05.md` (code
  review, its fix plan in §10 and the developer's answers in §11),
  `critique_2026-09-22.md` and `literature_2026-09-22.md` (the aims, the
  code size and the six benchmarks to cite; the `C22-…` items),
  `critique_fable/` (evidence for the structural `S…` items, 2026-09-10) and
  `proposals_for_dbpq.md` (helpers that could move to
  [`dbpq`](https://github.com/adrientaudiere/dbpq), and the dbpq defects
  worked around here).
- **`docs/reference_databases.md`** — where each reference database comes
  from (UNITE and EUKARYOME releases, URLs, DOIs), how it is converted, and
  what to change when a new release comes out.
- **`docs/README.md`** — index of `docs/`: what each document holds and how current it is (also the cross-validation diagnostics, the database procedure, the design diagram and the external references).

## Layout

```
_targets.yaml                 # named targets projects (dada2, assign_taxo, *_mini, per-dataset pairs, species)
pipelines/dada2.R             # ASV pipeline (Pauvert mock) -> store_dada2
pipelines/dada2_bio.R         # ASV pipeline of the other fastq datasets -> store_dada2_<dataset>
pipelines/load_bio.R          # published OTUs (air, wood) -> store_dada2_<dataset>
pipelines/assign_taxo.R       # methods x DBs x ITSx inputs -> store_assign_taxo[_<dataset>]
pipelines/species.R           # species by exact and near match -> store_species
make.R                        # sourcing it runs DB derivation + every dataset's projects
run_bio_night.R               # night runner: per dataset, its DADA2 store then its grid
R/run_project.R               # run_project(), dataset_projects(), production_datasets()
make_databases.R              # idempotent, atomic DB derivation (delegates to dbpq::)
config.R                      # shared constants (primers, threads, pinned settings, databases, datasets)
R/                            # helpers, one file per concern (see AGENTS.md)
analysis/                     # Quarto project: 00_setup.R + chapters 01-05, 07-09, sandbox/
logs/                         # run logs (not tracked)
figures/                      # save_fig() output of the chapters (not tracked)
tests/                        # Rscript-runnable testthat fixtures
docs/                         # design documents, database procedure, external references
tools/etat.R                  # writes the dated figures in front of ETAT.md
reviews/                      # review documents (code review, critiques, dbpq proposals)
Archives/                     # superseded code (cross-validation, pre-OD 34 sandbox; kept for reference)
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
  (`conda create -n itsxenv -c conda-forge -c bioconda itsx`), for the ITS
  extraction of the ASVs (ROADMAP Q4)
- `vsearch` (>= 2.31) and BLAST+ (`blastn`, `makeblastdb`) on `PATH`; install
  them from your package manager or from source, they are not shipped in this
  repository
