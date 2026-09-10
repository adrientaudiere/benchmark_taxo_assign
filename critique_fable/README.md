# critique_fable — structural review of `benchmar_assign_taxo`

Review date: 2026-09-10. Reviewer: Claude Fable 5.1 (session
`session_01StoovfuhsMZbyTeW8QZvuf`). Scope: project structure, pipeline
wiring, analysis notebook, documentation. Not reviewed: the scientific
content of the benchmark, nor the internals of the sister pqverse packages.

Every item here is mirrored as a ticket in `../ROADMAP.md` (section
"Structure & maintenance"). The ROADMAP is the todo list; these files hold
the evidence and the proposed fix so the ROADMAP bullets can stay short.

## Files

| File | What it covers |
|---|---|
| [01_bugs_and_traps.md](01_bugs_and_traps.md) | Confirmed bugs and silent traps found by reading the code and the on-disk artefacts (stores, logs, data). Start here. |
| [02_dependencies_and_paths.md](02_dependencies_and_paths.md) | How the sister packages are loaded (four different ways), hard-coded absolute paths, and a single-entry-point proposal. |
| [03_layout_and_files.md](03_layout_and_files.md) | Proposed directory tree, stray files at the root, git hygiene, targets multi-project config. |
| [04_pipeline_code.md](04_pipeline_code.md) | Duplication (`values_map` in three places), naming, store hygiene, test coverage. |
| [05_analysis_and_docs.md](05_analysis_and_docs.md) | `benchmark.qmd` split plan, Quarto project setup, and a table of statements in README/CLAUDE/CONTEXT/ROADMAP that no longer match the code. |
| [06_tidypq_adoption.md](06_tidypq_adoption.md) | Every hand-written phyloseq slot manipulation and the tidypq verb that replaces it (decision 13, items S5.1–S5.2). |

## How the evidence was gathered

- Read every tracked file except `Archives/` internals and the two
  `store_*/meta/meta` blobs (only their target-name column was inspected).
- Parsed `data/data_final/autometric_log_assign_taxo.txt` (47 512 rows) to
  cross-check the `benchmark_costs` target against the timestamps of the
  last pipeline run recorded in `store_assign_taxo/meta/meta`.
- Listed `data/data_raw/` to compare the reference-database files on disk
  with the names used in the scripts and the docs.
- Checked installed package versions (`ggplot2 4.0.3`, `targets 1.12.0`,
  `autometric 0.1.2`) where a behaviour depends on the version.

## Severity legend (same tags as the workspace ROADMAP)

`[Critical / easy]` etc. Priority: Critical, High, Medium, Low. Facility:
easy (< 1 day), moderate (a few days), hard (architectural / external).

"Confirmed" means the failure was reproduced from on-disk evidence in this
session. "Suspected" means the code reads as wrong but no run was done to
prove it.
