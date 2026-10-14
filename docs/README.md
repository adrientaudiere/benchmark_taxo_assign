# docs/

What each document is for, and how current it is. Read `objectives_design.md`
first: it is the source of truth for the design — a new design fact is defined
there and nowhere else. `ROADMAP.md` is the todo, `HISTORY.md` the archive of
what the code already did about it (§0).

| File | What it holds | State |
|---|---|---|
| `objectives_design.md` | The design, objective by objective: input data, negative controls, databases, methods, parameters and metrics, with the developer's decisions 1–27 (§4) and the files that still contradict them (§6) | **Current.** Source of truth |
| `hleap_2021_metrics.md` | Hleap et al. 2021 read from the paper and from their code, compared with comparpq; the per-unit truth of the two mocks and the scoring rules adopted here (§5–6) | **Current.** Source of truth for the scoring |
| `experimental_design.md` | Measured compute cost (§3), the ITS region of each dataset and what it constrains (§5.1, §5.4, §5.5), the cost structure and the levers (§8) | **Updated 2026-09-22** (item 0.10): §1, §2 and §4 read the current grid. **Its measured costs of §3 were read on the former 44-computation grid** and describe less than half the current work |
| `reference_databases.md` | Where each reference database comes from, how it is converted, what to change for a new release | **Current** |
| `diagnostics_cv_and_databases.md` | Measurements of 2026-09-15/16 on the production stores: cross-validation reference design, database effects. Kept for the manuscript's Methods and Discussion | **Current**; the cross-validation design decision it documents is still open |
| `approach_schema.csv`, `approach_schema.svg`, `approach_schema_use1.svg` … `_use3.svg` | The approach as a graph: the three uses (OD 51), the predictions (OD 45), the technical questions, the assignment grid, the datasets with their controls and truth, the scoring rules and the metrics, with the dada2 overflow and the parked cross-validation. One CSV row per relation, with its status (built / to compute / parked), the uses it is drawn for (`all` or their numbers) and the decision or chapter behind it. `approach_schema.svg` draws every row; `approach_schema_use<k>.svg` draws the rows of use *k*, without the predictions and the output. `Rscript tools/approach_schema.R` redraws the four SVG | **Current** (2026-10-13); edit the CSV when a decision changes a relation |
| `mock_scoring_matrix.csv`, `mock_scoring_metrics.csv`, `mock_scoring_table.csv`, `mock_scoring.svg` | How a mock is scored: (1a) the truth matrix, the cell (TP, FP, FN, TN or left out) of each kind of unit — real with a truth, tied, `no_truth`, `foreign`, `fake_`, `external_` by database — at one rank; (1b) the metrics computed from the cells (OD 49: headline on the real units); (2) the extended mock table, the phyloseq of a mock input with its controls and the truth joined to it. Follows `comparpq::unit_cells()` and `R/score_assignments.R`; `Rscript tools/mock_scoring_schema.R` redraws the SVG | **Current** (2026-10-13); edit the CSV when a decision changes a scoring rule |
| `design_analysis.svg` | Diagram of the analysis, February 2025 | **Outdated** (predates the 2026-09 design); kept as a visual trace, replaced by `approach_schema.svg` |
| `references/` | External documents consulted, not read by any script | See its own README |

Deleted on 2026-09-22: `compute_budget.md` (its 2026-09-11 estimates were
replaced by the measurements of `experimental_design.md` §3; its cost structure
and lever table moved to §8 of that file), and the rendered `*.html` pages with
their `*_files/` library folders, which are regenerated from the Markdown with
`quarto render docs/<file>.md` and are ignored by git.
