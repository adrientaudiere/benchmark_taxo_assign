# 05 — Analysis notebook and documentation drift

## `analysis/benchmark.qmd` (1 381 lines) [S3.3]

Reading order of the file today:

| Lines | Content | State |
|---|---|---|
| 1–100 | setup, `load_all()`s, `tar_read()`, NA cleaning, `Gen_sp` columns | needed by everything below |
| 104–160 | consensus columns (Q3.1, Q3.3) | needed by metrics |
| 165–475 | NA proportion plots, `res_comp_tax` build, M1.1 enrichments, exploratory plots | half exploratory |
| 479–1075 | "Publication figures": Q1.2, Q1.3, Q1.4, Q3.2, Q3.3, M2.1, M2.2, Q2.1, Q2.2, Q2.5, Q3.5, Q1.6, Q1.7 | the manuscript |
| 1077–1200 | "En chantier": distance-to-true-community, UpSet plots, empty RMSE heading | sandbox |
| 1201–1381 | a *second* compute-cost analysis parsing the autometric log by regex (predates `benchmark_costs`), ends with `?metrica::confusion_matrix` | dead / duplicated |

Problems this causes in daily work:

- A render fails at the first `ggsave()` (B3) and at the sandbox section
  (`mat_com_true` uses `Gen_sp_true` then `distinct_species_taxo_mock$Species`,
  a column that no longer exists after the `rename(Gen_sp = Species)` at
  line 47). So the file is only ever run chunk by chunk, which means the
  "publication figures" are never regenerated from a clean session.
- The section order does not follow the ROADMAP (Q1 → Q2 → Q3 → M1 → M2);
  Q2 figures sit after M2.
- The regex cost analysis at the end and `benchmark_costs` compute the same
  quantities two ways; only one can be the manuscript's.

Proposed split, one Quarto project:

```
analysis/
├── _quarto.yml
├── 00_setup.R                 # libraries, load_pqverse(), theme, dir.create("figures")
├── 01_load_and_clean.qmd      # tar_read, NA cleaning, Gen_sp, consensus columns -> saves res_comp_tax.rds
├── 02_q1_methods.qmd          # Q1.2 Q1.3 Q1.4 Q1.6 Q1.7
├── 03_q2_databases.qmd        # Q2.1 Q2.2 Q2.5
├── 04_q3_consensus.qmd        # Q3.2 Q3.3 Q3.5
├── 05_m2_costs.qmd            # M2.1 M2.2 (after B1 is fixed)
├── 06_cross_validation.qmd    # D1a figures (currently a code snippet in the ROADMAP)
├── in_silico_simulation.qmd   # D1b (moved from the root)
└── sandbox/
    ├── distance_to_true_community.qmd   # lines 1077–1200
    └── autometric_regex_costs.qmd       # lines 1201–1381, or delete
```

`_quarto.yml` minimal content:

```yaml
project:
  type: default
  execute-dir: project
  render:
    - "0*.qmd"
    - "in_silico_simulation.qmd"
execute:
  freeze: auto
```

`execute-dir: project` makes `here::here()`, `tar_read(store = ...)` and
`ggsave("figures/...")` resolve from the project root without `setwd()`.
`freeze: auto` means chapter 02 re-renders without re-running chapter 01
unless 01 changed. Passing data between chapters through one
`rds`/`qs` file per chapter (`analysis/_cache/res_comp_tax.qs`) is the
simplest option; wrapping the notebook computations as targets
(`tarchetypes::tar_quarto()`) is the heavier option and not needed yet.

Smaller points in the notebook:

- `library("MiscMetabar")` then `load_all()` (critique 02).
- `here::i_am()` twice, `setwd()` once (critique 02).
- Line 774 still carries a commented `EUK_SSU_v1_9_3_Fungi` row in
  `db_meta`; `db_meta` belongs in `config.R` with `db_list` (critique 04).
- `conflicted::conflicts_prefer(dplyr::filter)` is declared twice (setup and
  line 1207). Keep one, in `00_setup.R`.
- The figure sizes (`width = 16, height = 7`) are repeated 26 times. A
  `save_fig(plot, name, width, height)` helper in `00_setup.R` that writes
  both pdf and png in one call halves the boilerplate and gives one place to
  change dpi.

## Documentation drift [S4.2]

Statements that no longer match the code, with the correct fact. All four
files are otherwise valuable; none of this is about deleting them.

| File | Says | Actually |
|---|---|---|
| README.md | "nine reference databases" | six active in `values_map` (four SSU variants commented out; decision 7). Nine was the pre-decision count. |
| README.md | "Four assignment methods … plus optional idtaxa" (via CLAUDE.md) | idtaxa is out (decision 1). |
| README.md, Quick start | `Rscript tests/test_combine_taxo_assignments.R` | fine today; becomes `tests/run_all.R` with S4.1. |
| README.md, Dependencies | vsearch on PATH | also BLAST+ on PATH; vsearch tarball tracked in repo. |
| CLAUDE.md | Section "`script_assign_taxo.R` → `store_assign_taxo`" describes the sequential `previous_target` chain as the current design | the chain is archived; the parallel script is the only one on `make.R`'s path. The section should describe `script_assign_taxo_parallel.R` and the sequential one should be one sentence under "Other notes". |
| CLAUDE.md | "`make.R` calls the first two in sequence" / "The cross-validation pipeline is run independently" | `make.R` runs DB derivation + all three pipelines. |
| CLAUDE.md | `cv_fold_tested` default 2, `cv_max_seq` default 100 | both are at their publication values (`cv_fold_number`, `NULL`). |
| CLAUDE.md | "`script_dada2.R` still has inline copies [of config.R] and should be migrated" (twice) | migrated (C3 ticked); it sources `config.R`. |
| CLAUDE.md | "`script_dada2.R` sources MiscMetabar's `R/` directly" | the directory no longer exists; the line is a silent no-op (B4). |
| CLAUDE.md | `EUK_SSU_99_*` variants, `v1_9_3` names, "4 methods × 9 DBs × 2 variants" for CV | v2 names; no 99% clustered variant is derived any more (`derive_clustered()` is defined but never called); 6 DBs. |
| CLAUDE.md | "`R/functions.R` used to mix three concerns…" paragraph | history, not guidance; move to Archives/README or drop. |
| CLAUDE.md | `mini_db` flag not mentioned | it is the single most consequential switch in `config.R` (B2). |
| CONTEXT.md | "nine reference-database variants" | six. |
| CONTEXT.md | "`make.R` runs the first two in sequence; `script_cross_val.R` is run independently" | see above. |
| CONTEXT.md | "`script_dada2.R` still has inline copies" | see above. |
| CONTEXT.md | `EUK_SSU_v1_9_3_Fungi_cut` in decision "SSU database reduced" | `EUK_SSU_v2_Fungi_cut`. |
| ROADMAP.md (old) | Snapshot table "4 methods × 7 DBs" | six DBs. |
| ROADMAP.md (old) | Q2 intro lists `v1_9_3` names; D1a says "9 DBs" | v2 names; six. |
| ROADMAP.md (old) | M2.1, M2.2, M2.3, Q2.5 ticked | built on stale cost data (B1); re-tick after S1.1. |
| ROADMAP.md (old) | C3 ticked but C1 open, C3's text still says "next time you re-execute" | C1 (byte-for-byte check of dbpq derivations) is still open and is the real blocker for a rerun. |
| config.R header | "script_dada2.R still has inline copies" | see above. |
| config.R line 60 | `cv_fold_tested <- cv_fold_number  # smoke-test default` | it is the publication value. |
| make.R | "Step 3: … (independent pipeline, run separately)" then runs it | drop the parenthesis or drop the step. |
| proposals_for_dbpq.md item 4 | "`derive_mini()` uses `head -n 10000`" | it now counts records with awk (`count > n {exit}`), so the truncation caveat is already fixed locally; update the proposal text. |

Suggested rule going forward: CLAUDE.md describes *what is on the
`make.R` path today*; anything historical goes to `Archives/README.md`.
That keeps CLAUDE.md short enough that drift is visible on a re-read.
