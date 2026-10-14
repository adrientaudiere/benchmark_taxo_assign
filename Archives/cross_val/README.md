The cross-validation axis, parked by OD 28, archived on 2026-10-07 (code
review §8) so that nobody reads its stale figures as current. Files keep
their project-relative layout: `R/cross_val.R`, `R/examples_cross_val.R`,
`R/cv_to_tidy.R`, `pipelines/cross_val.R`, `run_cv_night.R`, chapter
`analysis/06_cross_validation.qmd`, its test and two sandbox scripts.

`store_cross_val` (7-database grid, 2026-09-15/16) stays at the project root,
untouched; the `cross_val` and `cross_val_mini` projects left `_targets.yaml`.
Moved here on 2026-10-14 (audit of that day), as nothing outside the CV
called them: `R/cv_queries.R`, `R/create_fake_pq_from_refseq.R`, the grid
builder `build_cv_values_map()` (now `R/cv_values_map.R`), the `cv_*`
constants of `config.R` (now `config_cv.R`, sourced after `config.R`), their
tests and the six CV diagnostics of `analysis/sandbox/`. To restore the
axis, copy the tree back to the project root.

Known defects, not fixed (review §8): FP and FN swapped in
`compute_by_tax_level`, Kingdom dropped, NA answers never FN, cumulative
NA-masking for a non-ascending `min_bootstrap`. When OD 28 is revisited,
rewrite `cross_val()` around `comparpq::tc_metrics_unit()`, as the mocks are
scored, rather than restoring these files.
