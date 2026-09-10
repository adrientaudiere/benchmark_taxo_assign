# 06 — tidypq adoption: call-site inventory

Decision 13 (2026-09-10): the benchmark adopts `tidypq` in the notebook
(S5.1) and in the pipelines (S5.2). This file lists every place where a
phyloseq slot is manipulated by hand today and the tidypq verb that is the
natural replacement. Exported verbs as of 2026-09-10 (from
`pqverse_pkg/tidypq/NAMESPACE`):

```
pq_to_tidy  tax_table_to_df  taxa_prevalence
filter_samples_pq  filter_taxa_pq  filter_occurrences_pq  filter_tree_pq
mutate_samdata_pq  mutate_taxa_pq  mutate_occurrences_pq
select_samdata_pq  select_taxa_pq  rename_samples_pq  rename_taxa_pq
arrange_samples_pq  arrange_taxa_pq  slice_samples_pq  slice_taxa_pq
merge_markers_pq  (+ contaminant / chimera helpers, not relevant here)
```

Argument names below are the *intended* mapping; check each signature in
`pqverse_pkg/tidypq/R/` before editing, the package is still in development.

## S5.1 — `analysis/benchmark.qmd`

| Lines | Today | Replacement | Note |
|---|---|---|---|
| 55–59 | `colnames(d_all_taxo@tax_table)[grep("^Species_", ...)]` to pick the ranks fed to `taxtab_replace_pattern_by_NA()` | `tax_table_to_df(d_all_taxo) |> select(starts_with("Species_")) |> names()` | comparpq's `taxtab_replace_pattern_by_NA()` stays; only the column discovery changes. |
| 73–103 | `d_all_taxo@tax_table |> as.data.frame()`, build `Gen_sp_*` columns with a `for` loop, then `d_all_taxo@tax_table <- tax_table(as.matrix(taxtab_d_all))` | `mutate_taxa_pq(d_all_taxo, across(...))` or `mutate_taxa_pq()` inside a loop over `cl_names`; no slot assignment | The round trip through `as.matrix()` silently coerces every column to character and drops the taxa order guarantee; `mutate_taxa_pq()` is the safe path. |
| 172–173 | `colSums(is.na(d_all_taxo@tax_table)) / ntaxa(d_all_taxo)` | `tax_table_to_df(d_all_taxo) |> summarise(across(everything(), ~ mean(is.na(.x)))) |> pivot_longer(...)` | Same numbers, tidy output ready for the ggplot that follows. |
| 199, 207 | `sum(grepl("external|fake", taxa_names(d_all_taxo))) / ntaxa(d_all_taxo)` | `ntaxa(filter_taxa_pq(d_all_taxo, grepl("^fake_|^external_", taxa_names))) / ntaxa(d_all_taxo)` | Check how `filter_taxa_pq()` exposes taxa names (a `.taxa` pronoun or a column of `tax_table_to_df()`). |
| 247–248 | `colnames(d_all_taxo@tax_table)` membership check | `names(tax_table_to_df(d_all_taxo))` | Trivial. |
| 1090 | `d_all_taxo |> psmelt() |> select(OTU, Sample, Abundance, contains("Gen_sp"))` | `pq_to_tidy(d_all_taxo, ranks = starts_with("Gen_sp"))` | This is the workspace-wide rule: never `psmelt()`. Part of the sandbox section; rewrite when it moves to `analysis/sandbox/`. |

Also in scope while touching the notebook:

- `simplify_taxo()` (line 69) is a MiscMetabar function; keep it.
- `In_silico_simulation.qmd` lines 204–205 build a phyloseq from matrices
  with `phyloseq::otu_table()` / `sample_data()`; that is construction,
  not manipulation — no tidypq equivalent, leave as is.

## S5.2 — pipelines and helpers

| File / function | Today | Replacement | Impact |
|---|---|---|---|
| `R/combine_taxo_assignments.R` | `as.matrix(unclass(pq@tax_table))`, `setdiff()` of column names, `cbind`, then `base_pq@tax_table <- tax_table(...)` | `tax_table_to_df()` on base and on each assignment, `dplyr::bind_cols()` of the new columns, one `mutate_taxa_pq()` (or a tidypq constructor if one exists) to write back | `tests/test_combine_taxo_assignments.R` pins the column order and values; it must stay green unchanged. Rebuilds `d_all_taxo`. |
| `R/functions.R::create_fake_pq_from_refseq()` | builds `tax_table()` / `otu_table()` by hand from the sintax headers | keep the construction (no tidypq verb builds a phyloseq from a fasta), but replace `taxa_names(x) <- ...` triple with one `rename_taxa_pq()` after `phyloseq()` | Used by every CV fold; rebuilds `store_cross_val`. |
| `R/cross_val.R` lines 151, 160, 173–177, 192, 243, 394 | `fake_pq@tax_table`, `as_tibble(as.matrix(unclass(...)))` | `tax_table_to_df(fake_pq)` | Six sites, same expression; a one-line helper is enough. |
| `script_dada2.R` lines 166–186 | `rename_samples(sample_data(...))`, `rename_samples(otu_table(seqtab[...]))` with a duplicated de-duplication predicate | build `d_asv` first, then `rename_samples_pq(d_asv, new_names = ...)` once; drop duplicated samples with `filter_samples_pq()` before renaming | Rebuilds `store_dada2` from `sam_tab`; already invalidated by S0.3, so no extra cost if done in the same rerun. |
| `script_assign_taxo_parallel.R` | `read.csv(...) |> select(...) |> magrittr::set_rownames(...)` for `taxo_mock` | not a phyloseq operation; leave | – |

## Sequencing with the other items

1. S2.2 (loader) adds `tidypq` to `load_pqverse()`.
2. S5.2 lands with the production rerun that S0.1 already forces, so the
   stores are rebuilt once, not twice.
3. S5.1 lands with the S3.3 notebook split; each chapter is rewritten once.
4. `tests/test_combine_taxo_assignments.R` is the regression guard for
   S5.2; add `tests/test_create_fake_pq_from_refseq.R` (S4.1) before
   touching that function.
