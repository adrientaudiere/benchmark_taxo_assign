# Seed spread of the k-mer methods on Pauvert (ROADMAP §1, code review §1.9,
# developer 2026-10-06): how much do dada2 and sintax move when only their
# random seed changes? This is the noise floor of every small effect of the
# Results, and of the read-filtering choice (test_read_filtering.R).
#
#   input:     d_ITS1_asv_for_assignation of store_assign_taxo (read only)
#   databases: the two of the zoom (EUK_ITS_v2.1_Fungi_rep, Unite_all_20250219)
#              and Unite_all_20250219_Fungi, where §1.9 measured the spread
#   methods:   dada2 (set.seed(), 3 threads) and sintax (--randseed, 1 thread);
#              lca and blastn draw nothing, their spread is zero
#   seeds:     replicate 0 is the stored computation (the target's seed);
#              replicates 1 to 3 rerun it with seeds 1, 2 and 3
#
# Replicate 0 is scored from the stored computation by the same path as the
# reruns, and checked against the cached res_comp_tax of chapter 01 before
# anything is reported.
#
# One raw result per computation x seed in data/snapshots/seed_spread/, so a
# rerun only computes what is missing. Run from the project root:
#   Rscript analysis/sandbox/test_seed_spread.R

suppressMessages({
  source(here::here("analysis/00_setup.R"))
  source(here("analysis/sandbox/seed_filter_helpers.R"))
})
assert_dada2_patched()

out_dir <- here("data/snapshots/seed_spread")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
seeds <- 1:3
dbs <- c(zoom_test_dbs, "Unite_all_20250219_Fungi")
methods <- c("dada2", "sintax")

store <- here("store_assign_taxo")
input_pq <- tar_read(d_ITS1_asv_for_assignation, store = store)
truth <- build_mock_truth(
  tar_read(d_asv, store = here("store_dada2")),
  pauvert_cfg
)

grid <- tidyr::expand_grid(method = methods, db = dbs, seed = c(0L, seeds)) |>
  mutate(
    compute_name = compute_name_for(method, db, "ITS1_asv"),
    file = file.path(out_dir, paste0(compute_name, "__seed", seed, ".rds"))
  )

score_replicate <- function(i) {
  g <- grid[i, ]
  pq <- derive_columns(input_pq, readRDS(g$file), g$method, g$db)
  cols <- pauvert_vm$full_name[pauvert_vm$compute_name == g$compute_name]
  score_columns(pq, cols, truth) |>
    mutate(seed = g$seed)
}

# ---- replicate 0: the stored computations, checked first ---------------------
for (i in which(grid$seed == 0L)) {
  if (!file.exists(grid$file[i])) {
    saveRDS(tar_read_raw(grid$compute_name[i], store = store), grid$file[i])
  }
}
scores_0 <- bind_rows(lapply(which(grid$seed == 0L), score_replicate))
# Replicate 0 must give the cached scores of chapter 01, value for value.
cached <- cache_read("res_comp_tax") |>
  filter(dataset == "pauvert", region == "ITS1", unit == "ASV") |>
  select(method_db, tax_level, metrics, external_scoring, cached = values)
check <- inner_join(
  scores_0,
  cached,
  by = c("method_db", "tax_level", "metrics", "external_scoring")
)
stopifnot(
  "replicate 0 rows missing from res_comp_tax" = nrow(check) == nrow(scores_0),
  "replicate 0 differs from res_comp_tax" = isTRUE(all.equal(
    check$values,
    check$cached
  ))
)
message("✓ replicate 0 reproduces res_comp_tax (", nrow(check), " values)")

# ---- replicates 1 to 3 -------------------------------------------------------
for (i in which(grid$seed != 0L)) {
  g <- grid[i, ]
  if (file.exists(g$file)) {
    next
  }
  message(format(Sys.time(), "%H:%M"), " ", g$compute_name, " seed ", g$seed)
  saveRDS(run_computation(input_pq, g$method, g$db, g$seed), g$file)
}
scores <- bind_rows(
  scores_0,
  bind_rows(lapply(which(grid$seed != 0L), score_replicate))
)
write.csv(scores, file.path(out_dir, "scores.csv"), row.names = FALSE)

# ---- spread ------------------------------------------------------------------
spread <- scores |>
  default_rows() |>
  endpoint_rows() |>
  group_by(method, db, tax_level, metrics) |>
  summarise(
    stored = values[seed == 0L],
    min = min(values),
    max = max(values),
    range = max - min,
    sd = sd(values),
    .groups = "drop"
  ) |>
  arrange(tax_level, metrics, method, db)
write.csv(spread, file.path(out_dir, "spread_default.csv"), row.names = FALSE)
cat("\nSeed spread at default settings (4 replicates: stored + 3 seeds)\n")
print(mutate(spread, across(stored:sd, \(x) round(x, 4))), n = Inf, width = 200)

# Every threshold, to see whether the spread depends on it.
spread_all <- scores |>
  endpoint_rows() |>
  group_by(method, db, min_bootstrap, tax_level, metrics) |>
  summarise(range = max(values) - min(values), .groups = "drop")
write.csv(
  spread_all,
  file.path(out_dir, "spread_all_thresholds.csv"),
  row.names = FALSE
)
cat("\nLargest range over the thresholds 0.4 / 0.5 / 0.6, per endpoint\n")
print(
  spread_all |>
    group_by(tax_level, metrics, method) |>
    summarise(max_range = round(max(range), 4), .groups = "drop"),
  n = Inf
)

# Units whose call changes between replicates, at default settings, Genus.
calls <- bind_rows(lapply(which(grid$seed %in% c(0L, seeds)), function(i) {
  g <- grid[i, ]
  pq <- derive_columns(input_pq, readRDS(g$file), g$method, g$db)
  col <- paste0(
    "Genus_",
    full_name_for(g$method, g$db, 0.5, input = "ITS1_asv")
  )
  tibble(
    method = g$method,
    db = g$db,
    seed = g$seed,
    taxon = taxa_names(pq),
    genus = as.character(pq@tax_table[, col])
  )
}))
changed <- calls |>
  group_by(method, db, taxon) |>
  summarise(n_calls = n_distinct(genus, na.rm = FALSE), .groups = "drop") |>
  group_by(method, db) |>
  summarise(units = n(), changing_genus = sum(n_calls > 1), .groups = "drop")
cat("\nUnits whose Genus call (bootstrap 0.5) differs between replicates\n")
print(changed)
