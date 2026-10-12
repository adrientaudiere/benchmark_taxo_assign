source("analysis/00_setup.R")
source(here("R/mock_truth.R"))
source(here("R/score_assignments.R"))
source(here("R/assignment_helpers.R"))

# Evidence for ROADMAP C22-P21 (2026-09-30): Q3 voter diversity, at zero
# compute on the Pauvert ITS1 ASV input (OD 34 store).
#
# Hypothesis (reviews/literature_2026-09-22.md §4.3, critique addendum): a vote helps
# through error independence, so it gains from method-family and provider
# diversity, not from column count. Test: the `rel_majority` vote of chapter 01
# (Q3 base strategy, nb_agree_threshold 1, not strict) rerun on subsets of its
# 36 default columns (4 methods x 9 databases), each scored like chapter 01
# and compared with the best and the mean of its own members.
#
# Families: composition (k-mer) = dada2, sintax; search = blastn, lca.
# Providers: UNITE (`Unite_all` nested in `Unite_s_all`) and EUKARYOME.
# Metrics on the real units only (C22-P8): F1 = 2 TP / (2 TP + FP_real + FN).
# Check: the 36-column subset must reproduce chapter 01's
# `ITS1_asv_rel_majority_consensus` (analysis/_cache/res_comp_tax.rds).

input <- "ITS1_asv"
d <- tar_read_raw(paste0("d_all_taxo_", input), store = store_assign_taxo) |>
  clean_taxo() |>
  add_gen_sp()

vm <- values_map_for("pauvert")
members <- vm |>
  filter(.data$input == .env$input) |>
  filter(
    (method %in% c("dada2", "sintax") & min_bootstrap == 0.5) |
      (method == "lca" & lca_cutoff == 0.51) |
      (method == "blastn" & vote_algorithm == "rel_majority" & min_id == 95)
  ) |>
  transmute(
    full_name,
    algo = method,
    db = assignment_db(full_name),
    provider = case_when(
      startsWith(db, "Unite_s_all") ~ "Unite_s_all",
      startsWith(db, "Unite_all") ~ "Unite_all",
      startsWith(db, "EUK_ITS") ~ "EUK"
    )
  )
stopifnot(nrow(members) == 36, !anyNA(members$provider))

pick <- function(
  algos = single_methods,
  providers = c("Unite_s_all", "Unite_all", "EUK")
) {
  members$full_name[members$algo %in% algos & members$provider %in% providers]
}

designed <- list(
  all_36 = pick(),
  fam_composition = pick(c("dada2", "sintax")),
  fam_search = pick(c("blastn", "lca")),
  cross_dada2_blastn = pick(c("dada2", "blastn")),
  cross_dada2_lca = pick(c("dada2", "lca")),
  cross_sintax_blastn = pick(c("sintax", "blastn")),
  cross_sintax_lca = pick(c("sintax", "lca")),
  prov_unite_pair = pick(providers = c("Unite_s_all", "Unite_all")),
  prov_Unite_s_all_EUK = pick(providers = c("Unite_s_all", "EUK")),
  prov_Unite_all_EUK = pick(providers = c("Unite_all", "EUK")),
  prov_EUK_only = pick(providers = "EUK"),
  prov_Unite_all_only = pick(providers = "Unite_all"),
  prov_Unite_s_all_only = pick(providers = "Unite_s_all")
)

# Random subsets of the same sizes, to read the effect of column count alone.
n_draws <- 30
random <- withr::with_seed(20260930, {
  unlist(
    lapply(c(12, 18, 24), function(k) {
      stats::setNames(
        lapply(seq_len(n_draws), \(i) sample(members$full_name, k)),
        sprintf("random_%02d_%02d", k, seq_len(n_draws))
      )
    }),
    recursive = FALSE
  )
})
subsets <- c(designed, random)

# One rel_majority vote per subset, on a phyloseq holding only its columns.
otu <- otu_table(d)
taxtab <- as(tax_table(d), "matrix")
vote <- function(cols, label) {
  rank_cols <- as.vector(outer(tax_order, cols, paste, sep = "_"))
  sub <- phyloseq(otu, tax_table(taxtab[, rank_cols, drop = FALSE]))
  new_names <- paste0(tax_order, "_", label)
  out <- resolve_taxo_conflict(
    sub,
    pattern_tax_ranks = tax_order,
    method = "rel_majority",
    strict = FALSE,
    nb_agree_threshold = 1L,
    new_names = new_names,
    replace_collapsed_rank_by_NA = TRUE
  )
  as(tax_table(out), "matrix")[, new_names, drop = FALSE]
}
labels <- paste0("vote_", names(subsets))
votes <- do.call(cbind, Map(vote, subsets, labels))
tax_table(d) <- tax_table(cbind(taxtab, votes))

cfg <- mock_sanger_tables[mock_sanger_tables$dataset == "pauvert", ]
truth <- build_mock_truth(tar_read(d_asv, store = here("store_dada2")), cfg)
ranks_df <- build_ranks_df(d, labels, ranks = tax_order)
scores <- score_assignments(
  d,
  ranks_df,
  truth,
  db_meta,
  external_truth = external_control_truth(taxa_names(d)),
  truth_ranks = truth_ranks,
  consensus_scorings = "aside",
  no_truth = truth_rules[[lead_truth_rule]]
)

real_only <- function(x) {
  x |>
    filter(tax_level %in% c("Genus", "Species", "Gen_sp")) |>
    mutate(tax_level = ifelse(tax_level == "Genus", "Genus", "Species")) |>
    select(method_db, tax_level, metrics, values) |>
    distinct() |>
    pivot_wider(names_from = metrics, values_from = values) |>
    mutate(
      FP_real = misassign_seq * n_real,
      PPV_real = TP / (TP + FP_real),
      F1_real = 2 * TP / (2 * TP + FP_real + FN)
    ) |>
    select(method_db, tax_level, TP, FP_real, FN, PPV_real, TPR, F1_real)
}

vote_scores <- real_only(scores) |>
  mutate(subset = sub("^vote_", "", method_db))
single <- readRDS("analysis/_cache/res_comp_tax.rds") |>
  keep_truth_rule() |>
  filter(
    dataset == "pauvert",
    unit == "ASV",
    method_db %in% members$full_name
  ) |>
  real_only()

# Check against chapter 01.
chapter01 <- readRDS("analysis/_cache/res_comp_tax.rds") |>
  keep_truth_rule() |>
  filter(
    dataset == "pauvert",
    unit == "ASV",
    method_db == paste0(input, "_rel_majority_consensus"),
    external_scoring == "aside"
  ) |>
  real_only()
mine <- filter(vote_scores, subset == "all_36") |> arrange(tax_level)
theirs <- arrange(chapter01, tax_level)
stopifnot(
  nrow(theirs) == 2,
  isTRUE(all.equal(mine$F1_real, theirs$F1_real)),
  isTRUE(all.equal(mine$TP, theirs$TP))
)
stopifnot(
  isTRUE(all.equal(mine$FP_real, theirs$FP_real)),
  isTRUE(all.equal(mine$FN, theirs$FN))
)
message(
  "all_36 reproduces chapter 01's rel_majority consensus: OK"
)

members_of <- tibble::enframe(subsets, name = "subset", value = "full_name") |>
  tidyr::unnest(full_name)
summary <- members_of |>
  inner_join(
    single,
    by = c("full_name" = "method_db"),
    relationship = "many-to-many"
  ) |>
  group_by(subset, tax_level) |>
  summarise(
    n_col = n(),
    best_single = max(F1_real),
    mean_single = mean(F1_real),
    .groups = "drop"
  ) |>
  inner_join(
    select(vote_scores, subset, tax_level, PPV_real, TPR, F1_real),
    by = c("subset", "tax_level")
  ) |>
  mutate(
    gain_vs_best = F1_real - best_single,
    gain_vs_mean = F1_real - mean_single
  )

show <- function(x) {
  x |>
    mutate(across(where(is.numeric) & !any_of("n_col"), \(v) round(v, 3))) |>
    as.data.frame() |>
    print()
}
for (lvl in c("Genus", "Species")) {
  message("\n== ", lvl, ": designed subsets (F1 on real units)")
  summary |>
    filter(tax_level == lvl, !startsWith(subset, "random_")) |>
    select(
      subset,
      n_col,
      F1_real,
      best_single,
      mean_single,
      gain_vs_best,
      PPV_real,
      TPR
    ) |>
    show()
  message("\n== ", lvl, ": random subsets, ", n_draws, " draws per size")
  summary |>
    filter(tax_level == lvl, startsWith(subset, "random_")) |>
    group_by(n_col) |>
    summarise(
      F1_mean = mean(F1_real),
      F1_min = min(F1_real),
      F1_max = max(F1_real),
      gain_vs_best_mean = mean(gain_vs_best),
      share_beating_best = mean(gain_vs_best > 0),
      .groups = "drop"
    ) |>
    show()
}

message("\n== Best single columns (F1 on real units)")
single |>
  group_by(tax_level) |>
  slice_max(F1_real, n = 3) |>
  select(method_db, tax_level, F1_real) |>
  show()
