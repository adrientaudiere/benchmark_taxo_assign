source("analysis/00_setup.R")
source(here("R/mock_truth.R"))
source(here("R/assignment_helpers.R"))
source(here("R/score_assignments.R"))

# Evidence for ROADMAP SA4 (2026-09-30): what the `*_Incertae_sedis` labels of
# the external-control truth do to `ext_correct`, and what each fix would
# change. Pauvert ITS1 ASV, default settings, the three `full` databases (the
# only ones where `ext_correct` is read below Kingdom; on `Fungi_rep` it is
# Kingdom only, which holds no placeholder).
#
# comparpq::tc_metrics_unit() computes
#   ext_correct = #(controls given their own label) / #(all external controls),
# so a control whose truth is a placeholder (erased on the assignment side by
# clean_taxo()) can never be correct but stays in the denominator. Variants:
#   A  before 2026-09-30: raw truth, denominator = all controls;
#   B  the fix proposed on 2026-09-23: truth truncated from its first
#      placeholder down, denominator unchanged;
#   C  B, denominator = controls with a truth at that rank (the rule of the
#      mock units, left out below their truth depth);
#   D  placeholder set to NA at its own rank only (real names below it kept),
#      denominator = controls with a truth at that rank.

input <- "ITS1_asv"
d <- tar_read_raw(paste0("d_all_taxo_", input), store = store_assign_taxo) |>
  clean_taxo() |>
  add_gen_sp()
tt <- as.data.frame(as(tax_table(d), "matrix"), stringsAsFactors = FALSE)
ext <- grepl("^external_", rownames(tt))
placeholder <- function(v) {
  !is.na(v) & grepl("_incertae_sedis$|_Incertae_sedis$|^unclassified|^$", v)
}
# The raw UNITE labels: external_control_truth() has set the placeholders to
# NA since decision D (2026-09-30), so it is called here with a pattern that
# matches nothing.
placeholder_label_d <- placeholder_label
placeholder_label <- "a^"
truth <- external_control_truth(rownames(tt)[ext])
placeholder_label <- placeholder_label_d
stopifnot(identical(
  external_control_truth(rownames(tt)[ext])$Genus,
  ifelse(placeholder(truth$Genus), NA_character_, truth$Genus)
))
stopifnot(identical(truth$taxon, rownames(tt)[ext]))

ranks <- truth_ranks

message("External controls: ", nrow(truth))
message(
  "Controls with a placeholder at each rank, and with a real name below one"
)
first_ph <- apply(
  vapply(ranks, \(r) placeholder(truth[[r]]), logical(nrow(truth))),
  1,
  \(x) if (any(x)) which(x)[1] else NA_integer_
)
data.frame(
  rank = ranks,
  placeholder = vapply(ranks, \(r) sum(placeholder(truth[[r]])), integer(1)),
  real_below_placeholder = vapply(
    seq_along(ranks),
    \(i) {
      sum(
        !is.na(first_ph) &
          first_ph < i &
          !is.na(truth[[ranks[i]]]) &
          !placeholder(truth[[ranks[i]]])
      )
    },
    integer(1)
  ),
  no_truth = vapply(ranks, \(r) sum(is.na(truth[[r]])), integer(1))
) |>
  print()
message("Controls with at least one placeholder: ", sum(!is.na(first_ph)))

truth_B <- truth
truth_D <- truth
for (i in seq_along(ranks)) {
  r <- ranks[i]
  truth_B[[r]][!is.na(first_ph) & first_ph <= i] <- NA_character_
  truth_D[[r]][placeholder(truth[[r]])] <- NA_character_
}

# Assignment column of each rank: Species is scored on Gen_sp (chapter 01).
col_rank <- setNames(c(ranks[-length(ranks)], "Gen_sp"), ranks)
defaults <- values_map_for("pauvert") |>
  filter(.data$input == .env$input) |>
  filter(
    (method %in% c("dada2", "sintax") & min_bootstrap == 0.5) |
      (method == "lca" & lca_cutoff == 0.51) |
      (method == "blastn" & vote_algorithm == "rel_majority" & min_id == 95)
  ) |>
  transmute(full_name, algo = method, db = assignment_db(full_name)) |>
  inner_join(filter(db_meta, simplification == "full"), by = "db")
stopifnot(nrow(defaults) == 12)

rates <- lapply(seq_len(nrow(defaults)), function(k) {
  lapply(ranks, function(r) {
    v <- tt[ext, paste0(col_rank[[r]], "_", defaults$full_name[k])]
    hit <- function(tr) !is.na(v) & !is.na(tr[[r]]) & v == tr[[r]]
    tibble(
      full_name = defaults$full_name[k],
      algo = defaults$algo[k],
      db = defaults$db[k],
      rank = r,
      A = sum(hit(truth)) / nrow(truth),
      B = sum(hit(truth_B)) / nrow(truth),
      C = sum(hit(truth_B)) / sum(!is.na(truth_B[[r]])),
      D = sum(hit(truth_D)) / sum(!is.na(truth_D[[r]])),
      scorable_C = sum(!is.na(truth_B[[r]])),
      scorable_D = sum(!is.na(truth_D[[r]]))
    )
  }) |>
    bind_rows()
}) |>
  bind_rows()

# Variant D is chapter 01's ext_correct since decision D (2026-09-30).
ch01 <- readRDS("analysis/_cache/res_comp_tax.rds") |>
  keep_truth_rule() |>
  filter(
    dataset == "pauvert",
    unit == "ASV",
    metrics == "ext_correct",
    method_db %in% defaults$full_name
  ) |>
  mutate(
    rank = ifelse(tax_level %in% c("Gen_sp", "Species"), "Species", tax_level)
  ) |>
  select(full_name = method_db, rank, ch01 = values)
check <- inner_join(rates, ch01, by = c("full_name", "rank"))
stopifnot(nrow(check) == nrow(rates), all(abs(check$D - check$ch01) < 1e-9))
message(
  "\nVariant D reproduces chapter 01's ext_correct (",
  nrow(check),
  " values): OK"
)

message(
  "\next_correct, full databases, mean over the three (min-max), per method and rank"
)
rates |>
  group_by(algo, rank = factor(rank, levels = ranks)) |>
  summarise(
    across(c(A, B, C, D), \(v) {
      sprintf("%.2f (%.2f-%.2f)", mean(v), min(v), max(v))
    }),
    n_C = first(scorable_C),
    n_D = first(scorable_D),
    .groups = "drop"
  ) |>
  filter(rank != "Kingdom") |>
  arrange(rank, algo) |>
  as.data.frame() |>
  print()

message("\nLargest change per variant against A (any method, database, rank)")
rates |>
  summarise(across(c(B, C, D), \(v) round(max(abs(v - A)), 3))) |>
  print()

message("\nPlaceholders in the mock truth (Pauvert), for the one-rule claim")
cfg <- mock_sanger_tables[mock_sanger_tables$dataset == "pauvert", ]
mock <- build_mock_truth(tar_read(d_asv, store = here("store_dada2")), cfg)
print(vapply(ranks, \(r) sum(placeholder(mock[[r]])), integer(1)))
