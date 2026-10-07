suppressMessages({
  library(targets)
  library(dplyr)
  library(tidyr)
  library(phyloseq)
  library(MiscMetabar)
  library(comparpq)
})
source("config.R")
source("R/assignment_helpers.R")

# Evidence for ROADMAP SA1 (2026-09-30), measured again on the ITSx-only stores
# (OD 34). Former measure (2026-09-23, before OD 34, raw inputs):
# analysis/sandbox/check_controls_by_simplification.R, text in HISTORY.md §6.4.
#
# Question: on the three `Fungi` databases at Genus, where naming a control is
# a false positive by construction, which methods abstain on the controls, and
# does the abstention cost them accuracy on the real units (PPV / F1)?
#
# Part A, Pauvert: read off chapter 01's score table
# (analysis/_cache/res_comp_tax.rds, BENCHMARK_ANALYSIS_MINI=FALSE render), so
# the labels are cleaned exactly as the scoring cleans them. On the `Fungi`
# databases the scoring is "matrix": every control given a genus is an FP.
# PPV_real and F1_real leave the controls out (C22-P8).
# Part B, wood: no truth, so the tax_table itself, cleaned with clean_taxo();
# the same parser run on Pauvert must reproduce part A's control counts.
#
# Two views per method: the default setting (dada2 / sintax bootstrap 0.5, lca
# cutoff 0.51, blastn rel_majority at 95 % identity), and the mean over every
# sweep value, the view of the 2026-09-23 measure.

is_default <- function(algo, method_db) {
  case_when(
    algo %in% c("dada2", "sintax") ~ grepl("___0\\.5$", method_db),
    algo == "lca" ~ grepl("___0\\.51$", method_db),
    algo == "blastn" ~ grepl(
      "\\.\\.\\.rel_majority\\.\\.\\.100\\.\\.\\.95$",
      method_db
    ),
    .default = FALSE
  )
}

# Part A ----------------------------------------------------------------------

store_pauvert <- "store_assign_taxo"
n_ctrl <- function(target) {
  units <- taxa_names(tar_read_raw(target, store = store_pauvert))
  c(fake = sum(grepl("^fake_", units)), ext = sum(grepl("^external_", units)))
}
n_ctrl_unit <- list(
  ASV = n_ctrl("d_all_taxo_ITS1_asv"),
  OTU = n_ctrl("d_all_taxo_ITS1_otu")
)

res <- readRDS("analysis/_cache/res_comp_tax.rds") |>
  filter(
    dataset == "pauvert",
    tax_level == "Genus",
    external_scoring == "matrix",
    !is_consensus,
    unit %in% c("ASV", "OTU")
  ) |>
  select(method_db, unit, algo, db, metrics, values) |>
  pivot_wider(names_from = metrics, values_from = values) |>
  mutate(
    default = is_default(algo, method_db),
    n_fake = vapply(unit, \(u) n_ctrl_unit[[u]][["fake"]], numeric(1)),
    n_ext = vapply(unit, \(u) n_ctrl_unit[[u]][["ext"]], numeric(1)),
    real = 1 - NA_real,
    external = ctrl_assigned_ext / n_ext,
    fake = ctrl_assigned_fake / n_fake,
    FP_real = misassign_seq * n_real,
    PPV_real = TP / (TP + FP_real),
    F1_real = 2 * TP / (2 * TP + FP_real + FN)
  )

stopifnot(all(abs(res$fake - (1 - res$NA_fake)) < 1e-9))
stopifnot(all(res$db %in% db_meta$db[db_meta$simplification == "Fungi"]))

summarise_view <- function(x) {
  x |>
    group_by(unit, algo) |>
    summarise(
      n_col = n(),
      across(
        c(real, external, fake, PPV, F1_score, PPV_real, F1_real, TPR),
        \(v) round(mean(v), 3)
      ),
      .groups = "drop"
    ) |>
    arrange(unit, external) |>
    as.data.frame()
}

message("\nSA1, Pauvert ITS1, Fungi databases, Genus: default setting")
print(summarise_view(filter(res, default)))

message(
  "\nSA1, Pauvert ITS1, Fungi databases, Genus: mean over every sweep value"
)
print(summarise_view(res))

message("\nDefault setting, per database (ASV)")
res |>
  filter(default, unit == "ASV") |>
  transmute(
    algo,
    db,
    across(c(real, external, fake, PPV, PPV_real, F1_score, F1_real), \(v) {
      round(v, 3)
    })
  ) |>
  arrange(algo, db) |>
  as.data.frame() |>
  print()

# Part B ----------------------------------------------------------------------

control_rates <- function(target, store) {
  d <- tar_read_raw(target, store = store) |> clean_taxo()
  tt <- as.data.frame(tax_table(d), stringsAsFactors = FALSE)
  class <- ifelse(
    grepl("^external_", rownames(tt)),
    "external",
    ifelse(grepl("^fake_", rownames(tt)), "fake", "real")
  )
  cols <- grep("^Genus_.*__", colnames(tt), value = TRUE)
  tibble(
    col = cols,
    method_db = sub("^Genus_", "", cols),
    algo = sub("__.*$", "", sub("^Genus_[^_]+_[a-z]+_", "", cols)),
    db = sub("___.*$", "", sub("^.*?__", "", cols)),
    real = vapply(
      cols,
      \(x) mean(!is.na(tt[[x]][class == "real"])),
      numeric(1)
    ),
    external = vapply(
      cols,
      \(x) mean(!is.na(tt[[x]][class == "external"])),
      numeric(1)
    ),
    fake = vapply(cols, \(x) mean(!is.na(tt[[x]][class == "fake"])), numeric(1))
  ) |>
    filter(db %in% db_meta$db[db_meta$simplification == "Fungi"]) |>
    mutate(default = is_default(algo, method_db))
}

# The parser against the score table, on Pauvert ASV.
check <- control_rates("d_all_taxo_ITS1_asv", store_pauvert) |>
  inner_join(
    filter(res, unit == "ASV") |>
      select(method_db, ext_a = external, fake_a = fake),
    by = "method_db"
  )
stopifnot(nrow(check) == sum(res$unit == "ASV"))
stopifnot(all(abs(check$external - check$ext_a) < 1e-9))
stopifnot(all(abs(check$fake - check$fake_a) < 1e-9))
message(
  "\nParser reproduces the score table on Pauvert ASV (",
  nrow(check),
  " columns): OK"
)

wood <- control_rates("d_all_taxo_ITS2_asv", "store_assign_taxo_korhonen_wood")
message(
  "\nSA1, wood ITS2 ASV (no truth: real = assignment rate), Fungi databases, Genus"
)
for (view in c("default setting", "mean over every sweep value")) {
  message("-- ", view)
  wood |>
    filter(if (view == "default setting") default else TRUE) |>
    group_by(algo) |>
    summarise(
      n_col = n(),
      across(c(real, external, fake), \(v) round(mean(v), 3)),
      .groups = "drop"
    ) |>
    arrange(external) |>
    as.data.frame() |>
    print()
}
