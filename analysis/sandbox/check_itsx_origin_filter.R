suppressMessages({
  library(targets)
  library(dplyr)
  library(phyloseq)
  library(MiscMetabar)
  library(comparpq)
})
source("config.R")
source("R/assignment_helpers.R")
source("R/mock_truth.R")

# Evidence for ROADMAP SA7 (2026-09-30), on the ITSx-only stores (OD 34).
#
# Question: `itsx_pq(add_origin = TRUE)` writes the putative origin ITSx gives
# each unit (`ITSx_origin`, the name of the profile set that found it: "Fungi",
# "Tracheophyta", ...; NA when ITSx found nothing). Would a cheap filter
# keeping only the units ITSx calls "Fungi" already remove the non-fungal
# false positives, i.e. the `external_` controls a method calls Fungi at
# Kingdom? And what does the same filter cost on the real units?
#
# Kingdom columns at the default setting (dada2 / sintax bootstrap 0.5, lca
# cutoff 0.51, blastn rel_majority at 95 % identity), labels cleaned with
# clean_taxo() as the scoring cleans them. Calling an external control Fungi is
# an error on every database (they are all non-fungal), so the nine databases
# are read, per simplification. Pauvert ITS1 ASV (real units checked against
# the mock truth) and wood ITS2 ASV (no truth).

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

read_input <- function(target, store) {
  d <- tar_read_raw(target, store = store) |> clean_taxo()
  tt <- as.data.frame(tax_table(d), stringsAsFactors = FALSE)
  tt$class <- ifelse(
    grepl("^external_", rownames(tt)),
    "external",
    ifelse(grepl("^fake_", rownames(tt)), "fake", "real")
  )
  tt$itsx_fungi <- !is.na(tt$ITSx_origin) & tt$ITSx_origin == "Fungi"
  tt
}

# One row per method x database: share of the external controls called Fungi
# at Kingdom, and what is left of them once the ITSx filter is applied.
kingdom_calls <- function(tt) {
  cols <- grep("^Kingdom_.*__", colnames(tt), value = TRUE)
  ext <- tt$class == "external"
  tibble(
    method_db = sub("^Kingdom_", "", cols),
    algo = sub("__.*$", "", sub("^Kingdom_[^_]+_[a-z]+_", "", cols)),
    db = sub("___.*$", "", sub("^.*?__", "", cols)),
    ext_fungi = vapply(cols, \(x) sum(ext & tt[[x]] %in% "Fungi"), integer(1)),
    ext_fungi_after_filter = vapply(
      cols,
      \(x) sum(ext & tt[[x]] %in% "Fungi" & tt$itsx_fungi),
      integer(1)
    )
  ) |>
    filter(is_default(algo, method_db)) |>
    left_join(db_meta, by = "db") |>
    select(algo, simplification, db, ext_fungi, ext_fungi_after_filter) |>
    arrange(simplification, algo, db)
}

origin_table <- function(tt, label) {
  message("\n", label, ": ITSx origin per unit class")
  print(table(
    class = tt$class,
    origin = ifelse(
      is.na(tt$ITSx_origin),
      "none (ITSx found nothing)",
      ifelse(tt$itsx_fungi, "Fungi", "non-fungal profile")
    )
  ))
}

# Pauvert ----------------------------------------------------------------------

pauvert <- read_input("d_all_taxo_ITS1_asv", "store_assign_taxo")
origin_table(pauvert, "Pauvert ITS1 ASV")

ext_truth <- external_control_truth(rownames(pauvert)[
  pauvert$class == "external"
])
itsx_f <- ext_truth$taxon %in% rownames(pauvert)[pauvert$itsx_fungi]
message(
  "\nTrue kingdom of the external controls ITSx calls Fungi (",
  sum(itsx_f),
  "):"
)
print(table(ext_truth$Kingdom[itsx_f], useNA = "ifany"))
stopifnot(!any(ext_truth$Kingdom %in% "Fungi"))

message(
  "\nExternal controls (n = ",
  sum(pauvert$class == "external"),
  ") called Fungi at Kingdom, before and after the ITSx filter"
)
kingdom_calls(pauvert) |> as.data.frame() |> print()

# Cost on the real units: the mock is fungal, so every real unit ITSx does not
# call Fungi is a fungal unit the filter would drop.
cfg <- mock_sanger_tables[mock_sanger_tables$dataset == "pauvert", ]
truth <- build_mock_truth(tar_read(d_asv, store = "store_dada2"), cfg)
real <- pauvert[pauvert$class == "real", ]
real_truth <- truth[match(rownames(real), truth$taxon), ]
message("\nReal Pauvert ASVs, ITSx origin against the truth kingdom")
print(table(
  truth_kingdom = real_truth$Kingdom,
  origin = ifelse(is.na(real$ITSx_origin), "none", real$ITSx_origin),
  useNA = "ifany"
))

dropped <- !real$itsx_fungi
message("\nThe ", sum(dropped), " real ASVs the filter drops, by truth genus:")
print(sort(
  table(real_truth$Genus[dropped], useNA = "ifany"),
  decreasing = TRUE
))

# Correct genus lost per method: the true positives among the dropped units.
genus_cols <- grep("^Genus_.*__", colnames(real), value = TRUE)
genus_cols <- genus_cols[is_default(
  sub("__.*$", "", sub("^Genus_[^_]+_[a-z]+_", "", genus_cols)),
  sub("^Genus_", "", genus_cols)
)]
correct <- function(x) {
  v <- real[[x]]
  !is.na(v) &
    !is.na(real_truth$Genus) &
    (v == real_truth$Genus |
      v %in% real_truth$Genus_accepted & v == real_truth$Genus_accepted)
}
tp_lost <- tibble(
  algo = sub("__.*$", "", sub("^Genus_[^_]+_[a-z]+_", "", genus_cols)),
  db = sub("___.*$", "", sub("^.*?__", "", genus_cols)),
  tp = vapply(genus_cols, \(x) sum(correct(x)), integer(1)),
  tp_lost = vapply(genus_cols, \(x) sum(correct(x) & dropped), integer(1))
) |>
  left_join(db_meta, by = "db") |>
  group_by(algo, simplification) |>
  summarise(tp = mean(tp), tp_lost = mean(tp_lost), .groups = "drop") |>
  mutate(share_lost = round(tp_lost / tp, 3))
message(
  "\nCorrect genera (mean over the three databases) the filter would remove"
)
print(as.data.frame(tp_lost))

# Wood -------------------------------------------------------------------------

wood <- read_input("d_all_taxo_ITS2_asv", "store_assign_taxo_korhonen_wood")
origin_table(wood, "Wood ITS2 ASV")
message(
  "\nExternal controls (n = ",
  sum(wood$class == "external"),
  ") called Fungi at Kingdom, before and after the ITSx filter"
)
kingdom_calls(wood) |> as.data.frame() |> print()

# What the methods call the real wood units ITSx does not call Fungi.
k_full <- grep(
  "^Kingdom_.*__(Unite_s_all_20250219|EUK_ITS_v2.1)___",
  colnames(wood),
  value = TRUE
)
k_full <- k_full[is_default(
  sub("__.*$", "", sub("^Kingdom_[^_]+_[a-z]+_", "", k_full)),
  sub("^Kingdom_", "", k_full)
)]
wood_nf <- wood[wood$class == "real" & !wood$itsx_fungi, ]
message(
  "\nWood real units ITSx does not call Fungi (",
  nrow(wood_nf),
  "): Kingdom given on the two full databases"
)
for (x in k_full) {
  message("-- ", sub("^Kingdom_", "", x))
  print(table(wood_nf[[x]], useNA = "ifany"))
}
