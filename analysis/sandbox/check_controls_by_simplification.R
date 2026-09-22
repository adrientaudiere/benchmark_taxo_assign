suppressMessages({
  library(targets)
  library(dplyr)
  library(phyloseq)
})
source("config.R")

# Evidence for ROADMAP section 7 (SA1-SA3): rate of units given a genus, per
# method and per simplification, on the three unit classes separately.
#
# Never pool the simplifications. Naming an `external_` control is an error on
# a `Fungi` database only; on `full` it is correct, on `Fungi+rep` correct at
# Kingdom only (decisions 6, 15, 22). Pooling also mixes the `itsx_*` columns,
# which exist for the six non-full databases only, with the raw ones, which
# exist for all nine: a pooled mean then compares different database sets.

d <- tar_read(d_all_taxo, store = "store_assign_taxo")
tt <- as.data.frame(tax_table(d), stringsAsFactors = FALSE)
unit <- ifelse(
  grepl("^external_", rownames(tt)),
  "external",
  ifelse(grepl("^fake_", rownames(tt)), "fake", "mock")
)

rate <- function(cols, keep) {
  vapply(cols, \(x) mean(!is.na(tt[[x]][keep])), numeric(1))
}

cols <- colnames(tt)[grepl("__", colnames(tt))]
info <- tibble(
  col = cols,
  rank = sub("_.*$", "", cols),
  method = sub("__.*$", "", sub("^[^_]+_", "", cols)),
  db = sub("___.*$", "", sub("^[^_]+_[a-z0-9_]+?__", "", cols)),
  mock = rate(cols, unit == "mock"),
  fake = rate(cols, unit == "fake"),
  external = rate(cols, unit == "external")
) |>
  left_join(db_meta, by = "db")

genus <- filter(info, rank == "Genus")

message(
  "SA1: Fungi databases at Genus, where naming an external control is an error"
)
genus |>
  filter(simplification == "Fungi") |>
  group_by(method) |>
  summarise(
    n = n(),
    mock = round(mean(mock), 3),
    external = round(mean(external), 3),
    fake = round(mean(fake), 3),
    .groups = "drop"
  ) |>
  arrange(external) |>
  as.data.frame() |>
  print()

message(
  "\nSA2 / SA3: each ITSx method against its raw twin, non-full databases only"
)
pairs <- c("blastn", "dada2", "lca", "sintax")
genus |>
  filter(
    simplification != "full",
    method %in% c(pairs, paste0("itsx_", pairs))
  ) |>
  mutate(
    family = sub("^itsx_", "", method),
    input = ifelse(grepl("^itsx_", method), "itsx", "raw")
  ) |>
  group_by(family, input) |>
  summarise(
    n = n(),
    mock = round(mean(mock), 3),
    external = round(mean(external), 3),
    fake = round(mean(fake), 3),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()

# The pooling trap of 2026-09-22: the same comparison unmatched, kept as a
# warning. blastn averages over 9 databases and itsx_blastn over 6, so the
# `full` rows (where naming an external control is correct) inflate blastn
# alone and manufacture a gap that matching removes.
message(
  "\nThe same two columns pooled over every database (wrong, kept as a warning)"
)
genus |>
  filter(method %in% c("blastn", "itsx_blastn")) |>
  group_by(method) |>
  summarise(n = n(), external = round(mean(external), 3), .groups = "drop") |>
  as.data.frame() |>
  print()

matched <- genus |>
  filter(method %in% c("blastn", "itsx_blastn"), simplification != "full") |>
  group_by(method) |>
  summarise(external = mean(external), .groups = "drop")
stopifnot(abs(diff(matched$external)) < 1e-9)
message("\nmatched, blastn and itsx_blastn agree on the external controls: OK")
