suppressMessages({
  library(targets)
  library(dplyr)
  library(phyloseq)
})
source("config.R")

# S1.3: does blastn assign the mock ASVs at min_cover = 80 on the production
# nine-database grid? Rates are computed on the 195 mock ASVs alone, the
# `fake_` and `external_` controls reported apart (C22-P8).

d <- tar_read(d_all_taxo, store = "store_assign_taxo")
tt <- as.data.frame(tax_table(d), stringsAsFactors = FALSE)
unit <- ifelse(
  grepl("^external_", rownames(tt)),
  "external",
  ifelse(grepl("^fake_", rownames(tt)), "fake", "mock")
)

rate_for <- function(cols, keep) {
  vapply(cols, \(x) mean(!is.na(tt[[x]][keep])), numeric(1))
}

cols <- colnames(tt)[grepl("__", colnames(tt))]
info <- tibble(
  col = cols,
  rank = sub("_.*$", "", cols),
  method = sub("__.*$", "", sub("^[^_]+_", "", cols)),
  db = sub("___.*$", "", sub("^[^_]+_[a-z0-9_]+?__", "", cols)),
  mock = rate_for(cols, unit == "mock"),
  fake = rate_for(cols, unit == "fake"),
  external = rate_for(cols, unit == "external")
)

message(
  "blastn, min_cover = ",
  blastn_min_cover,
  ", rate assigned on the ",
  sum(unit == "mock"),
  " mock ASVs"
)
info |>
  filter(method == "blastn", rank %in% c("Kingdom", "Genus", "Species")) |>
  group_by(rank, db) |>
  summarise(n = n(), min = min(mock), max = max(mock), .groups = "drop") |>
  arrange(match(rank, c("Kingdom", "Genus", "Species")), db) |>
  as.data.frame() |>
  print()

message("\nblastn at Genus, by vote and min_id (mock ASVs)")
info |>
  filter(method == "blastn", rank == "Genus") |>
  mutate(
    vote = sub("^.*___0\\.5\\.\\.\\.([a-z_]+)\\.\\.\\..*$", "\\1", col),
    min_id = sub("^.*\\.\\.\\.", "", col)
  ) |>
  group_by(vote, min_id) |>
  summarise(mean = mean(mock), min = min(mock), .groups = "drop") |>
  as.data.frame() |>
  print()

message("\nall methods at Genus: mock ASVs against the two control classes")
info |>
  filter(rank == "Genus") |>
  group_by(method) |>
  summarise(
    n_col = n(),
    mock = mean(mock),
    fake = mean(fake),
    external = mean(external),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()

n_empty <- sum(info$method == "blastn" & info$mock == 0)
message(
  "\nblastn columns assigning zero mock ASV: ",
  n_empty,
  " of ",
  sum(info$method == "blastn")
)
stopifnot(n_empty == 0)
