# How many OTUs of each biological dataset survive an abundance x occurrence
# filter (developer, 2026-09-26)? The benchmark should measure assignment on
# the common OTUs rather than on a random draw, so the draw of OD 30b is to be
# replaced by one filter applied to every dataset: `n_seq > X & n_occur >= Y`,
# where n_seq is the total read count of an OTU and n_occur the number of
# samples it occurs in. This script only counts; the developer picks X and Y.
#
# Wood is presence / absence (Korhonen et al. 2024 publish no read counts), so
# only its occurrence threshold applies: its n_seq column is NA.
#
# Run from the project root (reads the 155 MB GSSP table, about a minute):
#   Rscript analysis/sandbox/check_bio_otu_filters.R

suppressMessages(library(dplyr))

# Air: every published OTU, reads and occurrences from the full table.
air_counts <- data.table::fread(
  "data/data_raw/bio_gssp_air/otu.table.csv",
  sep = ";",
  data.table = FALSE
)
air_mat <- as.matrix(air_counts[, -1])
air <- tibble(
  otu = colnames(air_mat),
  n_seq = colSums(air_mat),
  n_occur = colSums(air_mat > 0)
)

# Wood: both deposits pooled on the hash, as load_korhonen_wood() does, before
# its "present in at least one sample" rule.
source("R/bio_loaders.R")
wood_pq <- load_korhonen_wood(
  "data/data_raw/bio_korhonen_wood",
  min_occur = 1
)
wood_mat <- as(phyloseq::otu_table(wood_pq), "matrix")
wood <- tibble(
  otu = colnames(wood_mat),
  n_seq = NA_real_,
  n_occur = colSums(wood_mat > 0)
)

seq_thresholds <- c(0, 10, 100, 1000, 10000)
occur_thresholds <- c(1, 2, 5, 10, 20, 50)

grid_count <- function(tab, name) {
  n_samples <- if (name == "air") nrow(air_mat) else nrow(wood_mat)
  tidyr::expand_grid(
    n_seq_gt = seq_thresholds,
    n_occur_ge = occur_thresholds
  ) |>
    rowwise() |>
    mutate(
      kept = if (all(is.na(tab$n_seq))) {
        if (n_seq_gt == 0) {
          sum(tab$n_occur >= n_occur_ge)
        } else {
          NA_integer_
        }
      } else {
        sum(tab$n_seq > n_seq_gt & tab$n_occur >= n_occur_ge)
      }
    ) |>
    ungroup() |>
    mutate(dataset = name, n_otu = nrow(tab), n_samples = n_samples)
}

res <- bind_rows(grid_count(air, "air"), grid_count(wood, "wood"))
for (d in c("air", "wood")) {
  sub <- res[res$dataset == d, ]
  cat(sprintf(
    "\n== %s: %d OTUs, %d samples. OTUs kept (rows: n_seq >, columns: n_occur >=)\n",
    d,
    sub$n_otu[1],
    sub$n_samples[1]
  ))
  wide <- tidyr::pivot_wider(
    sub[, c("n_seq_gt", "n_occur_ge", "kept")],
    names_from = n_occur_ge,
    values_from = kept
  )
  print(as.data.frame(wide), row.names = FALSE)
}
cat("\nOccurrence quantiles (air, wood):\n")
print(rbind(
  air = quantile(air$n_occur, c(0.5, 0.75, 0.9, 0.99)),
  wood = quantile(wood$n_occur, c(0.5, 0.75, 0.9, 0.99))
))
utils::write.csv(
  res,
  "data/data_final/bio_otu_filter_grid.csv",
  row.names = FALSE
)
