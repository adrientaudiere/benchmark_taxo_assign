# Worked examples for cross_val(). Not sourced by any pipeline. Run pieces
# interactively after `source(here::here("R/cross_val.R"))`.
#
# Each CV_* run can take minutes to hours depending on `fold_number` and the
# database size — start with a high `fold_number` and a low `fold_tested` to
# keep each fold small.

CV_unite <- cross_val(
  ref_fasta = "data/data_raw/refseq/sintax_format/Unite.fasta",
  method = "sintax",
  verbose = TRUE,
  fold_number = 1000,
  fold_tested = 5,
  min_bootstrap = 0.5,
  compute_by_tax_level = TRUE,
  nproc = 4
)

CV_eukF <- cross_val(
  ref_fasta = "data/data_raw/refseq/sintax_format/EUK_ITS_v1_9_3_Fungi.fasta",
  method = "sintax",
  verbose = TRUE,
  fold_number = 1000,
  fold_tested = 5,
  min_bootstrap = 0.5,
  compute_by_tax_level = TRUE,
  nproc = 4
)

CV_unite_minboot <- cross_val(
  ref_fasta = "data/data_raw/refseq/sintax_format/Unite.fasta",
  method = "sintax",
  verbose = TRUE,
  fold_number = 10000,
  fold_tested = 2,
  min_bootstrap = c(0.5, 0.6),
  compute_by_tax_level = TRUE,
  nproc = 4
)

CV_unite_blast <- cross_val(
  ref_fasta = "data/data_raw/refseq/sintax_format/Unite.fasta",
  method = "blastn",
  verbose = TRUE,
  fold_number = 500,
  fold_tested = 2,
  compute_by_tax_level = TRUE,
  nproc = 4
)

# Aggregate metrics by rank.
CV_unite$metrics |>
  ggplot() +
  geom_bar(aes(y = as.numeric(mean),
               x = factor(name, levels = rev(c("Kingdom", "Phylum", "Order",
                                               "Class", "Family", "Genus",
                                               "Species")))),
           stat = "identity") +
  geom_errorbar(aes(x = name, ymin = mean - sd, ymax = mean + sd),
                width = .2, position = position_dodge(.9)) +
  coord_flip() +
  facet_grid(~ metric, scales = "free")

# Class-level metrics (excluding Incertae_sedis).
CV_unite$metrics_by_taxonomic_level |>
  filter(taxonomic_rank == "Class") |>
  filter(!grepl("*_Incertae_sedis", taxonomic_value, ignore.case = TRUE)) |>
  ggplot() +
  geom_bar(aes(y = as.numeric(mean), x = taxonomic_value), stat = "identity") +
  geom_errorbar(aes(x = taxonomic_value, ymin = mean - sd, ymax = mean + sd),
                width = .2, position = position_dodge(.9)) +
  coord_flip() +
  facet_grid(~ metric, scales = "free")

# Taxa with poor F1 (database/method weak spots).
CV_unite$metrics_by_taxonomic_level |>
  filter(metric == "F1_score") |>
  filter(!grepl("*_Incertae_sedis", taxonomic_value)) |>
  filter(mean > 0, mean < 0.7) |>
  filter(bootstrap == 0.5) |>
  filter(taxonomic_rank %in% c("Phylum", "Order", "Class", "Family", "Genus")) |>
  ggplot() +
  geom_bar(aes(y = as.numeric(mean), x = taxonomic_value, fill = taxonomic_rank),
           stat = "identity") +
  geom_errorbar(aes(x = taxonomic_value, ymin = mean - sd, ymax = mean + sd),
                width = .2, position = position_dodge(.9)) +
  coord_flip() +
  facet_grid(~ taxonomic_rank)

# Taxa with high false-negative rate.
CV_unite$metrics_by_taxonomic_level |>
  filter(metric == "FN") |>
  filter(!grepl("*_Incertae_sedis", taxonomic_value)) |>
  filter(mean > 0, mean < 1) |>
  filter(bootstrap == 0.5) |>
  filter(taxonomic_rank %in% c("Phylum", "Order", "Class", "Family", "Genus")) |>
  ggplot() +
  geom_bar(aes(y = as.numeric(mean), x = taxonomic_value, fill = taxonomic_rank),
           stat = "identity") +
  geom_errorbar(aes(x = taxonomic_value, ymin = mean - sd, ymax = mean + sd),
                width = .2, position = position_dodge(.9)) +
  coord_flip() +
  facet_grid(~ taxonomic_rank)

# Taxa with high false-positive rate.
CV_unite$metrics_by_taxonomic_level |>
  filter(metric == "FP") |>
  filter(!grepl("*_Incertae_sedis", taxonomic_value)) |>
  filter(mean > 0, mean < 1) |>
  filter(bootstrap == 0.5) |>
  filter(taxonomic_rank %in% c("Phylum", "Order", "Class", "Family", "Genus")) |>
  ggplot() +
  geom_bar(aes(y = as.numeric(mean), x = taxonomic_value, fill = taxonomic_rank),
           stat = "identity") +
  geom_errorbar(aes(x = taxonomic_value, ymin = mean - sd, ymax = mean + sd),
                width = .2, position = position_dodge(.9)) +
  coord_flip() +
  facet_grid(~ taxonomic_rank)

# Effect of bootstrap threshold (0.5 vs 0.6) on per-taxon metrics.
CV_unite_minboot$metrics_by_taxonomic_level |>
  filter(taxonomic_rank %in% c("Order", "Class", "Family")) |>
  filter(!grepl("*_Incertae_sedis", taxonomic_value)) |>
  filter(!is.na(sd)) |>
  pivot_wider(names_prefix = "mean_",
              names_from = bootstrap,
              values_from = mean) |>
  mutate(diff = mean_0.5 - mean_0.6) |>
  filter(diff != 0) |>
  ggplot() +
  geom_bar(aes(y = as.numeric(diff), x = taxonomic_value),
           stat = "identity", position = position_dodge()) +
  coord_flip() +
  facet_grid(taxonomic_rank ~ metric, scales = "free")
