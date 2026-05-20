# Thin script: reconstructs `values_map` from script_assign_taxo_parallel.R
# so that benchmark.qmd has the column-name machinery available.
# Do NOT source this for tar_make — it is qmd-only.

library("conflicted")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")

here::i_am("analysis/benchmark.qmd")
source(here("config.R"))

methods <- tidyr::expand_grid(
  method = c("dada2", "sintax", "lca"),
  min_bootstrap = c(0.4, 0.5, 0.6)
) |>
  full_join(tidyr::expand_grid(
    method = c("blastn"),
    vote_algorithm = c("rel_majority", "abs_majority", "unanimity"),
    nb_voting = 100,
    min_bootstrap = 0.5
  ))

values_map <-
  tidyr::expand_grid(
    methods,
    db = c(
      "Unite",
      "Unite_Fungi",
      "EUK_ITS_v1_9_3",
      "EUK_ITS_v1_9_3_Fungi",
      "EUK_ITS_v1_9_3_Fungi_cut",
      "EUK_SSU_v1_9_3",
      "EUK_SSU_v1_9_3_Fungi",
      "EUK_SSU_v1_9_3_cut",
      "EUK_SSU_v1_9_3_Fungi_cut"
    )
  ) |>
  mutate(cutadapted_db = ifelse(grepl("cut", db), "cut", "")) |>
  mutate(db_filter = ifelse(grepl("Fungi", db), "Fungi", "")) |>
  mutate(db_name = db) |>
  mutate(do_clean_pq = ifelse(method == "dada2", TRUE, FALSE)) |>
  mutate(db_path = paste0(ifelse(
    method == "dada2",
    paste0("data/data_raw/refseq/dada2_format/", db_name),
    paste0("data/data_raw/refseq/sintax_format/", db_name)
  ), ".fasta")) |>
  mutate(full_name = gsub(
    "...NA...NA", "",
    paste0(method, "__", db, "___",
           min_bootstrap, "...",
           vote_algorithm, "...",
           nb_voting)
  ))