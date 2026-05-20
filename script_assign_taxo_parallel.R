library("conflicted")
library("MiscMetabar")
devtools::load_all("/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/MiscMetabar/")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
library("autometric")

if (tar_active()) {
  log_start(path = "data/data_final/autometric_log_assign_taxo.txt", seconds = 1)
}

here::i_am("script_assign_taxo_parallel.R")
source(here("config.R"))
source(here("R/combine_taxo_assignments.R"))
source("/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/comparpq/R/compare_taxo.R")
source("/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/comparpq/R/fake_creation.R")

tar_option_set(
  seed = targets_seed,
  controller = crew::crew_controller_group(
    crew::crew_controller_local(name = "dada2_ctrl",  workers = 1,         seconds_idle = 60),
    crew::crew_controller_local(name = "fast_ctrl",   workers = n_workers, seconds_idle = 60)
  )
)

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
      "EUK_ITS_v1_9_3_Fungi",       # Q2.3: files exist in both formats
      "EUK_ITS_v1_9_3_Fungi_cut",   # Q2.3
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
  mutate(controller = ifelse(method == "dada2", "dada2_ctrl", "fast_ctrl")) |>
  mutate(full_name = gsub(
    "...NA...NA", "",
    paste0(method, "__", db, "___",
           min_bootstrap, "...",
           vote_algorithm, "...",
           nb_voting)
  ))

assignment_targets <- tarchetypes::tar_eval(
  tar_target(
    full_name,
    {
      # Tag the autometric log with this target name so benchmark_costs can
      # split runtime/memory per (method, db) row.
      autometric::log_phase_set(full_name)
      add_new_taxonomy_pq(
        d_asv_for_assignation,
        method = method,
        ref_fasta = db_path,
        suffix = paste0("_", full_name),
        min_bootstrap = min_bootstrap,
        vote_algorithm = vote_algorithm,
        nb_voting = nb_voting
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = controller)
    )
  ),
  values_map
)

tar_plan(
  tar_target(d_asv, tar_read(d_asv, store = here::here("store_dada2"))),
  tar_target(
    file_taxo_mock,
    here(taxo_mock_csv),
    format = "file"
  ),
  tar_target(
    taxo_mock,
    read.csv(file_taxo_mock) |>
      select(any_of(c(
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
      ))) |>
      magrittr::set_rownames(read.csv(file_taxo_mock)$MockStrain)
  ),
  tar_target(
    d_asv_for_assignation_fake,
    add_shuffle_seq_pq(d_asv, prop_fake = prop_fake)
  ),
  tar_target(
    d_asv_for_assignation,
    add_external_seq_pq(
      d_asv_for_assignation_fake,
      Biostrings::readDNAStringSet(here(fake_ref_fasta))
    )
  ),
  assignment_targets,
  tarchetypes::tar_combine(
    d_all_taxo,
    assignment_targets,
    command = combine_taxo_assignments(d_asv_for_assignation, !!!.x)
  ),
  tar_target(
    benchmark_costs,
    {
      # Force aggregation to run after the combine, so the log is complete.
      d_all_taxo
      log_df <- autometric::log_read(
        here("data/data_final/autometric_log_assign_taxo.txt")
      )
      log_df |>
        dplyr::filter(!is.na(phase), phase != "") |>
        dplyr::group_by(phase) |>
        dplyr::summarise(
          wall_time_s      = as.numeric(max(time) - min(time)),
          peak_resident_mb = max(resident, na.rm = TRUE),
          mean_cpu_pct     = mean(cpu, na.rm = TRUE),
          n_samples        = dplyr::n(),
          .groups = "drop"
        ) |>
        dplyr::inner_join(values_map, by = c("phase" = "full_name"))
    }
  )
)
