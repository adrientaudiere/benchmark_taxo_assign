# Taxonomic assignment: every (method x database x parameter) combination runs
# in parallel on `d_asv_for_assignation`, then `d_all_taxo` stacks the results.
#
# Run with:
#   Sys.setenv(TAR_PROJECT = "assign_taxo"); targets::tar_make()
# Smoke test on the mini_* databases: TAR_PROJECT = "assign_taxo_mini"
# (own store store_assign_taxo_mini, see _targets.yaml and config.R).

library("conflicted")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
library("dplyr")

here::i_am("pipelines/assign_taxo.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/values_map.R"))
source(here("R/combine_taxo_assignments.R"))
source(here("R/autometric_helpers.R"))
load_pqverse(c("MiscMetabar", "comparpq"))

tar_option_set(
  seed = targets_seed,
  # Workers attach these CRAN packages; the pqverse checkouts are loaded inside
  # each target with load_pqverse() (they are not installed packages).
  packages = c("here", "phyloseq", "dplyr", "tidyr", "tibble"),
  controller = crew::crew_controller_group(
    crew::crew_controller_local(name = "dada2_ctrl", workers = 1,         seconds_idle = 60),
    crew::crew_controller_local(name = "fast_ctrl",  workers = n_workers, seconds_idle = 60)
  )
)

values_map <- build_values_map(dbs = db_list, mini_db = mini_db)

# One autometric log file per target and per run (see R/autometric_helpers.R).
autometric_dir_assign <- here("data/data_final/autometric/assign_taxo")

assignment_targets <- tarchetypes::tar_eval(
  tar_target(
    full_name,
    {
      load_pqverse(c("MiscMetabar", "comparpq"))
      with_autometric(
        full_name,
        add_new_taxonomy_pq(
          d_asv_for_assignation,
          method = method,
          ref_fasta = db_path,
          suffix = paste0("_", full_name),
          min_bootstrap = min_bootstrap,
          vote_algorithm = vote_algorithm,
          nb_voting = nb_voting
        ),
        dir = autometric_dir_assign
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = controller)
    )
  ),
  values_map
)

tar_plan(
  # File dependency on the DADA2 store so a rebuilt d_asv invalidates
  # everything downstream (ROADMAP S1.2).
  tar_target(d_asv_file, here("store_dada2/objects/d_asv"), format = "file",
             deployment = "main"),
  tar_target(d_asv, readRDS(d_asv_file), deployment = "main"),
  tar_target(
    file_taxo_mock,
    here(taxo_mock_csv),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    taxo_mock,
    read.csv(file_taxo_mock) |>
      select(any_of(c(
        "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
      ))) |>
      magrittr::set_rownames(read.csv(file_taxo_mock)$MockStrain),
    deployment = "main"
  ),
  # Negative controls: shuffled ASVs (fake_*) then external non-Fungi
  # sequences (external_*); both feed the TN denominator of tc_metrics_mock().
  tar_target(
    d_asv_shuffled,
    add_shuffle_seq_pq(d_asv, prop_fake = prop_fake),
    deployment = "main"
  ),
  tar_target(
    d_asv_for_assignation,
    add_external_seq_pq(
      d_asv_shuffled,
      Biostrings::readDNAStringSet(here(fake_ref_fasta))
    ),
    deployment = "main"
  ),
  assignment_targets,
  # Production refuses an assignment that added no column (ROADMAP S1.3);
  # smoke tests on the mini_* databases tolerate it (blastn often has no hit).
  tarchetypes::tar_combine(
    d_all_taxo,
    assignment_targets,
    command = combine_taxo_assignments(d_asv_for_assignation, !!!.x,
                                       allow_empty = mini_db),
    deployment = "main"
  ),
  tar_target(
    benchmark_costs,
    {
      d_all_taxo # aggregate only once every assignment has run
      read_autometric_dir(autometric_dir_assign) |>
        summarise_autometric_costs() |>
        dplyr::inner_join(values_map, by = c("phase" = "full_name"))
    },
    deployment = "main"
  ),
  # Versions of the pqverse checkouts and of R used for this run (decision 10).
  tar_target(
    session_info,
    {
      d_all_taxo
      list(
        pqverse = pqverse_versions(),
        session = sessioninfo::session_info()
      )
    },
    deployment = "main"
  )
)
