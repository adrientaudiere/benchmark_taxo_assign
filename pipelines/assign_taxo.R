# Taxonomic assignment: one computation per (method x database x input) runs
# in parallel, every (parameter) row of the benchmark is derived from it, then
# `d_all_taxo` stacks the rows (ROADMAP S8.2, docs/compute_budget.md).
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
source(here("R/itsx.R"))
source(here("R/assign_compute.R"))
source(here("R/otu_taxonomy.R"))
load_pqverse(c("MiscMetabar", "comparpq", "tidypq"))  # tidypq: combine_taxo_assignments() in the main process

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

# One computation per method × database × input (R/assign_compute.R): 44
# computations for the 132 assignments. dada2 runs on its own worker with
# assign_threads_dada2 threads, sintax / lca / blastn on fast_ctrl.
compute_map <- dplyr::distinct(
  values_map, compute_name, method, db, preprocess, input_pq, ref_file, controller
)

compute_targets <- tarchetypes::tar_eval(
  tar_target(
    compute_name,
    {
      load_pqverse(c("MiscMetabar", "comparpq"))
      with_autometric(
        compute_name,
        compute_assignment(
          input_pq_sym,
          method = method,
          ref_fasta = ref_file_sym,
          nproc = if (method == "dada2") assign_threads_dada2 else assign_threads_fast
        ),
        dir = autometric_dir_assign
      )
    },
    resources = tar_resources(
      crew = tar_resources_crew(controller = controller)
    )
  ),
  dplyr::mutate(
    compute_map,
    ref_file_sym = rlang::syms(ref_file),
    input_pq_sym = rlang::syms(input_pq)
  )
)

# The 132 benchmark rows, with the target names used before S8.2: thresholds,
# votes or lca copies applied to their computation, in the main process.
assignment_targets <- tarchetypes::tar_eval(
  tar_target(
    full_name,
    derive_assignment(
      input_pq_sym,
      compute_sym,
      method = method,
      suffix = paste0("_", full_name),
      min_bootstrap = min_bootstrap,
      vote_algorithm = vote_algorithm,
      nb_voting = nb_voting,
      min_cover = min_cover
    ),
    deployment = "main"
  ),
  dplyr::mutate(
    values_map,
    input_pq_sym = rlang::syms(input_pq),
    compute_sym = rlang::syms(compute_name)
  )
)

# One file target per reference fasta (ROADMAP S7.9): a database rebuilt in
# place, under the same name, invalidates every assignment that reads it.
ref_file_targets <- tarchetypes::tar_eval(
  tar_target(ref_file, here(db_path), format = "file", deployment = "main"),
  dplyr::distinct(values_map, ref_file, db_path)
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
  # ITS1 of every ASV with ITSx (ROADMAP Q4, R/itsx.R), used only to trim the
  # flanking regions. An ASV whose ITS1 equals the ITS1 of a more abundant ASV
  # is removed from both inputs (decision 20), so the raw and ITSx inputs keep
  # the same taxa; undetected ASVs keep their full sequence. Assignments
  # reading the ITSx input are the "itsx_" rows of values_map.
  tar_target(
    itsx_asv,
    run_itsx(
      d_asv@refseq,
      region = itsx_region,
      organism_groups = itsx_organism_groups,
      cpu = n_threads,
      prelude = itsx_conda_prelude
    ),
    deployment = "main"
  ),
  tar_target(itsx_dropped_taxa, itsx_duplicated_taxa(d_asv, itsx_asv), deployment = "main"),
  tar_target(
    d_asv_common,
    phyloseq::prune_taxa(setdiff(phyloseq::taxa_names(d_asv), itsx_dropped_taxa), d_asv),
    deployment = "main"
  ),
  tar_target(d_asv_itsx, itsx_replace_refseq(d_asv_common, itsx_asv), deployment = "main"),
  # Post-clustering (ROADMAP Q5): membership of the ASVs of both inputs in
  # 97 % vsearch clusters. No assignment is rerun on OTUs: chapter 01 builds
  # each OTU's taxonomy from its member ASVs with resolve_otu_taxonomy().
  tar_target(otu_clusters, otu_membership(d_asv_common), deployment = "main"),
  # Negative controls: shuffled ASVs (fake_*) then external non-Fungi
  # sequences (external_*); both feed the TN denominator of tc_metrics_mock().
  # fake_ref_file is a file target so a regenerated fake reference
  # (make_databases.R::derive_fake_ref()) invalidates every assignment.
  # with_seed(): both inputs shuffle the same ASVs, so they carry the same
  # fake_1..n and external_* taxa (combine_taxo_assignments() requires it).
  tar_target(fake_ref_file, here(fake_ref_fasta), format = "file",
             deployment = "main"),
  tar_target(
    d_asv_shuffled,
    withr::with_seed(targets_seed, add_shuffle_seq_pq(d_asv_common, prop_fake = prop_fake)),
    deployment = "main"
  ),
  tar_target(
    d_asv_for_assignation,
    add_external_seq_pq(
      d_asv_shuffled,
      Biostrings::readDNAStringSet(fake_ref_file)
    ),
    deployment = "main"
  ),
  tar_target(
    d_asv_itsx_shuffled,
    withr::with_seed(targets_seed, add_shuffle_seq_pq(d_asv_itsx, prop_fake = prop_fake)),
    deployment = "main"
  ),
  tar_target(
    d_asv_itsx_for_assignation,
    add_external_seq_pq(
      d_asv_itsx_shuffled,
      Biostrings::readDNAStringSet(fake_ref_file)
    ),
    deployment = "main"
  ),
  ref_file_targets,
  compute_targets,
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
      # Costs are logged per computation (compute_name). right_join keeps one
      # row per benchmark row: the rows derived from one computation share its
      # cost, and `phase` carries the full_name the chapters join on.
      # Computations shorter than the 1 s sampling interval get NA costs.
      read_autometric_dir(autometric_dir_assign) |>
        summarise_autometric_costs() |>
        dplyr::right_join(values_map, by = c("phase" = "compute_name")) |>
        dplyr::mutate(compute_name = phase, phase = full_name)
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
