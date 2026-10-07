# Taxonomic assignment: one computation per (method x database x input) runs
# in parallel, every (parameter) row of the benchmark is derived from it, then
# `d_all_taxo` stacks the rows (ROADMAP S8.2, docs/experimental_design.md §8).
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
source(here("R/itsx_inputs.R"))
source(here("R/assign_compute.R"))
source(here("R/otu_taxonomy.R"))
source(here("R/dada2_patched.R"))
# dada2 overflow fix (HISTORY.md §2, 2026-09-30): refuse to run with the
# Bioconductor build. Outside any target, so no target depends on it.
assert_dada2_patched()
load_pqverse(c("MiscMetabar", "comparpq", "tidypq")) # tidypq: combine_taxo_assignments() in the main process

tar_option_set(
  seed = targets_seed,
  # Workers attach these CRAN packages; the pqverse checkouts are loaded inside
  # each target with load_pqverse() (they are not installed packages).
  packages = c("here", "phyloseq", "dplyr", "tidyr", "tibble"),
  # config.R::single_worker: one worker for every computation, relaunched
  # after each task (tasks_max = 1) so that the memory of a dada2 computation
  # is given back before the next one starts.
  controller = if (single_worker) {
    crew::crew_controller_local(
      name = "single_ctrl",
      workers = 1,
      tasks_max = 1,
      seconds_idle = 60
    )
  } else {
    crew::crew_controller_group(
      crew::crew_controller_local(
        name = "dada2_ctrl",
        workers = 1,
        seconds_idle = 60
      ),
      crew::crew_controller_local(
        name = "fast_ctrl",
        workers = n_workers,
        seconds_idle = 60
      )
    )
  }
)

# One input per ITSx region of the dataset and per unit (OD 34):
# config.R::itsx_regions x itsx_units ("otu" only when the dataset has an OTU
# input, OD 30c). The raw ASVs and OTUs are no longer assigned.
values_map <- build_values_map(dbs = db_list, mini_db = mini_db)
inputs_map <- dplyr::distinct(values_map, input, region, unit, input_pq)

# One autometric log file per target and per run (see R/autometric_helpers.R),
# in one folder per dataset. The phase names are the same for every dataset
# (`compute_dada2__EUK_ITS_v2.1`, ...), and summarise_autometric_costs() keeps
# the newest run of each phase, so a shared folder lets one dataset's costs
# replace another's: the in silico run of 2026-09-22 wrote 97 files into the
# mock's folder before this was split (2026-09-23). The production project
# keeps its folder name, so its targets are not invalidated.
# Keyed on the project name: the production folders keep their names
# (assign_taxo, assign_taxo_<dataset>), and a *_mini project gets its own
# instead of writing into the production one (code review §3.1).
autometric_dir_assign <- here(
  "data/data_final/autometric",
  Sys.getenv("TAR_PROJECT", unset = "assign_taxo")
)

# One computation per method × database × input (R/assign_compute.R): 36 per
# input with the 9 databases of config.R (4 methods), 189 assignment rows per
# input (21 per database). Threads per method: config.R::assign_threads (M2.4).
compute_map <- dplyr::distinct(
  values_map,
  compute_name,
  method,
  db,
  input,
  input_pq,
  ref_file,
  controller
)

compute_targets <- tarchetypes::tar_eval(
  tar_target(
    compute_name,
    {
      dada2_build_on_worker # every computation waits for the worker check
      load_pqverse(c("MiscMetabar", "comparpq"))
      with_autometric(
        compute_name,
        compute_assignment(
          input_pq_sym,
          method = method,
          ref_fasta = ref_file_sym,
          nproc = assign_threads[[method]]
        ),
        dir = autometric_dir_assign
      )
    },
    resources = if (single_worker) {
      tar_resources()
    } else {
      tar_resources(crew = tar_resources_crew(controller = controller))
    }
  ),
  dplyr::mutate(
    compute_map,
    ref_file_sym = rlang::syms(ref_file),
    input_pq_sym = rlang::syms(input_pq)
  )
)

# The benchmark rows: thresholds, lca cutoffs or blastn filters and votes
# applied to their computation, in the main process. Each input keeps its rows
# apart: the inputs carry different taxa (merged units, fakes), and
# combine_taxo_assignments() needs the taxa of its base phyloseq.
assignment_targets_for <- function(map) {
  tarchetypes::tar_eval(
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
        min_cover = min_cover,
        min_id = min_id,
        lca_cutoff = lca_cutoff
      ),
      deployment = "main"
    ),
    dplyr::mutate(
      map,
      input_pq_sym = rlang::syms(input_pq),
      compute_sym = rlang::syms(compute_name)
    )
  )
}
assignment_targets <- lapply(
  stats::setNames(inputs_map$input, inputs_map$input),
  \(i) assignment_targets_for(dplyr::filter(values_map, input == i))
)

# One file target per reference fasta (ROADMAP S7.9): a database rebuilt in
# place, under the same name, invalidates every assignment that reads it.
ref_file_targets <- tarchetypes::tar_eval(
  tar_target(ref_file, here(db_path), format = "file", deployment = "main"),
  dplyr::distinct(values_map, ref_file, db_path)
)


# The ITSx inputs of one region (OD 34, R/itsx_inputs.R): the real ASVs
# extracted by MiscMetabar::itsx_pq() (ASVs identical once trimmed are merged
# into the most abundant, undetected ASVs keep their full sequence, the ITSx
# origin goes to `ITSx_origin`, the detection to `ITSx_detected`), the external controls extracted alone, then
# the controls added: `fake_` shuffled from the trimmed units, `external_`
# trimmed, their origin written back. With an OTU input, the units are the
# 97 % vsearch clusters of the trimmed ASVs (asv2otu(), no MUMU), and
# `otu_clusters_<region>` records the membership of each ASV (Q5,
# OTU_derived).
regions_map <- tibble::tibble(region = itsx_regions) |>
  dplyr::mutate(
    real_asv = rlang::syms(paste0("d_", region, "_asv_real")),
    controls = rlang::syms(paste0("itsx_controls_", region)),
    for_asv = rlang::syms(paste0("d_", region, "_asv_for_assignation")),
    real_otu = rlang::syms(paste0("d_", region, "_otu_real")),
    clusters = rlang::syms(paste0("otu_clusters_", region)),
    for_otu = rlang::syms(paste0("d_", region, "_otu_for_assignation"))
  )
asv_input_targets <- tarchetypes::tar_eval(
  list(
    tar_target(
      real_asv,
      MiscMetabar::itsx_pq(
        d_asv,
        region = region,
        organism_groups = itsx_organism_groups,
        keep_undetected = TRUE,
        add_origin = TRUE,
        add_detected = TRUE,
        duplicated_seqs = "merge",
        nproc = n_threads,
        args_before_itsx = itsx_conda_prelude,
        verbose = FALSE
      ),
      deployment = "main"
    ),
    tar_target(
      controls,
      itsx_controls(
        external_seqs,
        region = region,
        organism_groups = itsx_organism_groups,
        nproc = n_threads,
        args_before_itsx = itsx_conda_prelude
      ),
      deployment = "main"
    ),
    tar_target(
      for_asv,
      with_controls(
        real_asv,
        controls,
        prop_fake = prop_fake,
        seed = targets_seed
      ),
      deployment = "main"
    )
  ),
  values = regions_map
)
otu_input_targets <- if (run_otu_input) {
  tarchetypes::tar_eval(
    list(
      tar_target(
        real_otu,
        MiscMetabar::asv2otu(real_asv, method = "vsearch", tax_adjust = 0),
        deployment = "main"
      ),
      tar_target(clusters, otu_membership(real_asv), deployment = "main"),
      tar_target(
        for_otu,
        with_controls(
          real_otu,
          controls,
          prop_fake = prop_fake,
          seed = targets_seed
        ),
        deployment = "main"
      )
    ),
    values = regions_map
  )
} else {
  list()
}

# One object per input: its phyloseq with every assignment column. Production
# refuses an assignment that added no column (ROADMAP S1.3); smoke tests on
# the mini_* databases tolerate it (blastn often has no hit).
d_all_taxo_names <- paste0("d_all_taxo_", inputs_map$input)
d_all_taxo_targets <- lapply(seq_len(nrow(inputs_map)), function(k) {
  tarchetypes::tar_combine_raw(
    d_all_taxo_names[[k]],
    assignment_targets[[inputs_map$input[[k]]]],
    command = substitute(
      combine_taxo_assignments(BASE, !!!.x, allow_empty = mini_db),
      list(BASE = as.name(inputs_map$input_pq[[k]]))
    ),
    deployment = "main"
  )
})
# Expression that only names every d_all_taxo_<input>, so the targets below
# run once all of them are built.
wait_all_inputs <- as.call(c(as.name("{"), lapply(d_all_taxo_names, as.name)))

tar_plan(
  # assert_dada2_patched() above runs in the main process, but dada2 runs on
  # the crew worker, which gets the patched build only through .Rprofile
  # (code review §3.4): check it there, before any computation.
  tar_target(dada2_build_on_worker, {
    assert_dada2_patched()
    find.package("dada2")
  }),
  # File dependency on the DADA2 store so a rebuilt d_asv invalidates
  # everything downstream (ROADMAP S1.2).
  tar_target(
    d_asv_file,
    here(dada2_store, "objects/d_asv"),
    format = "file",
    deployment = "main"
  ),
  tar_target(d_asv, readRDS(d_asv_file), deployment = "main"),
  # fake_ref_file is a file target so a regenerated fake reference
  # (make_databases.R::derive_fake_ref()) invalidates every assignment.
  tar_target(
    fake_ref_file,
    here(fake_ref_fasta),
    format = "file",
    deployment = "main"
  ),
  tar_target(
    external_seqs,
    Biostrings::readDNAStringSet(fake_ref_file),
    deployment = "main"
  ),
  asv_input_targets,
  otu_input_targets,
  ref_file_targets,
  compute_targets,
  assignment_targets,
  d_all_taxo_targets,
  tar_target_raw(
    "benchmark_costs",
    substitute(
      {
        WAIT # aggregate only once every assignment has run
        # Costs are logged per computation (compute_name). right_join keeps
        # one row per benchmark row: the rows derived from one computation
        # share its cost, and `phase` carries the full_name the chapters join
        # on. add_tree_costs(): memory and CPU of the external programs too
        # (vsearch, BLAST+), which autometric does not see (code review §1.3).
        # Computations shorter than the 1 s sampling interval get NA
        # autometric columns but keep `wall_s` and the CPU seconds.
        read_autometric_dir(autometric_dir_assign) |>
          summarise_autometric_costs() |>
          add_tree_costs(autometric_dir_assign) |>
          dplyr::right_join(values_map, by = c("phase" = "compute_name")) |>
          dplyr::mutate(compute_name = phase, phase = full_name)
      },
      list(WAIT = wait_all_inputs)
    ),
    deployment = "main"
  ),
  # Versions of the pqverse checkouts and of R used for this run (decision 10).
  tar_target_raw(
    "session_info",
    substitute(
      {
        WAIT
        list(
          pqverse = pqverse_versions(),
          project = project_version(),
          session = sessioninfo::session_info()
        )
      },
      list(WAIT = wait_all_inputs)
    ),
    deployment = "main"
  )
)
