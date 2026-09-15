# k-fold cross-validation of every (method x database) pair on the database
# itself, standard (held-out sequences removed) and leaked variants.
#
# WARNING: the full run (cv_fold_number folds x 4 methods x length(db_list) DBs
# x 2 variants) takes hours. The cross_val_mini project (mini_* databases, own
# store) uses smoke-test values set in config.R: 2 folds, 200 sequences.
#
# The queries are trimmed to the ITS1F-ITS2 amplicon with cutadapt (config.R
# primers, ROADMAP D1a): cutadapt must run in the cutadaptenv conda env.
#
# Run with:
#   Sys.setenv(TAR_PROJECT = "cross_val"); targets::tar_make()

library("conflicted")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
library("dplyr")

here::i_am("pipelines/cross_val.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/values_map.R"))
source(here("R/cross_val.R"))
source(here("R/cv_queries.R"))
source(here("R/create_fake_pq_from_refseq.R"))
source(here("R/cv_to_tidy.R"))
load_pqverse(c("MiscMetabar", "comparpq", "tidypq", "dbpq"))

tar_option_set(
  seed = targets_seed,
  packages = c("here", "phyloseq", "dplyr", "tidyr", "tibble", "Biostrings"),
  controller = crew::crew_controller_local(workers = n_workers, seconds_idle = 60)
)

cv_values_map <- build_cv_values_map(dbs = db_list, mini_db = mini_db)

# Runs cross_val() and forwards blastn-specific args only when they are not NA.
# vote_algorithm / nb_voting are NA for non-blastn methods (tar_eval substitutes
# them literally from cv_values_map).
run_cv <- function(method, db_path, fold_number, fold_tested, min_bootstrap,
                   remove_tested, seed, vote_algorithm = NULL, nb_voting = NULL,
                   max_seq = NULL, min_cover = NULL, nproc = 1,
                   primer_fw = NULL, primer_rev = NULL, primer_min_overlap = NULL,
                   oversample = 2, cutadapt_prelude = NULL, query_fasta = NULL,
                   id_fasta = NULL, query_id_fasta = NULL) {
  extra_args <- list()
  if (!is.null(vote_algorithm) && !is.na(vote_algorithm)) {
    extra_args$vote_algorithm <- vote_algorithm
    extra_args$nb_voting <- as.integer(nb_voting)
  }
  # cross_val() passes `...` unfiltered to every method: min_cover only for blastn.
  if (method == "blastn" && !is.null(min_cover) && !is.na(min_cover)) {
    extra_args$min_cover <- min_cover
  }
  do.call(cross_val, c(
    list(
      ref_fasta = db_path,
      fold_number = fold_number,
      fold_tested = fold_tested,
      method = method,
      min_bootstrap = min_bootstrap,
      remove_tested_sequences = remove_tested,
      seed = seed,
      max_seq = max_seq,
      nproc = nproc,
      primer_fw = primer_fw,
      primer_rev = primer_rev,
      primer_min_overlap = primer_min_overlap,
      oversample = oversample,
      cutadapt_prelude = cutadapt_prelude,
      query_fasta = query_fasta,
      id_fasta = id_fasta,
      query_id_fasta = query_id_fasta
    ),
    extra_args
  ))
}

cv_targets <- tarchetypes::tar_eval(
  tar_target(
    full_name,
    {
      load_pqverse(c("MiscMetabar", "comparpq", "tidypq", "dbpq"))
      cv_to_tidy(
        run_cv(
          method = method,
          db_path = ref_file_sym,
          fold_number = cv_fold_number,
          fold_tested = cv_fold_tested,
          min_bootstrap = min_bootstrap,
          remove_tested = remove_tested,
          seed = targets_seed,
          vote_algorithm = vote_algorithm,
          nb_voting = nb_voting,
          max_seq = cv_max_seq,
          min_cover = min_cover,
          nproc = cv_threads,
          primer_fw = fw_primer_sequences,
          primer_rev = rev_primer_sequences,
          primer_min_overlap = primer_min_overlap,
          oversample = cv_oversample,
          cutadapt_prelude = cutadapt_conda_prelude,
          query_fasta = query_ref_file_sym,
          id_fasta = id_ref_file_sym,
          query_id_fasta = query_id_ref_file_sym
        ),
        method = method,
        db = db,
        remove_tested = remove_tested,
        min_bootstrap = min_bootstrap
      )
    }
  ),
  values = dplyr::mutate(
    cv_values_map,
    ref_file_sym = rlang::syms(ref_file),
    query_ref_file_sym = rlang::syms(query_ref_file),
    id_ref_file_sym = rlang::syms(id_ref_file),
    query_id_ref_file_sym = rlang::syms(query_id_ref_file)
  )
)

# One file target per reference fasta (ROADMAP S7.9): a database rebuilt in
# place, under the same name, invalidates every CV target that reads it.
ref_file_targets <- tarchetypes::tar_eval(
  tar_target(ref_file, here(db_path), format = "file", deployment = "main"),
  # References and query sources (the _Fungi file of a _Fungi_cut database).
  values = dplyr::distinct(dplyr::bind_rows(
    dplyr::distinct(cv_values_map, ref_file, db_path),
    dplyr::distinct(cv_values_map, ref_file = query_ref_file, db_path = query_db_path),
    dplyr::distinct(cv_values_map, ref_file = id_ref_file, db_path = id_db_path),
    dplyr::distinct(cv_values_map, ref_file = query_id_ref_file, db_path = query_id_db_path)
  ))
)

tar_plan(
  ref_file_targets,
  cv_targets,
  tarchetypes::tar_combine(
    cv_results,
    cv_targets,
    command = dplyr::bind_rows(!!!.x),
    deployment = "main"
  ),
  tar_target(
    session_info,
    {
      cv_results
      list(
        pqverse = pqverse_versions(),
        session = sessioninfo::session_info()
      )
    },
    deployment = "main"
  )
)
