# k-fold cross-validation of every (method x database) pair on the database
# itself, standard (held-out sequences removed) and leaked variants.
#
# WARNING: the full run (cv_fold_number folds x 4 methods x 6 DBs x 2 variants)
# takes many hours. Lower cv_fold_tested / cv_max_seq in config.R for a smoke
# test, or use the cross_val_mini project (mini_* databases, own store).
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
source(here("R/create_fake_pq_from_refseq.R"))
source(here("R/cv_to_tidy.R"))
load_pqverse(c("MiscMetabar", "comparpq"))

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
                   max_seq = NULL) {
  extra_args <- list()
  if (!is.null(vote_algorithm) && !is.na(vote_algorithm)) {
    extra_args$vote_algorithm <- vote_algorithm
    extra_args$nb_voting <- as.integer(nb_voting)
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
      max_seq = max_seq
    ),
    extra_args
  ))
}

cv_targets <- tarchetypes::tar_eval(
  tar_target(
    full_name,
    {
      load_pqverse(c("MiscMetabar", "comparpq"))
      cv_to_tidy(
        run_cv(
          method = method,
          db_path = db_path,
          fold_number = cv_fold_number,
          fold_tested = cv_fold_tested,
          min_bootstrap = min_bootstrap,
          remove_tested = remove_tested,
          seed = targets_seed,
          vote_algorithm = vote_algorithm,
          nb_voting = nb_voting,
          max_seq = cv_max_seq
        ),
        method = method,
        db = db,
        remove_tested = remove_tested,
        min_bootstrap = min_bootstrap
      )
    }
  ),
  values = cv_values_map
)

tar_plan(
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
