# WARNING: full cross-validation (cv_fold_number folds x 4 methods x 9 DBs x 2
# variants) takes many hours. Keep cv_fold_tested = 2 in config.R for smoke-
# testing; raise it to cv_fold_number only for the publication run.
#
# Run with:
#   tar_make(script = "script_cross_val.R", store = "store_cross_val")

library("conflicted")
library("MiscMetabar")
devtools::load_all("/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/MiscMetabar/")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
library("dplyr")

here::i_am("script_cross_val.R")
source(here("config.R"))
source(here("R/cross_val.R"))
source(here("R/functions.R"))
source(here("R/cv_to_tidy.R"))
source("/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/comparpq/R/compare_taxo.R")
source("/home/adrien/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg/comparpq/R/fake_creation.R")

tar_option_set(
  seed = targets_seed,
  controller = crew::crew_controller_local(workers = n_workers, seconds_idle = 60)
)

cv_methods <- dplyr::bind_rows(
  tidyr::expand_grid(
    method = c("dada2", "sintax", "lca"),
    min_bootstrap = 0.5
  ),
  tibble::tibble(
    method = "blastn",
    min_bootstrap = 0.5,
    vote_algorithm = "rel_majority",
    nb_voting = 100L
  )
)

cv_values_map <- tidyr::expand_grid(
  cv_methods,
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
  ),
  remove_tested = c(TRUE, FALSE)
) |>
  dplyr::mutate(
    db_path = paste0(ifelse(
      method == "dada2",
      paste0("data/data_raw/refseq/dada2_format/", db),
      paste0("data/data_raw/refseq/sintax_format/", db)
    ), ".fasta"),
    leaked_suffix = ifelse(remove_tested, "standard", "leaked"),
    full_name = paste0(
      "cv__", method, "__", db, "__", leaked_suffix
    )
  )

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
  ),
  values = cv_values_map
)

tar_plan(
  cv_targets,
  tarchetypes::tar_combine(
    cv_results,
    cv_targets,
    command = dplyr::bind_rows(!!!.x)
  )
)
