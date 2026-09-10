# Builders for the (method x database) grids consumed by the pipelines and the
# analysis notebook (ROADMAP S2.1). Replaces the three hand-copied
# `values_map` blocks (the former script_assign_taxo_parallel.R, R/values_map_for_qmd.R
# and script_cross_val.R).
#
# Requires `db_list` and `mini_db` from config.R. The `full_name` strings are
# the targets names in the stores: changing `full_name_for()` renames every
# target and forces a full rerun. tests/test_values_map.R pins them.

# Assignment methods and their parameter sweeps.
build_methods_grid <- function(min_bootstrap = c(0.4, 0.5, 0.6),
                               vote_algorithm = c("rel_majority", "abs_majority", "unanimity"),
                               nb_voting = 100) {
  dplyr::full_join(
    tidyr::expand_grid(
      method = c("dada2", "sintax", "lca"),
      min_bootstrap = min_bootstrap
    ),
    tidyr::expand_grid(
      method = "blastn",
      vote_algorithm = vote_algorithm,
      nb_voting = nb_voting,
      min_bootstrap = 0.5
    ),
    by = c("method", "min_bootstrap")
  )
}

# dada2 reads the dada2_format/ fasta, every other method the sintax_format/.
db_path_for <- function(method, db_name) {
  paste0(
    ifelse(
      method == "dada2",
      paste0("data/data_raw/refseq/dada2_format/", db_name),
      paste0("data/data_raw/refseq/sintax_format/", db_name)
    ),
    ".fasta"
  )
}

# Target name of one assignment. Blastn rows carry the vote suffix; the other
# methods stop after min_bootstrap. Vectorised.
full_name_for <- function(method, db, min_bootstrap, vote_algorithm = NA, nb_voting = NA) {
  base <- paste0(method, "__", db, "___", min_bootstrap)
  ifelse(
    is.na(vote_algorithm),
    base,
    paste0(base, "...", vote_algorithm, "...", nb_voting)
  )
}

# Grid for pipelines/assign_taxo.R and the notebook.
build_values_map <- function(dbs = db_list,
                             mini_db = FALSE,
                             methods = build_methods_grid()) {
  tidyr::expand_grid(methods, db = dbs) |>
    dplyr::mutate(
      cutadapted_db = ifelse(grepl("cut", db), "cut", ""),
      db_filter     = ifelse(grepl("Fungi", db), "Fungi", ""),
      db_name       = ifelse(mini_db, paste0("mini_", db), db),
      do_clean_pq   = method == "dada2",
      db_path       = db_path_for(method, db_name),
      controller    = ifelse(method == "dada2", "dada2_ctrl", "fast_ctrl"),
      full_name     = full_name_for(method, db, min_bootstrap, vote_algorithm, nb_voting)
    )
}

# Grid for pipelines/cross_val.R: one bootstrap value, blastn with rel_majority,
# and the leaked / standard variants.
build_cv_values_map <- function(dbs = db_list,
                                mini_db = FALSE,
                                remove_tested = c(TRUE, FALSE),
                                min_bootstrap = 0.5) {
  cv_methods <- dplyr::bind_rows(
    tidyr::expand_grid(
      method = c("dada2", "sintax", "lca"),
      min_bootstrap = min_bootstrap
    ),
    tibble::tibble(
      method = "blastn",
      min_bootstrap = min_bootstrap,
      vote_algorithm = "rel_majority",
      nb_voting = 100L
    )
  )
  tidyr::expand_grid(cv_methods, db = dbs, remove_tested = remove_tested) |>
    dplyr::mutate(
      db_name       = ifelse(mini_db, paste0("mini_", db), db),
      db_path       = db_path_for(method, db_name),
      leaked_suffix = ifelse(remove_tested, "standard", "leaked"),
      full_name     = paste0("cv__", method, "__", db, "__", leaked_suffix)
    )
}
