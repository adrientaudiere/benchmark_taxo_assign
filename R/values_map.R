# Builders for the (method x database) grids consumed by the pipelines and the
# analysis notebook (ROADMAP S2.1). Replaces the three hand-copied
# `values_map` blocks (the former script_assign_taxo_parallel.R, R/values_map_for_qmd.R
# and script_cross_val.R).
#
# Requires `db_list`, `itsx_db_list`, `mini_db` and `blastn_min_cover` from config.R. The
# `min_cover` column is set on blastn rows only (NA elsewhere, dropped by
# add_new_taxonomy_pq() for the other methods), so changing
# `blastn_min_cover` reruns the blastn targets only. The `full_name` strings are
# the targets names in the stores: changing `full_name_for()` renames every
# target and forces a full rerun. tests/test_values_map.R pins them.

# Assignment methods and their parameter sweeps.
build_methods_grid <- function(min_bootstrap = c(0.4, 0.5, 0.6),
                               vote_algorithm = c("rel_majority", "abs_majority", "unanimity"),
                               nb_voting = 100,
                               min_cover = blastn_min_cover) {
  dplyr::full_join(
    tidyr::expand_grid(
      method = c("dada2", "sintax", "lca"),
      min_bootstrap = min_bootstrap
    ),
    tidyr::expand_grid(
      method = "blastn",
      vote_algorithm = vote_algorithm,
      nb_voting = nb_voting,
      min_bootstrap = 0.5,
      min_cover = min_cover
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

# Name of the file target that tracks the reference fasta read by `method` for
# database `db` (ROADMAP S7.9). Like full_name, it ignores mini_db: the *_mini
# projects have their own stores. Vectorised.
ref_file_target_for <- function(method, db) {
  paste0("ref_", ifelse(method == "dada2", "dada2", "sintax"), "__", db)
}

# Target name of one assignment. Blastn rows carry the vote suffix; the other
# methods stop after min_bootstrap. A preprocessed input ("itsx", ROADMAP Q4)
# prefixes the method ("itsx_sintax__..."), so the chapters parse it as a
# method of its own and the names of the unprocessed targets do not change.
# Vectorised.
full_name_for <- function(method, db, min_bootstrap, vote_algorithm = NA, nb_voting = NA,
                          preprocess = "none") {
  base <- paste0(
    ifelse(preprocess == "none", "", paste0(preprocess, "_")),
    method, "__", db, "___", min_bootstrap
  )
  ifelse(
    is.na(vote_algorithm),
    base,
    paste0(base, "...", vote_algorithm, "...", nb_voting)
  )
}

# Name of the target that runs one method on one database and one input
# (ROADMAP S8.2). dada2 and sintax thresholds, blastn votes and the lca rows
# (lca ignores min_bootstrap) are applied afterwards to its result by the
# `full_name` targets, so one computation serves three rows. Vectorised.
compute_name_for <- function(method, db, preprocess = "none") {
  paste0(
    ifelse(preprocess == "none", "", paste0(preprocess, "_")),
    "compute_", method, "__", db
  )
}

# Name of the phyloseq target an assignment reads: the ASVs as denoised, or a
# preprocessed copy with the same negative controls. Vectorised.
input_pq_for <- function(preprocess) {
  ifelse(
    preprocess == "none",
    "d_asv_for_assignation",
    paste0("d_asv_", preprocess, "_for_assignation")
  )
}

# Grid for pipelines/assign_taxo.R and the notebook. Preprocessed inputs
# ("itsx") are assigned against `itsx_dbs` only (ROADMAP S8.5); every database
# gets the raw ASVs.
build_values_map <- function(dbs = db_list,
                             mini_db = FALSE,
                             methods = build_methods_grid(),
                             preprocess = c("none", "itsx"),
                             itsx_dbs = itsx_db_list) {
  tidyr::expand_grid(methods, db = dbs, preprocess = preprocess) |>
    dplyr::filter(preprocess == "none" | db %in% itsx_dbs) |>
    dplyr::mutate(
      cutadapted_db = ifelse(grepl("cut", db), "cut", ""),
      db_filter     = ifelse(grepl("Fungi", db), "Fungi", ""),
      # `if`, not ifelse(): with a scalar mini_db, ifelse() returns one value
      # that is recycled to every row, so every target used the first database
      # (bug present from 2026-05 to 2026-09-10, ROADMAP S1.5).
      db_name       = if (mini_db) paste0("mini_", db) else db,
      do_clean_pq   = method == "dada2",
      db_path       = db_path_for(method, db_name),
      ref_file      = ref_file_target_for(method, db),
      controller    = ifelse(method == "dada2", "dada2_ctrl", "fast_ctrl"),
      input_pq      = input_pq_for(preprocess),
      compute_name  = compute_name_for(method, db, preprocess),
      full_name     = full_name_for(method, db, min_bootstrap, vote_algorithm, nb_voting, preprocess)
    )
}

# Grid for pipelines/cross_val.R: one bootstrap value, blastn with rel_majority,
# and the leaked / standard variants.
build_cv_values_map <- function(dbs = db_list,
                                mini_db = FALSE,
                                remove_tested = c(TRUE, FALSE),
                                min_bootstrap = 0.5,
                                min_cover = blastn_min_cover) {
  cv_methods <- dplyr::bind_rows(
    tidyr::expand_grid(
      method = c("dada2", "sintax", "lca"),
      min_bootstrap = min_bootstrap
    ),
    tibble::tibble(
      method = "blastn",
      min_bootstrap = min_bootstrap,
      vote_algorithm = "rel_majority",
      nb_voting = 100L,
      min_cover = min_cover
    )
  )
  tidyr::expand_grid(cv_methods, db = dbs, remove_tested = remove_tested) |>
    dplyr::mutate(
      # `if`, not ifelse(): with a scalar mini_db, ifelse() returns one value
      # that is recycled to every row, so every target used the first database
      # (bug present from 2026-05 to 2026-09-10, ROADMAP S1.5).
      db_name       = if (mini_db) paste0("mini_", db) else db,
      db_path       = db_path_for(method, db_name),
      ref_file      = ref_file_target_for(method, db),
      leaked_suffix = ifelse(remove_tested, "standard", "leaked"),
      full_name     = paste0("cv__", method, "__", db, "__", leaked_suffix),
      # A _Fungi_cut database no longer holds the primer sites: its CV queries
      # are drawn and trimmed from its _Fungi source (ROADMAP D1a). Other
      # databases are their own query source.
      query_db       = ifelse(grepl("_Fungi_cut$", db), sub("_cut$", "", db), db),
      query_db_path  = db_path_for(method, if (mini_db) paste0("mini_", query_db) else query_db),
      query_ref_file = ref_file_target_for(method, query_db),
      # dada2_format/ headers carry no identifier: cross_val() names the
      # records with the identifiers of the sintax_format file of the same
      # database (same records, same order; ROADMAP B22). For the other methods
      # these columns point to their own reference file. rep(): the helpers
      # use ifelse() on `method`, which must have one value per row.
      id_db_path        = db_path_for(rep("sintax", length(db)), db_name),
      id_ref_file       = ref_file_target_for(rep("sintax", length(db)), db),
      query_id_db_path  = db_path_for(
        rep("sintax", length(db)),
        if (mini_db) paste0("mini_", query_db) else query_db
      ),
      query_id_ref_file = ref_file_target_for(rep("sintax", length(db)), query_db)
    )
}
