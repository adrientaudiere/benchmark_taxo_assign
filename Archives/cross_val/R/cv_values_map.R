# Grid of the cross-validation (Archives/cross_val/pipelines/cross_val.R),
# moved out of R/values_map.R on 2026-10-14 with the other CV helpers. It uses
# db_path_for() and ref_file_target_for() of R/values_map.R, and cv_db_list of
# config_cv.R: source config.R, config_cv.R and R/values_map.R first.

# Grid for pipelines/cross_val.R: one bootstrap value, blastn with rel_majority,
# and the leaked / standard variants. Databases: config_cv.R::cv_db_list.
build_cv_values_map <- function(
  dbs = cv_db_list,
  mini_db = FALSE,
  remove_tested = c(TRUE, FALSE),
  min_bootstrap = 0.5,
  min_cover = blastn_min_cover
) {
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
      db_name = if (mini_db) paste0("mini_", db) else db,
      db_path = db_path_for(method, db_name),
      ref_file = ref_file_target_for(method, db),
      leaked_suffix = ifelse(remove_tested, "standard", "leaked"),
      full_name = paste0("cv__", method, "__", db, "__", leaked_suffix),
      # A _Fungi_cut database no longer holds the primer sites: its CV queries
      # are drawn and trimmed from its _Fungi source (ROADMAP D1a). Other
      # databases are their own query source.
      # Unanchored on purpose: a trimmed variant may carry a suffix naming its
      # amplicon (`_Fungi_cut_full_ITS`), and an anchored `_Fungi_cut$` would
      # miss it, leaving the CV to draw its queries from a database that no
      # longer holds the primer sites — the very thing this rule prevents.
      query_db = ifelse(grepl("_Fungi_cut", db), sub("_cut.*$", "", db), db),
      query_db_path = db_path_for(
        method,
        if (mini_db) paste0("mini_", query_db) else query_db
      ),
      query_ref_file = ref_file_target_for(method, query_db),
      # dada2_format/ headers carry no identifier: cross_val() names the
      # records with the identifiers of the sintax_format file of the same
      # database (same records, same order; ROADMAP B22). For the other methods
      # these columns point to their own reference file. rep(): the helpers
      # use ifelse() on `method`, which must have one value per row.
      id_db_path = db_path_for(rep("sintax", length(db)), db_name),
      id_ref_file = ref_file_target_for(rep("sintax", length(db)), db),
      query_id_db_path = db_path_for(
        rep("sintax", length(db)),
        if (mini_db) paste0("mini_", query_db) else query_db
      ),
      query_id_ref_file = ref_file_target_for(
        rep("sintax", length(db)),
        query_db
      ),
      # crew controller of the target (ROADMAP B24, 2026-09-16). With
      # cv_reduce_reference = FALSE a dada2 target loads the whole database:
      # 23.5 GB on Unite_s_all_20250219, at least 45 GB on EUK_ITS_v2.1. One
      # dada2 worker at a time, the other methods on the fast controller, as
      # pipelines/assign_taxo.R does. `method` has one value per row here, so
      # ifelse() is safe (unlike the scalar case of S1.5 above).
      controller = ifelse(method == "dada2", "dada2_ctrl", "fast_ctrl")
    )
}
