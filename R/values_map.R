# Builders for the (method x database) grids consumed by the pipelines and the
# analysis notebook (ROADMAP S2.1). Replaces the three hand-copied
# `values_map` blocks (the former script_assign_taxo_parallel.R, R/values_map_for_qmd.R
# and script_cross_val.R).
#
# Requires `db_list`, `cv_db_list`, `itsx_regions`, `itsx_units`, `mini_db` and `blastn_min_cover` from config.R. The
# `min_cover` column is set on blastn rows only (NA elsewhere, dropped by
# add_new_taxonomy_pq() for the other methods), so changing
# `blastn_min_cover` reruns the blastn targets only. The `full_name` strings are
# the targets names in the stores: changing `full_name_for()` renames every
# target and forces a full rerun. tests/test_values_map.R pins them.

# Assignment methods and their parameter sweeps (docs/objectives_design.md
# §1.5): bootstrap threshold for dada2 and sintax, `lca_cutoff` for lca, vote
# x `min_id` for blastn. Every row of one method is derived from a single
# computation (R/assign_compute.R). blastn rows keep min_bootstrap = 0.5, a
# placeholder that only exists in their target name.
build_methods_grid <- function(
  min_bootstrap = c(0.4, 0.5, 0.6),
  lca_cutoff = c(0.51, 0.8, 1),
  vote_algorithm = c("rel_majority", "abs_majority", "unanimity"),
  nb_voting = 100,
  min_id = c(90, 92, 95, 97),
  min_cover = blastn_min_cover
) {
  dplyr::bind_rows(
    tidyr::expand_grid(
      method = c("dada2", "sintax"),
      min_bootstrap = min_bootstrap
    ),
    tidyr::expand_grid(method = "lca", lca_cutoff = lca_cutoff),
    tidyr::expand_grid(
      method = "blastn",
      min_bootstrap = 0.5,
      vote_algorithm = vote_algorithm,
      nb_voting = nb_voting,
      min_id = min_id,
      min_cover = min_cover
    )
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

# Target name of one assignment: "<input>_<method>__<db>___<parameter>". The
# input is "<region>_<unit>" (OD 34: "ITS1_asv", "ITS2_otu", "full_asv"...);
# the parameter is min_bootstrap, or lca_cutoff for lca rows
# ("ITS1_asv_lca__DB___0.8"). blastn rows add their vote, nb_voting and min_id
# ("ITS2_asv_blastn__DB___0.5...rel_majority...100...95"). The chapters read
# the input, method and database from the columns of values_map, joined on
# full_name, not from the name (C22-S3). `input = NA` gives the bare name.
# Vectorised.
full_name_for <- function(
  method,
  db,
  min_bootstrap,
  vote_algorithm = NA,
  nb_voting = NA,
  input = NA,
  min_id = NA,
  lca_cutoff = NA
) {
  base <- paste0(
    ifelse(is.na(input), "", paste0(input, "_")),
    method,
    "__",
    db,
    "___",
    ifelse(method == "lca", lca_cutoff, min_bootstrap)
  )
  ifelse(
    is.na(vote_algorithm),
    base,
    paste0(
      base,
      "...",
      vote_algorithm,
      "...",
      nb_voting,
      ifelse(is.na(min_id), "", paste0("...", min_id))
    )
  )
}

# Name of the target that runs one method on one database and one input
# (ROADMAP S8.2). dada2 and sintax thresholds, blastn votes and the lca rows
# (lca ignores min_bootstrap) are applied afterwards to its result by the
# `full_name` targets, so one computation serves three rows. Vectorised.
compute_name_for <- function(method, db, input = NA) {
  paste0(
    ifelse(is.na(input), "", paste0(input, "_")),
    "compute_",
    method,
    "__",
    db
  )
}

# Name of the phyloseq target an assignment reads: the units of one input
# ("ITS1_asv": the ITSx-extracted ITS1 of the ASVs; "ITS1_otu": their 97 %
# vsearch clusters) with their negative controls. Vectorised.
input_pq_for <- function(input) {
  paste0("d_", input, "_for_assignation")
}

# Grid for pipelines/assign_taxo.R and the notebook (OD 34): every ITSx region
# of the dataset (config.R::itsx_regions) x every unit (itsx_units: "asv", and
# "otu" when the dataset has an OTU input) x every database. The raw ASVs and
# OTUs are no longer assigned; the ITSx input covers the nine databases,
# the full releases included (S8.5 reversed).
build_values_map <- function(
  dbs = db_list,
  mini_db = FALSE,
  methods = build_methods_grid(),
  regions = itsx_regions,
  units = itsx_units
) {
  tidyr::expand_grid(methods, db = dbs, region = regions, unit = units) |>
    dplyr::mutate(
      input = paste0(region, "_", unit),
      # `if`, not ifelse(): with a scalar mini_db, ifelse() returns one value
      # that is recycled to every row, so every target used the first database
      # (bug present from 2026-05 to 2026-09-10, ROADMAP S1.5).
      db_name = if (mini_db) paste0("mini_", db) else db,
      db_path = db_path_for(method, db_name),
      ref_file = ref_file_target_for(method, db),
      controller = ifelse(method == "dada2", "dada2_ctrl", "fast_ctrl"),
      input_pq = input_pq_for(input),
      compute_name = compute_name_for(method, db, input),
      full_name = full_name_for(
        method,
        db,
        min_bootstrap,
        vote_algorithm,
        nb_voting,
        input,
        min_id,
        lca_cutoff
      )
    )
}

# Grid for pipelines/cross_val.R: one bootstrap value, blastn with rel_majority,
# and the leaked / standard variants. Databases: config.R::cv_db_list.
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
