# Shared by test_seed_spread.R and test_read_filtering.R (ROADMAP §1, code
# review §1.9 and §11 Q4, Q9, Q10): the columns of one computation, built as
# pipelines/assign_taxo.R builds them, and their scores on the Pauvert mock,
# built as analysis/01_load_and_clean.qmd builds them.
#
# Source after analysis/00_setup.R (config, values_map, the pqverse
# checkouts). Reads the stores, never writes to them.

source(here("R/mock_truth.R"))
source(here("R/score_assignments.R"))
source(here("R/assignment_helpers.R"))
source(here("R/assign_compute.R"))
source(here("R/combine_taxo_assignments.R"))
source(here("R/dada2_patched.R"))

# The two databases of the zoom on Pauvert (developer, 2026-10-06): the best
# EUKARYOME and the best UNITE at Genus F1, default settings.
zoom_test_dbs <- c("EUK_ITS_v2.1_Fungi_rep", "Unite_all_20250219")
pauvert_cfg <- mock_sanger_tables[mock_sanger_tables$dataset == "pauvert", ]
pauvert_vm <- values_map_for("pauvert") |>
  filter(input == "ITS1_asv")

# The seed targets gives a target of that name in the production project.
target_seed <- function(name) {
  targets::tar_seed_create(name, global_seed = targets_seed)
}

# compute_assignment() with the threads of the stored computations these
# tests compare with (3 for dada2, 1 for the others: the pipeline's setting
# until 2026-10-07, before config.R::assign_threads) and an explicit seed:
# set.seed() for dada2 (its bootstrap draws use R's generator), the vsearch
# --randseed for sintax (single thread). lca and blastn draw nothing.
run_computation <- function(input_pq, method, db, seed) {
  ref <- here(db_path_for(method, db))
  nproc <- c(dada2 = 3L, sintax = 1L, lca = 1L, blastn = 1L)[[method]]
  sintax_seed <<- seed
  set.seed(seed)
  elapsed <- system.time(
    out <- compute_assignment(
      input_pq,
      method = method,
      ref_fasta = ref,
      nproc = nproc
    )
  )[["elapsed"]]
  attr(out, "elapsed") <- elapsed
  out
}

# `input_pq` with every row of the computation (method x db) as its columns,
# named "<rank>_<full_name>" as in d_all_taxo_ITS1_asv.
derive_columns <- function(input_pq, computed, method, db, vm = pauvert_vm) {
  rows <- vm[vm$method == method & vm$db == db, ]
  derived <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]
    derive_assignment(
      input_pq,
      computed,
      method = method,
      suffix = paste0("_", r$full_name),
      min_bootstrap = r$min_bootstrap,
      vote_algorithm = r$vote_algorithm,
      nb_voting = r$nb_voting,
      min_cover = r$min_cover,
      min_id = r$min_id,
      lca_cutoff = r$lca_cutoff
    )
  })
  do.call(combine_taxo_assignments, c(list(input_pq), derived))
}

# Scores of the columns `cols` of `physeq` against `truth`, as chapter 01
# scores them (clean_taxo(), add_gen_sp(), tc_metrics_unit() per scoring
# group), with the method, database and parameters of each column.
score_columns <- function(physeq, cols, truth, vm = pauvert_vm) {
  d <- physeq |>
    clean_taxo() |>
    add_gen_sp()
  ranks_df <- build_ranks_df(d, cols, ranks = tax_order)
  score_assignments(
    d,
    ranks_df,
    truth,
    db_meta,
    external_truth = external_control_truth(taxa_names(d)),
    truth_ranks = truth_ranks
  ) |>
    mutate(tax_level = ifelse(tax_level == "Species", "Gen_sp", tax_level)) |>
    left_join(
      select(
        vm,
        full_name,
        method,
        db,
        min_bootstrap,
        lca_cutoff,
        vote_algorithm,
        min_id
      ),
      by = c("method_db" = "full_name")
    )
}

# The default-settings rows of score_columns() (bootstrap 0.5, lca 0.51,
# blastn rel_majority at 95 % identity), as filter_default_settings().
default_rows <- function(scores) {
  scores |>
    filter(
      (method %in% c("dada2", "sintax") & min_bootstrap == 0.5) |
        (method == "lca" & lca_cutoff == 0.51) |
        (method == "blastn" & vote_algorithm == "rel_majority" & min_id == 95)
    )
}

# The four co-primary endpoints of OD 37.
endpoint_rows <- function(scores) {
  scores |>
    filter(
      metrics %in% c("F1_score", "PPV"),
      tax_level %in% c("Genus", "Gen_sp")
    )
}
