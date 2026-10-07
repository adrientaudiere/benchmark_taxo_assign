# Species by exact and near match (OD 43): the raw hits of the species
# section, for every assignment input of every dataset x every database.
# Scoring is done afterwards; this project only computes and stores
# (R/species_exact.R).
#
#   queries: the inputs d_<region>_<unit>_for_assignation of every
#            assignment store, pooled into distinct sequences (q1..qN) with
#            the map back to (dataset, input, taxon); queries with a base other
#            than A, C, G or T are left out for the three methods
#   per database, one call on all pooled queries:
#     species_assign_<db>   dada2::assignSpecies() on species_format/<db>.fasta
#     species_vsearch_<db>  vsearch --usearch_global --id 1 --query_cov 1
#                           --maxaccepts 0 --maxrejects 32 on sintax_format/
#     species_vsearch_near_<db>  vsearch candidates at --id 0.9 --query_cov
#                           0.9, kept with at most one difference over the
#                           query, or 100 % over >= 90 % of it with the record
#                           stopping early (near matches,
#                           developer 2026-10-07)
#     species_blastn_<db>   blastn's own search (megablast, >= 90 % identity
#                           over >= 89 % of the query, 2 000 records per
#                           query), same pairs kept; its exact match is the
#                           pairs with no difference. The grid's blastn tables
#                           were used before 2026-10-07 (500 hits per query,
#                           no alignment coordinates).
#
# Inputs are read from the stores of the assignment projects, never written.
# A missing or errored input object stops the first target with the list of
# what is missing, so a failed assignment run upstream (e.g. the same night)
# never gives a partial result.
#
# Run with:
#   Sys.setenv(TAR_PROJECT = "species"); targets::tar_make()
# Smoke test: TAR_PROJECT = "species_mini" (Pauvert and Tedersoo against
# Unite_all_20250219_Fungi, own store). SPECIES_DATASETS = "pauvert,..."
# restricts the datasets of either project.

library("conflicted")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("dplyr")

here::i_am("pipelines/species.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/autometric_helpers.R"))
source(here("R/species_exact.R"))
source(here("R/dada2_patched.R"))
assert_dada2_patched()
load_pqverse("dbpq")

species_project <- Sys.getenv("TAR_PROJECT", unset = "species")
species_mini <- grepl("_mini$", species_project)

tar_option_set(
  seed = targets_seed,
  packages = c("here", "tibble", "dplyr", "phyloseq", "Biostrings")
)

# Datasets and their assignment stores: the Pauvert mock of the production
# project, then every row of config.R::bio_datasets.
species_datasets <- dplyr::bind_rows(
  tibble(
    dataset = "pauvert",
    store = "store_assign_taxo",
    itsx_region = "ITS1",
    otu_input = TRUE
  ),
  dplyr::transmute(
    bio_datasets,
    dataset,
    store = paste0("store_assign_taxo_", dataset),
    itsx_region,
    otu_input
  )
)
keep_datasets <- Sys.getenv(
  "SPECIES_DATASETS",
  unset = if (species_mini) "pauvert,tedersoo_illumina" else ""
)
if (nzchar(keep_datasets)) {
  keep_datasets <- strsplit(keep_datasets, ",", fixed = TRUE)[[1]]
  stopifnot(all(keep_datasets %in% species_datasets$dataset))
  species_datasets <- species_datasets[
    species_datasets$dataset %in% keep_datasets,
  ]
}
species_dbs <- if (species_mini) "Unite_all_20250219_Fungi" else db_list

# One row per input: dataset x ITSx region x unit.
species_inputs <- dplyr::bind_rows(lapply(
  seq_len(nrow(species_datasets)),
  function(i) {
    d <- species_datasets[i, ]
    tidyr::expand_grid(
      dataset = d$dataset,
      store = d$store,
      region = strsplit(d$itsx_region, ",", fixed = TRUE)[[1]],
      unit = if (d$otu_input) c("asv", "otu") else "asv"
    )
  }
)) |>
  dplyr::mutate(
    input = paste0(region, "_", unit),
    object = paste0("d_", input, "_for_assignation")
  )

# Every store object the project reads: the inputs. Stops with the list of the
# objects that are missing or errored.
species_store_objects <- function(inputs = species_inputs) {
  wanted <- dplyr::transmute(inputs, store, name = object)
  ok <- vapply(
    seq_len(nrow(wanted)),
    function(i) {
      store <- here(wanted$store[i])
      if (!dir.exists(store)) {
        return(FALSE)
      }
      meta <- targets::tar_meta(
        names = tidyselect::any_of(wanted$name[i]),
        fields = c("name", "error"),
        store = store
      )
      nrow(meta) == 1 &&
        (is.na(meta$error) || !nzchar(meta$error)) &&
        file.exists(file.path(store, "objects", wanted$name[i]))
    },
    logical(1)
  )
  if (!all(ok)) {
    missing <- paste0(wanted$store[!ok], "/", wanted$name[!ok])
    stop(
      "species: ",
      sum(!ok),
      " input objects are missing or errored (run their assignment ",
      "project first): ",
      paste(utils::head(missing, 10), collapse = ", "),
      if (sum(!ok) > 10) ", ..."
    )
  }
  file.path(here(wanted$store), "objects", wanted$name)
}

autometric_dir_species <- here("data/data_final/autometric", species_project)

db_map <- tibble(
  db = species_dbs,
  species_path = here(
    "data/data_raw/refseq/species_format",
    paste0(db, ".fasta")
  ),
  sintax_path = here("data/data_raw/refseq/sintax_format", paste0(db, ".fasta"))
)

per_db_targets <- tarchetypes::tar_map(
  values = db_map,
  names = "db",
  unlist = FALSE,
  tar_target(species_ref, species_path, format = "file"),
  tar_target(sintax_ref, sintax_path, format = "file"),
  tar_target(
    species_assign,
    with_autometric(
      paste0("species_assign__", db),
      assign_species_exact(species_queries$seqs, species_ref),
      dir = autometric_dir_species
    ) |>
      dplyr::mutate(db = db, .before = 1)
  ),
  tar_target(
    species_vsearch,
    with_autometric(
      paste0("species_vsearch__", db),
      # Every core, as vsearch in the main grid (config.R::assign_threads);
      # assignSpecies() has no thread option.
      vsearch_exact_hits(
        species_queries$seqs,
        sintax_ref,
        nproc = assign_threads[["lca"]]
      ),
      dir = autometric_dir_species
    ) |>
      dplyr::mutate(db = db, .before = 1)
  ),
  tar_target(
    species_vsearch_near,
    with_autometric(
      paste0("species_vsearch_near__", db),
      vsearch_near_hits(
        species_queries$seqs,
        sintax_ref,
        nproc = assign_threads[["lca"]]
      ),
      dir = autometric_dir_species
    ) |>
      dplyr::mutate(db = db, .before = 1)
  ),
  tar_target(
    species_blastn,
    with_autometric(
      paste0("species_blastn__", db),
      blastn_near_hits(
        species_queries$seqs,
        sintax_ref,
        nproc = assign_threads[["blastn"]]
      ),
      dir = autometric_dir_species
    ) |>
      dplyr::mutate(db = db, .before = 1)
  )
)

list(
  tar_target(
    species_input_files,
    species_store_objects(),
    format = "file"
  ),
  tar_target(species_query_inputs, {
    species_input_files
    stats::setNames(
      lapply(seq_len(nrow(species_inputs)), function(i) {
        targets::tar_read_raw(
          species_inputs$object[i],
          store = here(species_inputs$store[i])
        )@refseq
      }),
      paste0(species_inputs$dataset, "|", species_inputs$input)
    )
  }),
  tar_target(species_queries, pool_queries(species_query_inputs)),
  per_db_targets,
  tarchetypes::tar_combine(
    species_assign_all,
    per_db_targets[["species_assign"]]
  ),
  tarchetypes::tar_combine(
    species_vsearch_all,
    per_db_targets[["species_vsearch"]]
  ),
  tarchetypes::tar_combine(
    species_vsearch_near_all,
    per_db_targets[["species_vsearch_near"]]
  ),
  tarchetypes::tar_combine(
    species_blastn_all,
    per_db_targets[["species_blastn"]]
  ),
  tar_target(session_info, {
    species_blastn_all
    species_vsearch_near_all
    list(
      pqverse = pqverse_versions(),
      project = project_version(),
      session = sessioninfo::session_info()
    )
  })
)
