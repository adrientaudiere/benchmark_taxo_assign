# Pins the target-name formula of R/values_map.R and checks the database
# configuration of config.R. Renaming a target silently orphans its store
# object, so any change to the formula must be deliberate. The tests use
# synthetic database names so that moving to a new release (new names in
# config.R) does not require editing them.
# Run from the project root:
#   Rscript tests/test_values_map.R

library("testthat")
library("here")
library("dplyr")

here::i_am("tests/test_values_map.R")
source(here("config.R"))
source(here("R/values_map.R"))

# The formula used by the pipelines before S2.1 (kept verbatim as the oracle).
legacy_full_name <- function(
  method,
  db,
  min_bootstrap,
  vote_algorithm,
  nb_voting
) {
  gsub(
    "...NA...NA",
    "",
    paste0(
      method,
      "__",
      db,
      "___",
      min_bootstrap,
      "...",
      vote_algorithm,
      "...",
      nb_voting
    )
  )
}

test_dbs <- c("DB_A", "DB_B_v2.1_Fungi_cut")

# 21 rows per database x input: dada2 3, sintax 3, lca 3, blastn 3 votes x 4
# min_id (ROADMAP 0.4).
rows_per_db <- 21

test_that("build_values_map keeps the legacy names of dada2 and sintax, and names lca and blastn rows by their parameters", {
  vm <- build_values_map(dbs = test_dbs, mini_db = FALSE, itsx_dbs = test_dbs)
  expect_equal(nrow(vm), rows_per_db * length(test_dbs) * 3)
  expect_false(any(duplicated(vm$full_name)))
  raw <- vm[vm$preprocess == "none", ]
  bootstrap_rows <- raw[raw$method %in% c("dada2", "sintax"), ]
  expect_equal(
    bootstrap_rows$full_name,
    legacy_full_name(
      bootstrap_rows$method,
      bootstrap_rows$db,
      bootstrap_rows$min_bootstrap,
      NA,
      NA
    )
  )
  blastn_rows <- raw[raw$method == "blastn", ]
  expect_equal(
    blastn_rows$full_name,
    paste0(
      legacy_full_name(
        "blastn",
        blastn_rows$db,
        0.5,
        blastn_rows$vote_algorithm,
        blastn_rows$nb_voting
      ),
      "...",
      blastn_rows$min_id
    )
  )
  expect_true(all(
    c(
      "dada2__DB_A___0.5",
      "sintax__DB_B_v2.1_Fungi_cut___0.4",
      "lca__DB_A___0.8",
      "lca__DB_A___1",
      "blastn__DB_B_v2.1_Fungi_cut___0.5...rel_majority...100...95",
      "blastn__DB_A___0.5...unanimity...100...90"
    ) %in%
      raw$full_name
  ))
  expect_setequal(
    raw$lca_cutoff[raw$method == "lca"],
    rep(c(0.8, 0.9, 1), length(test_dbs))
  )
  expect_true(all(is.na(raw$min_bootstrap[raw$method == "lca"])))
  expect_setequal(unique(blastn_rows$min_id), c(90, 92, 95, 97))
  expect_equal(
    nrow(build_values_map(dbs = db_list)),
    rows_per_db * (2 * length(db_list) + length(itsx_db_list))
  )
})

test_that("OTU rows run on every database, prefix the method and read d_vs_for_assignation", {
  vm <- build_values_map(dbs = test_dbs, itsx_dbs = "DB_A")
  otu <- vm[vm$preprocess == "otu", ]
  raw <- vm[vm$preprocess == "none", ]
  expect_setequal(unique(otu$db), test_dbs)
  expect_setequal(otu$full_name, paste0("otu_", raw$full_name))
  expect_equal(unique(otu$input_pq), "d_vs_for_assignation")
  expect_true("otu_compute_lca__DB_B_v2.1_Fungi_cut" %in% otu$compute_name)
})

test_that("the ITSx input runs on itsx_dbs only, the Fungi databases of config.R", {
  vm <- build_values_map(dbs = test_dbs, itsx_dbs = "DB_B_v2.1_Fungi_cut")
  expect_equal(unique(vm$db[vm$preprocess == "itsx"]), "DB_B_v2.1_Fungi_cut")
  expect_setequal(unique(vm$db[vm$preprocess == "none"]), test_dbs)
  expect_equal(sum(vm$preprocess == "itsx"), rows_per_db)
  expect_true(all(itsx_db_list %in% db_list))
  expect_setequal(
    itsx_db_list,
    benchmark_dbs$db[benchmark_dbs$simplification %in% c("Fungi", "Fungi+rep")]
  )
})

test_that("ITSx rows prefix the method and read the ITSx phyloseq", {
  vm <- build_values_map(dbs = test_dbs, itsx_dbs = test_dbs)
  itsx <- vm[vm$preprocess == "itsx", ]
  raw <- vm[vm$preprocess == "none", ]
  expect_equal(nrow(itsx), nrow(raw))
  expect_setequal(itsx$full_name, paste0("itsx_", raw$full_name))
  expect_true(
    "itsx_blastn__DB_A___0.5...rel_majority...100...95" %in% itsx$full_name
  )
  expect_equal(unique(itsx$input_pq), "d_asv_itsx_for_assignation")
  expect_equal(unique(raw$input_pq), "d_asv_for_assignation")
  expect_equal(
    nrow(build_values_map(dbs = test_dbs, preprocess = "none")),
    nrow(raw)
  )
})

test_that("db_path follows the method and the mini_db flag", {
  vm <- build_values_map(dbs = "DB_A", mini_db = FALSE)
  expect_equal(
    unique(vm$db_path[vm$method == "dada2"]),
    "data/data_raw/refseq/dada2_format/DB_A.fasta"
  )
  expect_equal(
    unique(vm$db_path[vm$method == "sintax"]),
    "data/data_raw/refseq/sintax_format/DB_A.fasta"
  )
  vm_mini <- build_values_map(dbs = "DB_A", mini_db = TRUE)
  expect_true(all(grepl("/mini_DB_A\\.fasta$", vm_mini$db_path)))
  expect_equal(vm_mini$full_name, vm$full_name)
})

test_that("every database gets its own reference file, in mini and full mode", {
  for (mini in c(FALSE, TRUE)) {
    prefix <- if (mini) "mini_" else ""
    vm <- build_values_map(dbs = db_list, mini_db = mini)
    expect_equal(
      dplyr::n_distinct(vm$db_path[vm$method == "sintax"]),
      length(db_list)
    )
    expect_equal(
      dplyr::n_distinct(vm$db_path[vm$method == "dada2"]),
      length(db_list)
    )
    expect_true(all(mapply(
      \(d, p) endsWith(p, paste0("/", prefix, d, ".fasta")),
      vm$db,
      vm$db_path
    )))
    cv <- build_cv_values_map(dbs = db_list, mini_db = mini)
    expect_true(all(mapply(
      \(d, p) endsWith(p, paste0("/", prefix, d, ".fasta")),
      cv$db,
      cv$db_path
    )))
  }
})

test_that("each reference file has one file target, named without the mini prefix", {
  for (mini in c(FALSE, TRUE)) {
    vm <- build_values_map(dbs = test_dbs, mini_db = mini)
    pairs <- dplyr::distinct(vm, ref_file, db_path)
    expect_equal(nrow(pairs), 2 * length(test_dbs))
    expect_false(any(duplicated(pairs$ref_file)))
    expect_false(any(duplicated(pairs$db_path)))
    cv <- build_cv_values_map(dbs = test_dbs, mini_db = mini)
    expect_setequal(
      dplyr::distinct(cv, ref_file, db_path)$ref_file,
      pairs$ref_file
    )
    expect_setequal(
      dplyr::distinct(cv, ref_file, db_path)$db_path,
      pairs$db_path
    )
  }
  vm_mini <- build_values_map(dbs = test_dbs, mini_db = TRUE)
  expect_true("ref_sintax__DB_B_v2.1_Fungi_cut" %in% vm_mini$ref_file)
  expect_equal(
    unique(vm_mini$ref_file[vm_mini$method == "dada2" & vm_mini$db == "DB_A"]),
    "ref_dada2__DB_A"
  )
  expect_false(any(vm_mini$ref_file %in% vm_mini$full_name))
})

test_that("one compute target per method × database × input serves all its rows", {
  vm <- build_values_map(dbs = test_dbs, itsx_dbs = test_dbs)
  counts <- table(vm$compute_name)
  expect_equal(length(counts), 4 * length(test_dbs) * 3)
  expect_true(all(counts[grepl("compute_blastn__", names(counts))] == 12))
  expect_true(all(counts[!grepl("compute_blastn__", names(counts))] == 3))
  expect_true(all(
    c(
      "compute_dada2__DB_A",
      "itsx_compute_blastn__DB_B_v2.1_Fungi_cut"
    ) %in%
      names(counts)
  ))
  expect_false(any(vm$compute_name %in% c(vm$full_name, vm$ref_file)))
  per_compute <- dplyr::distinct(
    vm,
    compute_name,
    method,
    db,
    preprocess,
    input_pq,
    ref_file,
    controller
  )
  expect_equal(nrow(per_compute), length(counts))
})

test_that("min_cover is set on blastn rows only, from config.R by default", {
  vm <- build_values_map(dbs = test_dbs)
  expect_true(all(vm$min_cover[vm$method == "blastn"] == blastn_min_cover))
  expect_true(all(is.na(vm$min_cover[vm$method != "blastn"])))
  vm_70 <- build_values_map(
    dbs = "DB_A",
    methods = build_methods_grid(min_cover = 70)
  )
  expect_equal(unique(vm_70$min_cover[vm_70$method == "blastn"]), 70)
  expect_equal(vm_70$full_name, build_values_map(dbs = "DB_A")$full_name)
  cv <- build_cv_values_map(dbs = test_dbs)
  expect_true(all(cv$min_cover[cv$method == "blastn"] == blastn_min_cover))
  expect_true(all(is.na(cv$min_cover[cv$method != "blastn"])))
})

test_that("build_cv_values_map reproduces the legacy CV target names", {
  cv <- build_cv_values_map(dbs = test_dbs, mini_db = FALSE)
  expect_equal(nrow(cv), 4 * length(test_dbs) * 2)
  expect_false(any(duplicated(cv$full_name)))
  expect_true(all(
    c(
      "cv__dada2__DB_A__standard",
      "cv__blastn__DB_B_v2.1_Fungi_cut__leaked"
    ) %in%
      cv$full_name
  ))
  expect_equal(unique(cv$vote_algorithm[cv$method == "blastn"]), "rel_majority")
  expect_true(all(is.na(cv$vote_algorithm[cv$method != "blastn"])))
})

test_that("CV queries of a _Fungi_cut database come from its _Fungi source", {
  cv <- build_cv_values_map(
    dbs = c("DB_A", "DB_B_v2.1_Fungi", "DB_B_v2.1_Fungi_cut")
  )
  cut <- cv[cv$db == "DB_B_v2.1_Fungi_cut", ]
  expect_true(all(cut$query_db == "DB_B_v2.1_Fungi"))
  expect_equal(
    unique(cut$query_ref_file[cut$method == "sintax"]),
    "ref_sintax__DB_B_v2.1_Fungi"
  )
  expect_equal(
    unique(cut$query_db_path[cut$method == "dada2"]),
    "data/data_raw/refseq/dada2_format/DB_B_v2.1_Fungi.fasta"
  )
  other <- cv[cv$db != "DB_B_v2.1_Fungi_cut", ]
  expect_equal(other$query_ref_file, other$ref_file)
  expect_equal(other$query_db_path, other$db_path)
  mini <- build_cv_values_map(dbs = "DB_B_v2.1_Fungi_cut", mini_db = TRUE)
  expect_true(all(grepl("/mini_DB_B_v2.1_Fungi\\.fasta$", mini$query_db_path)))
  cut_sources <- sub(
    "_cut$",
    "",
    zoom_dbs$db[zoom_dbs$simplification == "Fungi+cut"]
  )
  expect_true(all(cut_sources %in% cv_db_list))
})

test_that("CV records are named by the sintax_format identifiers of the same database", {
  cv <- build_cv_values_map(
    dbs = c("DB_A", "DB_B_v2.1_Fungi", "DB_B_v2.1_Fungi_cut")
  )
  expect_equal(
    cv$id_db_path,
    paste0("data/data_raw/refseq/sintax_format/", cv$db, ".fasta")
  )
  expect_equal(cv$id_ref_file, paste0("ref_sintax__", cv$db))
  expect_equal(
    cv$query_id_db_path,
    paste0("data/data_raw/refseq/sintax_format/", cv$query_db, ".fasta")
  )
  expect_equal(cv$query_id_ref_file, paste0("ref_sintax__", cv$query_db))
  sintax <- cv[cv$method == "sintax", ]
  expect_equal(sintax$id_db_path, sintax$db_path)
  mini <- build_cv_values_map(
    dbs = c("DB_A", "DB_B_v2.1_Fungi_cut"),
    mini_db = TRUE
  )
  expect_equal(
    mini$id_db_path,
    paste0("data/data_raw/refseq/sintax_format/mini_", mini$db, ".fasta")
  )
})

test_that("db_meta covers every database of db_list", {
  expect_setequal(db_meta$db, db_list)
})

test_that("config: sources, benchmarked databases and special picks are consistent", {
  expect_false(any(duplicated(reference_sources$source)))
  expect_true(all(reference_sources$provider %in% c("unite", "eukaryome")))
  expect_true(all(grepl("^https://", reference_sources$url)))
  expect_true(all(benchmark_dbs$source %in% reference_sources$source))
  expect_true(all(benchmark_dbs$simplification %in% names(db_suffix)))
  expect_true(all(zoom_dbs$simplification %in% names(db_suffix)))
  expect_false(any(duplicated(db_list)))
  expect_length(db_list, 9)
  expect_false(any(zoom_dbs$db %in% db_list))
  expect_true(all(cv_db_list %in% c(db_list, zoom_dbs$db)))
  expect_true(all(
    unique(benchmark_dbs$source[
      benchmark_dbs$simplification == "Fungi+rep"
    ]) %in%
      names(rep_kingdom_min_records)
  ))
  expect_true(all(bio_datasets$itsx_region %in% c("ITS1", "ITS2", "full")))
  expect_true(all(c(preference_db, seed_taxonomy_db) %in% db_list))
  expect_true(fake_ref_source %in% reference_sources$source)
  expect_true(
    reference_sources$provider[reference_sources$source == fake_ref_source] ==
      "unite"
  )
  expect_true(
    is.numeric(blastn_min_cover) &&
      blastn_min_cover > 0 &&
      blastn_min_cover <= 100
  )
})
