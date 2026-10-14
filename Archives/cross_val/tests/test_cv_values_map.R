# The grid of the cross-validation (R/cv_values_map.R, moved out of
# R/values_map.R on 2026-10-14). Like every file of Archives/cross_val/, it
# keeps its project-relative place: restored, it runs from the project root as
#   Rscript tests/test_cv_values_map.R

library("testthat")
library("here")
library("dplyr")

here::i_am("tests/test_cv_values_map.R")
source(here("config.R"))
source(here("config_cv.R"))
source(here("R/values_map.R"))
source(here("R/cv_values_map.R"))

test_dbs <- c("DB_A", "DB_B_v2.1_Fungi_cut")

test_that("the CV grid reads the same reference files as the main grid", {
  for (mini in c(FALSE, TRUE)) {
    prefix <- if (mini) "mini_" else ""
    cv <- build_cv_values_map(dbs = db_list, mini_db = mini)
    expect_true(all(mapply(
      \(d, p) endsWith(p, paste0("/", prefix, d, ".fasta")),
      cv$db,
      cv$db_path
    )))
    pairs <- dplyr::distinct(
      build_values_map(dbs = test_dbs, mini_db = mini),
      ref_file,
      db_path
    )
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
  cv <- build_cv_values_map(dbs = test_dbs)
  expect_true(all(cv$min_cover[cv$method == "blastn"] == blastn_min_cover))
  expect_true(all(is.na(cv$min_cover[cv$method != "blastn"])))
  expect_true(all(cv_db_list %in% c(db_list, zoom_dbs$db)))
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
