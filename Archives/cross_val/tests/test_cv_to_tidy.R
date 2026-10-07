# cv_to_tidy() turns one cross_val() result into the long tibble bound into
# cv_results. Run from the project root:
#   Rscript tests/test_cv_to_tidy.R

library("testthat")
library("here")

here::i_am("tests/test_cv_to_tidy.R")
source(here("R/cv_to_tidy.R"))

fake_cv <- list(
  metrics = tibble::tibble(
    name   = c("Kingdom", "Genus", "Kingdom", "Genus"),
    metric = c("good_classifications", "good_classifications", "NA_classif", "NA_classif"),
    mean   = c(0.95, 0.60, 0.02, 0.30),
    sd     = c(0.01, 0.05, 0.01, 0.04)
  ),
  good_classifications = data.frame(Kingdom = 0.95, Genus = 0.60)
)

test_that("cv_to_tidy adds the run metadata and renames name -> tax_level", {
  out <- cv_to_tidy(fake_cv, method = "sintax", db = "Unite",
                    remove_tested = TRUE, min_bootstrap = 0.5)
  expect_equal(nrow(out), 4)
  expect_true(all(c("tax_level", "metric", "mean", "sd",
                    "method", "db", "remove_tested", "min_bootstrap") %in% names(out)))
  expect_false("name" %in% names(out))
  expect_equal(unique(out$method), "sintax")
  expect_equal(unique(out$db), "Unite")
  expect_true(all(out$remove_tested))
  expect_equal(unique(out$min_bootstrap), 0.5)
})

test_that("cv_to_tidy refuses a vector-bootstrap cross_val() result", {
  multi <- fake_cv
  multi$good_classifications <- data.frame(bootstrap = c(0.4, 0.5), Kingdom = 1, Genus = 1)
  expect_error(cv_to_tidy(multi, "sintax", "Unite", TRUE), "single-min_bootstrap")
})
