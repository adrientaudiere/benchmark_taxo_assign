# Converts single-min_bootstrap cross_val() output to a long analysis-ready tibble.
# Stops if given a vector-bootstrap result — pass a scalar min_bootstrap to cross_val().

cv_to_tidy <- function(cv_result, method, db, remove_tested, min_bootstrap = 0.5) {
  if ("bootstrap" %in% names(cv_result$good_classifications)) {
    stop(paste0(
      "cv_to_tidy() requires a single-min_bootstrap cross_val() result. ",
      "Pass a scalar (not a vector) min_bootstrap to cross_val()."
    ))
  }
  cv_result$metrics |>
    dplyr::rename(tax_level = name) |>
    dplyr::mutate(
      method        = method,
      db            = db,
      remove_tested = remove_tested,
      min_bootstrap = min_bootstrap
    )
}
