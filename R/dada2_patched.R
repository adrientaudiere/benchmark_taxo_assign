# The project runs dada2 1.40.0 patched for the overflow of assignTaxonomy()
# above 65 536 distinct taxonomies (developer, 2026-09-30; HISTORY.md §2). The
# build lives in r_lib_patched/ (patches/install_dada2_patched.R) and
# .Rprofile puts it first on the library path of every R session started at
# the project root. Called at the top of the pipeline scripts, outside any
# target, so that it guards the runs without being a target dependency.
assert_dada2_patched <- function() {
  patched <- utils::packageDescription("dada2")$Patched
  if (is.null(patched) || !grepl("taxonomy_overflow", patched)) {
    stop(
      "dada2 is not the patched build (loaded from ",
      find.package("dada2"),
      "). Run `Rscript patches/install_dada2_patched.R` from the project ",
      "root, and start R from the project root so that .Rprofile applies."
    )
  }
  invisible(TRUE)
}
