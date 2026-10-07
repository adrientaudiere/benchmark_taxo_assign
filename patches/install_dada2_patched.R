# Build dada2 1.40.0 with the overflow fix into r_lib_patched/ at the project
# root (git-ignored). Decision of the developer, 2026-09-30: this project uses
# the patched dada2 for every run (HISTORY.md §2, "the dada2 collapse on the
# full EUKARYOME is a 32-bit overflow in dada2").
#
# The fix (patches/dada2_1.40.0_taxonomy_overflow.patch, one line of
# src/taxonomy.cpp): get_best_genus() read the k-mer profile of taxonomy g at
# g * 65 536 computed in 32 bits, so every reference with more than 65 536
# distinct taxonomies had its taxonomies scored with the wrong profiles.
#
# .Rprofile puts r_lib_patched/ first on the library path of every R session
# started at the project root (crew workers included), and
# R/dada2_patched.R::assert_dada2_patched() stops a pipeline that would load
# another dada2. Run from the project root (about 2 min):
#   Rscript patches/install_dada2_patched.R

version <- "1.40.0"
url <- sprintf(
  "https://bioconductor.org/packages/3.23/bioc/src/contrib/dada2_%s.tar.gz",
  version
)
patch_file <- normalizePath("patches/dada2_1.40.0_taxonomy_overflow.patch")
lib <- file.path(getwd(), "r_lib_patched")
dir.create(lib, showWarnings = FALSE)

work <- tempfile("dada2_build_")
dir.create(work)
tarball <- file.path(work, basename(url))
utils::download.file(url, tarball, quiet = TRUE)
utils::untar(tarball, exdir = work)
src <- file.path(work, "dada2")

status <- system2(
  "patch",
  c("-p1", "-d", shQuote(src), "-i", shQuote(patch_file)),
  stdout = FALSE
)
if (status != 0) {
  stop("The patch did not apply to dada2 ", version, ".")
}
fixed <- readLines(file.path(src, "src/taxonomy.cpp"))
stopifnot(
  sum(grepl("lgk_probability[(size_t) g * n_kmers]", fixed, fixed = TRUE)) == 1
)

# Mark the build, so that assert_dada2_patched() can tell it from the
# Bioconductor one.
description <- file.path(src, "DESCRIPTION")
write(
  "Patched: taxonomy_overflow (benchmar_assign_taxo patches/dada2_1.40.0_taxonomy_overflow.patch)",
  description,
  append = TRUE
)

status <- system2(
  file.path(R.home("bin"), "R"),
  c(
    "CMD",
    "INSTALL",
    "--no-docs",
    "--no-html",
    "-l",
    shQuote(lib),
    shQuote(src)
  ),
  stdout = FALSE,
  stderr = FALSE
)
if (status != 0) {
  stop("R CMD INSTALL failed; rerun it by hand on ", src, " to see why.")
}
unlink(work, recursive = TRUE)
message(
  "dada2 ",
  utils::packageDescription("dada2", lib.loc = lib)$Version,
  " installed in ",
  lib,
  " (",
  utils::packageDescription("dada2", lib.loc = lib)$Patched,
  ")"
)
