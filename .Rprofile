# dada2 1.40.0 patched for the overflow of assignTaxonomy() above 65 536
# distinct taxonomies (developer, 2026-09-30; HISTORY.md §2). Built by
# patches/install_dada2_patched.R into r_lib_patched/; put first on the
# library path of every R session started here, crew workers included.
if (dir.exists("r_lib_patched")) {
  .libPaths(c(normalizePath("r_lib_patched"), .libPaths()))
}
