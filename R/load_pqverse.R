# Single entry point for the sister pqverse packages (ROADMAP S2.2, decision
# 10: the manuscript uses the development checkouts, not CRAN releases).
#
# Set PQVERSE_PKG_DIR in .Renviron to point at the directory that holds the
# MiscMetabar / comparpq / dbpq / greenAlgoR / tidypq checkouts. Defaults to the
# layout of the pqverse workspace.

pqverse_pkg_dir <- function() {
  normalizePath(
    Sys.getenv("PQVERSE_PKG_DIR",
               unset = "~/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg"),
    mustWork = TRUE
  )
}

# Load the given packages from their checkouts with pkgload (what
# devtools::load_all() calls). Returns the root directory invisibly.
load_pqverse <- function(pkgs = c("MiscMetabar", "comparpq", "dbpq", "tidypq"),
                         quiet = TRUE) {
  root <- pqverse_pkg_dir()
  for (pkg in pkgs) {
    pkgload::load_all(file.path(root, pkg), quiet = quiet)
  }
  invisible(root)
}

# Version and git commit of each checkout, for the `session_info` target and
# the manuscript's Methods section.
pqverse_versions <- function(pkgs = c("MiscMetabar", "comparpq", "dbpq",
                                      "greenAlgoR", "tidypq")) {
  root <- pqverse_pkg_dir()
  rows <- lapply(pkgs, function(pkg) {
    path <- file.path(root, pkg)
    commit <- tryCatch(
      system2("git", c("-C", shQuote(path), "rev-parse", "HEAD"),
              stdout = TRUE, stderr = FALSE),
      error = function(e) NA_character_,
      warning = function(w) NA_character_
    )
    tibble::tibble(
      package = pkg,
      version = as.character(read.dcf(file.path(path, "DESCRIPTION"),
                                      fields = "Version")),
      git_commit = if (length(commit) == 1) commit else NA_character_
    )
  })
  dplyr::bind_rows(rows)
}
