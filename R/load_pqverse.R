# Single entry point for the sister pqverse packages (ROADMAP S2.2, decision
# 10: the manuscript uses the development checkouts, not CRAN releases).
#
# Set PQVERSE_PKG_DIR in .Renviron to point at the directory that holds the
# MiscMetabar / comparpq / dbpq / greenAlgoR / tidypq checkouts. Defaults to the
# layout of the pqverse workspace.

pqverse_pkg_dir <- function() {
  normalizePath(
    Sys.getenv(
      "PQVERSE_PKG_DIR",
      unset = "~/Nextcloud/IdEst/Projets/pqverse/pqverse_pkg"
    ),
    mustWork = TRUE
  )
}

# Load the given packages from their checkouts with pkgload (what
# devtools::load_all() calls). Returns the root directory invisibly.
load_pqverse <- function(
  pkgs = c("MiscMetabar", "comparpq", "dbpq", "tidypq"),
  quiet = TRUE
) {
  root <- pqverse_pkg_dir()
  for (pkg in pkgs) {
    pkgload::load_all(file.path(root, pkg), quiet = quiet)
  }
  invisible(root)
}

# Output lines of `git -C path args`, or NA when git fails.
git_lines <- function(path, args) {
  tryCatch(
    system2("git", c("-C", shQuote(path), args), stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_,
    warning = function(w) NA_character_
  )
}

# Git state of one checkout: its commit, whether tracked files differ from it
# (`dirty`), and a hash of that difference (`diff_hash`, NA when clean), so
# that two runs on the same dirty commit can be told apart
# (reviews/code_review_2026-10-05.md §1.7). Untracked files do not count.
git_state <- function(path) {
  commit <- git_lines(path, c("rev-parse", "HEAD"))
  diff <- git_lines(path, c("diff", "HEAD"))
  known <- length(commit) == 1 && !is.na(commit)
  tibble::tibble(
    git_commit = if (known) commit else NA_character_,
    dirty = if (known) length(diff) > 0 else NA,
    diff_hash = if (known && length(diff) > 0) {
      rlang::hash(diff)
    } else {
      NA_character_
    }
  )
}

# Version and git state of each checkout, for the `session_info` target and
# the manuscript's Methods section.
pqverse_versions <- function(
  pkgs = c("MiscMetabar", "comparpq", "dbpq", "greenAlgoR", "tidypq")
) {
  root <- pqverse_pkg_dir()
  rows <- lapply(pkgs, function(pkg) {
    path <- file.path(root, pkg)
    dplyr::bind_cols(
      tibble::tibble(
        package = pkg,
        version = as.character(read.dcf(
          file.path(path, "DESCRIPTION"),
          fields = "Version"
        ))
      ),
      git_state(path)
    )
  })
  dplyr::bind_rows(rows)
}

# Git state of this project, for the same purpose (scoring happens in the
# Quarto chapters, which no session_info target covers).
project_version <- function(path = here::here()) {
  dplyr::bind_cols(
    tibble::tibble(package = "benchmar_assign_taxo"),
    git_state(path)
  )
}
