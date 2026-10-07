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

# Git state of the code of one checkout: its commit, whether the files under
# `code` differ from it (`dirty`: a modified, deleted or untracked file), and
# a hash of that difference (`diff_hash`, NA when clean: the diff of the
# tracked files and the content of the untracked ones), so that two runs on
# the same dirty commit can be told apart (reviews/code_review_2026-10-05.md
# §1.7). Only `code` counts: the stores' tracked meta files change with every
# run and would make every run read as dirty.
git_state <- function(path, code) {
  commit <- git_lines(path, c("rev-parse", "HEAD"))
  known <- length(commit) == 1 && !is.na(commit)
  code <- code[file.exists(file.path(path, code))]
  if (!known || length(code) == 0) {
    return(tibble::tibble(
      git_commit = if (known) commit else NA_character_,
      dirty = NA,
      diff_hash = NA_character_
    ))
  }
  status <- git_lines(
    path,
    c("status", "--porcelain", "--untracked-files=all", "--", code)
  )
  dirty <- length(status) > 0
  diff_hash <- NA_character_
  if (dirty) {
    diff <- git_lines(path, c("diff", "HEAD", "--", code))
    untracked <- sub("^\\?\\? ", "", grep("^\\?\\? ", status, value = TRUE))
    contents <- vapply(
      file.path(path, untracked),
      \(f) paste(readLines(f, warn = FALSE), collapse = "\n"),
      character(1)
    )
    diff_hash <- rlang::hash(list(diff, untracked, unname(contents)))
  }
  tibble::tibble(git_commit = commit, dirty = dirty, diff_hash = diff_hash)
}

# What counts as code in a pqverse checkout and in this project.
package_code <- c("R", "src", "inst", "DESCRIPTION", "NAMESPACE")
project_code <- c(
  "R",
  "pipelines",
  "analysis",
  "tools",
  "tests",
  "config.R",
  "make_databases.R",
  "make.R",
  "run_bio_night.R",
  "_targets.yaml",
  ".Rprofile",
  "patches"
)

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
      git_state(path, package_code)
    )
  })
  dplyr::bind_rows(rows)
}

# Git state of this project, for the same purpose (scoring happens in the
# Quarto chapters, which no session_info target covers).
project_version <- function(path = here::here()) {
  dplyr::bind_cols(
    tibble::tibble(package = "benchmar_assign_taxo"),
    git_state(path, project_code)
  )
}
