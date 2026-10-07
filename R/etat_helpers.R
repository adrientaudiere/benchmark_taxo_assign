# Dated figures of the project for ETAT.md, written by tools/etat.R
# (code review §6, §11 Q7). The other documents state decisions and point to
# ETAT.md for counts, so the counts no longer go stale in prose.
#
# Each run puts today's block in front of the previous ones; a second run on
# the same day replaces that day's block.

# One row per targets project of `_targets.yaml`: targets in the pipeline,
# built in the store, outdated, errored, orphaned (in the store but no longer
# in the pipeline), and for the assignment projects the number of
# computations (`*_compute_*` targets) and of assignment rows (full_name
# targets, which carry "___").
project_counts <- function(project, store) {
  old <- Sys.getenv("TAR_PROJECT", unset = NA)
  on.exit(
    if (is.na(old)) {
      Sys.unsetenv("TAR_PROJECT")
    } else {
      Sys.setenv(TAR_PROJECT = old)
    }
  )
  Sys.setenv(TAR_PROJECT = project)
  manifest <- targets::tar_manifest(fields = "name")$name
  meta <- if (dir.exists(store)) {
    m <- targets::tar_meta(store = store, fields = c("name", "type", "error"))
    m[!m$type %in% c("function", "object"), ]
  } else {
    data.frame(name = character(), type = character(), error = character())
  }
  current <- meta[meta$name %in% manifest, ]
  tibble::tibble(
    project = project,
    targets = length(manifest),
    built = nrow(current),
    outdated = length(targets::tar_outdated(reporter = "silent")),
    errored = sum(!is.na(current$error)),
    orphans = sum(!meta$name %in% manifest),
    computations = sum(grepl("_compute_", manifest)),
    assignment_rows = sum(grepl("___", manifest))
  )
}

# Units of each assignment input of a store: the real units
# (`d_<input>_real`) and the controls added to them
# (`d_<input>_for_assignation` minus the real units).
input_units <- function(project, store) {
  if (!dir.exists(store)) {
    return(NULL)
  }
  names <- targets::tar_meta(store = store, fields = "name")$name
  real <- grep("^d_.+_real$", names, value = TRUE)
  if (length(real) == 0) {
    return(NULL)
  }
  ntaxa_of <- function(name) {
    if (!name %in% names) {
      return(NA_integer_)
    }
    phyloseq::ntaxa(targets::tar_read_raw(name, store = store))
  }
  input <- sub("^d_(.+)_real$", "\\1", real)
  n_real <- vapply(real, ntaxa_of, integer(1))
  n_all <- vapply(
    paste0("d_", input, "_for_assignation"),
    ntaxa_of,
    integer(1)
  )
  tibble::tibble(
    project = project,
    input = input,
    real_units = unname(n_real),
    controls = unname(n_all - n_real)
  )
}

# Test files and expectations, counted in the sources (not by running them).
test_counts <- function(dir) {
  files <- list.files(dir, pattern = "^test_.*\\.R$", full.names = TRUE)
  code <- unlist(lapply(files, readLines))
  c(
    files = length(files),
    tests = sum(grepl("test_that\\(", code)),
    expectations = sum(grepl("\\bexpect_[a-z_]+\\(", code))
  )
}

# The markdown block of one date.
etat_block <- function(date, projects, units, tests) {
  c(
    paste0("## ", date),
    "",
    "### Targets projects",
    "",
    knitr::kable(projects, format = "pipe"),
    "",
    "### Assignment inputs (units)",
    "",
    if (is.null(units) || nrow(units) == 0) {
      "No input built."
    } else {
      knitr::kable(units, format = "pipe")
    },
    "",
    "### Tests",
    "",
    paste0(
      tests[["files"]],
      " test files, ",
      tests[["tests"]],
      " `test_that()` blocks, ",
      tests[["expectations"]],
      " `expect_*()` calls (counted in the sources)."
    ),
    ""
  )
}

# `lines` of ETAT.md with `block` in front of the dated blocks, after the
# header; an existing block of the same date is replaced.
prepend_etat <- function(lines, block, date) {
  starts <- grep("^## [0-9]{4}-[0-9]{2}-[0-9]{2}$", lines)
  header <- if (length(starts) == 0) lines else lines[seq_len(starts[1] - 1)]
  dated <- if (length(starts) == 0) {
    character()
  } else {
    lines[starts[1]:length(lines)]
  }
  if (length(dated) > 0 && dated[1] == paste0("## ", date)) {
    next_start <- grep("^## [0-9]{4}-[0-9]{2}-[0-9]{2}$", dated)[-1]
    dated <- if (length(next_start) == 0) {
      character()
    } else {
      dated[next_start[1]:length(dated)]
    }
  }
  c(header, block, dated)
}
