### Draws how a mock community is scored, and renders docs/mock_scoring.svg
### with Graphviz `dot` (on PATH). Run from the project root:
###
###   Rscript tools/mock_scoring_schema.R
###
### Three tables, one CSV each:
### - 1a, docs/mock_scoring_matrix.csv: the truth matrix, the cell a unit
###   falls in at one rank for one assignment column. A cell starts with TP,
###   FP, FN, TN, "out" (left out of the matrix at that rank) or "—" (no such
###   case), and may add a note after a colon;
### - 1b, docs/mock_scoring_metrics.csv: the metrics computed from the cells;
### - 2, docs/mock_scoring_table.csv: the extended mock table, one row per
###   kind of unit (real units by truth status, then the controls), one
###   column per column of the phyloseq after R/itsx_inputs.R::with_controls()
###   and of the truth joined to it. Equal neighbouring cells of a column are
###   merged within a group.
### The rules are those of comparpq::unit_cells()
### (comparpq/R/compare_taxo_unit.R) and R/score_assignments.R; edit the CSV
### when a decision changes them. `source` names the decision behind a row.

library("here")
here::i_am("tools/mock_scoring_schema.R")

matrix_path <- here("docs/mock_scoring_matrix.csv")
metrics_path <- here("docs/mock_scoring_metrics.csv")
table_path <- here("docs/mock_scoring_table.csv")
svg_path <- here("docs/mock_scoring.svg")
wrap_width <- 24

cell_colors <- c(
  TP = "#A1D99B",
  FP = "#FC9272",
  FN = "#FDD0A2",
  TN = "#9ECAE1",
  out = "#E0E0E0",
  "—" = "#FFFFFF"
)
group_colors <- c(
  "real units" = "#EFF3FF",
  "fake_" = "#FEE6CE",
  "external_" = "#E5F5E0"
)
header_color <- "#D9D9D9"
# Column headers of the matrix (1a), in the order of its CSV columns.
matrix_headers <- c(
  group = "Unit",
  unit = "",
  when = "When",
  named_truth = "Named its truth at the rank (or the accepted name, Genus and Species)",
  named_other = "Named another taxon (any name when the unit has no truth)",
  left_na = "Left NA",
  reported_by = "Also counted by",
  source = "Source"
)
matrix_notes <- c(
  "One cell per unit, per assignment column and per rank: TP + FP + FN + TN = the units scored at that rank.",
  "OD 49: the headline PPV, FDR and F1 count the real units only (FP_real); the FP and TN of the controls enter TNR, MCC, ACC and the *_with_ctrl values.",
  "OD 44: every mock is scored under both truth rules; no_truth_aside leads, q2b is the stricter reading beside it.",
  "A consensus column has no database: it is scored with the external_ controls in the matrix and aside."
)
# Spanning headers of the extended table (2): a title and the number of CSV
# columns it covers, left to right.
table_groups <- c(
  "Unit" = 2,
  "phyloseq of the input, after with_controls() (R/itsx_inputs.R)" = 5,
  "Truth of the mock, joined by taxon: mock_truth_table() (R/mock_truth.R)" = 4,
  "external_control_truth()" = 1,
  "comparpq::tc_cells_unit()" = 1
)
table_notes <- c(
  "no_truth and foreign are not prefixes: they come from matching each real unit against the Sanger reads of the strains of its mock (mock_truth_table()).",
  "A unit gets the truth of a strain at identity ≥ 99.5 % with the shorter sequence covered at ≥ 99.5 %, or an overlap ≥ 100 bp at ≥ 99.5 % (OD 36); with no hit at 90 % or more it is foreign.",
  "The controls are not in the truth table, and get 0 reads in every sample: every control metric counts units, not reads."
)

read_table <- function(path) {
  read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = character(0),
    encoding = "UTF-8"
  )
}
matrix_rows <- read_table(matrix_path)
metric_rows <- read_table(metrics_path)
table_rows <- read_table(table_path)

# ---- Checks -----------------------------------------------------------------
if (!identical(names(matrix_rows), names(matrix_headers))) {
  stop(
    matrix_path,
    " must have the columns ",
    toString(names(matrix_headers))
  )
}
metric_columns <- c("metric", "formula", "counted_on", "source")
if (!identical(names(metric_rows), metric_columns)) {
  stop(metrics_path, " must have the columns ", toString(metric_columns))
}
if (sum(table_groups) != ncol(table_rows)) {
  stop(
    "`table_groups` covers ",
    sum(table_groups),
    " columns, ",
    table_path,
    " has ",
    ncol(table_rows)
  )
}
unknown_groups <- setdiff(
  c(matrix_rows$group, table_rows$group),
  names(group_colors)
)
if (length(unknown_groups) > 0) {
  stop(
    "Unknown groups (add them to `group_colors`): ",
    toString(unknown_groups)
  )
}
cell_code <- function(x) {
  trimws(sub(":.*$", "", x))
}
cells <- unlist(matrix_rows[c("named_truth", "named_other", "left_na")])
unknown_cells <- setdiff(cell_code(cells), names(cell_colors))
if (length(unknown_cells) > 0) {
  stop(
    "A matrix cell must start with one of ",
    toString(names(cell_colors)),
    "; found ",
    toString(unknown_cells)
  )
}

# ---- HTML helpers -----------------------------------------------------------
# Text for a Graphviz HTML label: wrapped, then escaped line by line.
html_text <- function(x, width = wrap_width) {
  lines <- strwrap(x, width = width)
  if (length(lines) == 0) {
    return("")
  }
  lines <- gsub("&", "&amp;", lines, fixed = TRUE)
  lines <- gsub("<", "&lt;", lines, fixed = TRUE)
  lines <- gsub(">", "&gt;", lines, fixed = TRUE)
  paste(lines, collapse = "<BR ALIGN=\"LEFT\"/>")
}

td <- function(
  text,
  bg = "#FFFFFF",
  colspan = 1,
  rowspan = 1,
  bold = FALSE,
  size = 11,
  color = "black",
  width = wrap_width,
  border = 1
) {
  body <- html_text(text, width = width)
  if (bold && nzchar(body)) {
    body <- paste0("<B>", body, "</B>")
  }
  # Graphviz refuses an empty <FONT>: an empty cell gets none.
  if (nzchar(body)) {
    body <- sprintf(
      "<FONT POINT-SIZE=\"%d\" COLOR=\"%s\">%s</FONT>",
      as.integer(size),
      color,
      body
    )
  }
  sprintf(
    "<TD BGCOLOR=\"%s\" COLSPAN=\"%d\" ROWSPAN=\"%d\" BORDER=\"%d\" ALIGN=\"LEFT\" BALIGN=\"LEFT\" VALIGN=\"MIDDLE\">%s</TD>",
    bg,
    as.integer(colspan),
    as.integer(rowspan),
    as.integer(border),
    body
  )
}

# A title row and note rows spanning the whole table, without borders.
title_row <- function(text, ncol) {
  sprintf(
    "<TR>%s</TR>",
    td(text, colspan = ncol, bold = TRUE, size = 16, width = 200, border = 0)
  )
}
note_rows <- function(notes, ncol) {
  vapply(
    notes,
    function(n) {
      sprintf(
        "<TR>%s</TR>",
        td(n, colspan = ncol, size = 10, width = 160, border = 0)
      )
    },
    character(1),
    USE.NAMES = FALSE
  )
}

# The row span of each cell of `x` when equal neighbours sharing `key` are
# merged: the run length on the first cell of a run, 0 on the others.
run_spans <- function(x, key) {
  runs <- rle(paste(key, x, sep = "\r"))
  spans <- integer(length(x))
  spans[cumsum(runs$lengths) - runs$lengths + 1] <- runs$lengths
  spans
}

html_table <- function(rows) {
  paste0(
    "<<TABLE BORDER=\"0\" CELLBORDER=\"1\" CELLSPACING=\"0\" CELLPADDING=\"5\" COLOR=\"#9E9E9E\">",
    paste(rows, collapse = ""),
    "</TABLE>>"
  )
}

# ---- 1a. Truth matrix -------------------------------------------------------
matrix_cell <- function(value) {
  code <- cell_code(value)
  note <- if (grepl(":", value)) trimws(sub("^[^:]*:", "", value)) else ""
  body <- paste0("<B>", html_text(code), "</B>")
  if (nzchar(note)) {
    body <- paste0(
      body,
      "<BR ALIGN=\"LEFT\"/><FONT POINT-SIZE=\"9\">",
      html_text(note, width = 20),
      "</FONT>"
    )
  }
  sprintf(
    "<TD BGCOLOR=\"%s\" ALIGN=\"LEFT\" BALIGN=\"LEFT\"><FONT POINT-SIZE=\"11\">%s</FONT></TD>",
    cell_colors[[code]],
    body
  )
}

matrix_body <- function(m) {
  group_span <- run_spans(m$group, "")
  unit_span <- run_spans(m$unit, m$group)
  vapply(
    seq_len(nrow(m)),
    function(i) {
      bg <- group_colors[[m$group[i]]]
      paste0(
        "<TR>",
        if (group_span[i] > 0) {
          td(
            m$group[i],
            bg = bg,
            rowspan = group_span[i],
            bold = TRUE,
            width = 10
          )
        },
        if (unit_span[i] > 0) {
          td(m$unit[i], bg = bg, rowspan = unit_span[i])
        },
        td(m$when[i], bg = bg, size = 10),
        matrix_cell(m$named_truth[i]),
        matrix_cell(m$named_other[i]),
        matrix_cell(m$left_na[i]),
        td(m$reported_by[i], size = 10),
        td(m$source[i], size = 9, color = "#616161", width = 14),
        "</TR>"
      )
    },
    character(1)
  )
}

matrix_label <- html_table(c(
  title_row(
    "1a. Truth matrix: the cell of a unit at one rank, for one assignment column (comparpq::unit_cells())",
    ncol(matrix_rows)
  ),
  paste0(
    "<TR>",
    td(matrix_headers[["group"]], bg = header_color, colspan = 2, bold = TRUE),
    paste(
      vapply(
        matrix_headers[-(1:2)],
        \(h) td(h, bg = header_color, bold = TRUE),
        character(1)
      ),
      collapse = ""
    ),
    "</TR>"
  ),
  matrix_body(matrix_rows),
  note_rows(matrix_notes, ncol(matrix_rows))
))

# ---- 1b. Metrics ------------------------------------------------------------
metrics_label <- html_table(c(
  title_row(
    "1b. Metrics computed from the cells, per assignment column and per rank",
    ncol(metric_rows)
  ),
  paste0(
    "<TR>",
    paste(
      vapply(
        c("Metric", "Formula", "Counted on", "Source"),
        \(h) td(h, bg = header_color, bold = TRUE),
        character(1)
      ),
      collapse = ""
    ),
    "</TR>"
  ),
  vapply(
    seq_len(nrow(metric_rows)),
    function(i) {
      paste0(
        "<TR>",
        td(metric_rows$metric[i], bold = TRUE, width = 22),
        td(metric_rows$formula[i], width = 30),
        td(metric_rows$counted_on[i], size = 10, width = 32),
        td(metric_rows$source[i], size = 9, color = "#616161", width = 14),
        "</TR>"
      )
    },
    character(1)
  )
))

# ---- 2. Extended mock table -------------------------------------------------
table_body <- function(t) {
  spans <- lapply(t, \(column) run_spans(column, t$group))
  vapply(
    seq_len(nrow(t)),
    function(i) {
      bg <- group_colors[[t$group[i]]]
      cells <- vapply(
        seq_along(t),
        function(j) {
          if (spans[[j]][i] == 0) {
            return("")
          }
          td(
            t[[j]][i],
            bg = bg,
            rowspan = spans[[j]][i],
            bold = j <= 2,
            size = if (j <= 2) 11 else 10,
            width = if (j == 1) 10 else 20
          )
        },
        character(1)
      )
      paste0("<TR>", paste(cells, collapse = ""), "</TR>")
    },
    character(1)
  )
}

table_label <- html_table(c(
  title_row(
    "2. Extended mock table: the units of one mock input, its controls and their truth",
    ncol(table_rows)
  ),
  paste0(
    "<TR>",
    paste(
      mapply(
        function(title, n) {
          td(title, bg = "#BDBDBD", colspan = n, bold = TRUE, width = 30 * n)
        },
        names(table_groups),
        table_groups
      ),
      collapse = ""
    ),
    "</TR>"
  ),
  paste0(
    "<TR>",
    td("", bg = header_color, colspan = 2),
    paste(
      vapply(
        names(table_rows)[-(1:2)],
        \(h) td(h, bg = header_color, bold = TRUE, width = 20),
        character(1)
      ),
      collapse = ""
    ),
    "</TR>"
  ),
  table_body(table_rows),
  note_rows(table_notes, ncol(table_rows))
))

# ---- DOT and SVG ------------------------------------------------------------
dot <- c(
  "digraph mock_scoring {",
  "  graph [rankdir=TB, nodesep=0.5, ranksep=0.5, fontname=\"Helvetica\"];",
  "  node [shape=plaintext, fontname=\"Helvetica\"];",
  paste0("  truth_matrix [label=", matrix_label, "];"),
  paste0("  metrics [label=", metrics_label, "];"),
  paste0("  extended_table [label=", table_label, "];"),
  "  { rank=same; truth_matrix; metrics; }",
  "  truth_matrix -> extended_table [style=invis];",
  "}"
)

if (!nzchar(Sys.which("dot"))) {
  stop("Graphviz `dot` is not on PATH.")
}
dot_path <- tempfile(fileext = ".dot")
writeLines(dot, dot_path, useBytes = TRUE)
out <- suppressWarnings(system2(
  "dot",
  c("-Tsvg", shQuote(dot_path), "-o", shQuote(svg_path)),
  stdout = TRUE,
  stderr = TRUE
))
status <- attr(out, "status")
if (!is.null(status) && status != 0) {
  stop("dot failed:\n", paste(out, collapse = "\n"))
}

message(
  "mock_scoring_schema: ",
  nrow(matrix_rows),
  " matrix rows, ",
  nrow(metric_rows),
  " metrics, ",
  nrow(table_rows),
  " kinds of unit -> ",
  svg_path
)
