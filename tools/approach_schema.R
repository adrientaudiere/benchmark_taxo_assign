### Draws the approach of the benchmark as a graph: reads the relations of
### docs/approach_schema.csv and renders, with Graphviz `dot` (on PATH),
### docs/approach_schema.svg (every row) and one docs/approach_schema_use<k>.svg
### per scientific question (the rows whose `uses` names it, without the
### types of `hidden_types`). Run from the project root:
###
###   Rscript tools/approach_schema.R
###
### One row of the CSV is one relation: `from` (of type `from_type`)
### `relation` `to` (of type `to_type`). A node is named by its label, so a
### label changed in one row must be changed in every row. `status` is
### "built" (read in a chapter, solid line), "to compute" (an open ROADMAP
### item, dashed line) or "parked" (docs/parked.md, grey dotted line); a node
### whose relations all share a status other than "built" gets that border
### too. `uses` lists the scientific questions a row is drawn for: "all", or
### their numbers separated by spaces ("1 2"); a scientific question is
### labelled "Use <k>: ...", and its own rows name it alone. `source` names
### the decision or chapter behind the relation and is not drawn. A new type
### needs a row in `node_types` below, a new status a row in `statuses`.

library("here")
here::i_am("tools/approach_schema.R")

csv_path <- here("docs/approach_schema.csv")
svg_path <- here("docs/approach_schema.svg")
use_svg_path <- function(k) {
  here(sprintf("docs/approach_schema_use%s.svg", k))
}
# Left out of the per-use drawings only; the full drawing keeps every row.
hidden_types <- c("prediction", "output")
show_relations <- TRUE
wrap_width <- 26

# One row per node type: the cluster that holds it (types sharing a cluster
# are drawn together), its fill colour, and the column it is pinned to (NA:
# placed by `dot` from its relations). Pinned columns keep the reading order
# questions -> design and data -> metrics -> output.
node_types <- data.frame(
  type = c(
    "scientific question",
    "prediction",
    "technical question",
    "design",
    "method",
    "parameter",
    "database",
    "input",
    "dataset",
    "control",
    "truth",
    "rule",
    "metric",
    "output"
  ),
  cluster = c(
    "uses",
    "predictions",
    "questions",
    "grid",
    "grid",
    "grid",
    "grid",
    "grid",
    "data",
    "data",
    "data",
    "rules",
    "metrics",
    "output"
  ),
  fill = c(
    "#F4A582",
    "#FDDBC7",
    "#FEE08B",
    "#C7E9C0",
    "#A1D99B",
    "#E5F5E0",
    "#9ECAE1",
    "#DEEBF7",
    "#BCBDDC",
    "#DADAEB",
    "#F2F0F7",
    "#E0E0E0",
    "#FFFFB3",
    "#FBB4AE"
  ),
  column = c(1, 1, 2, NA, NA, NA, NA, NA, NA, NA, NA, NA, 7, 8),
  title = c(
    "Scientific questions",
    "Predictions",
    "Technical questions",
    "Design",
    "Methods",
    "Parameters",
    "Reference databases",
    "Inputs",
    "Datasets",
    "Negative controls",
    "Truth",
    "Rules",
    "Metrics",
    "Output"
  )
)
# Cluster titles, and the colour of the relations that leave each cluster.
cluster_titles <- c(
  uses = "Scientific questions: three uses (OD 51)",
  predictions = "Predictions (OD 45)",
  questions = "Technical questions",
  grid = "Assignment grid",
  data = "Datasets, controls and truth",
  rules = "Scoring and reading rules",
  metrics = "Metrics",
  output = "Output"
)
use_cluster_title <- "Scientific question (OD 51)"
cluster_edge_colors <- c(
  uses = "#C0504D",
  predictions = "#D98C5F",
  questions = "#B8860B",
  grid = "#3A8F5B",
  data = "#6A51A3",
  rules = "#525252",
  metrics = "#8C8C00",
  output = "#C0504D"
)
# One row per status: line style, colour (NA: the colour of the cluster the
# relation leaves) and legend text.
statuses <- data.frame(
  status = c("built", "to compute", "parked"),
  line = c("solid", "dashed", "dotted"),
  color = c(NA, NA, "#9E9E9E"),
  legend = c(
    "built: read in a chapter",
    "to compute: open ROADMAP item",
    "parked: docs/parked.md"
  )
)

edges <- read.csv(
  csv_path,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character(0),
  encoding = "UTF-8"
)

# ---- Checks -----------------------------------------------------------------
required <- c(
  "from_type",
  "from",
  "relation",
  "to_type",
  "to",
  "status",
  "uses"
)
missing_cols <- setdiff(c(required, "source"), names(edges))
if (length(missing_cols) > 0) {
  stop("Missing columns in ", csv_path, ": ", toString(missing_cols))
}
empty <- which(rowSums(trimws(as.matrix(edges[required])) == "") > 0)
if (length(empty) > 0) {
  stop("Empty field in CSV rows ", toString(empty + 1))
}
unknown_types <- setdiff(c(edges$from_type, edges$to_type), node_types$type)
if (length(unknown_types) > 0) {
  stop(
    "Unknown node types (add them to `node_types`): ",
    toString(unknown_types)
  )
}
unknown_status <- setdiff(edges$status, statuses$status)
if (length(unknown_status) > 0) {
  stop(
    "Unknown status: ",
    toString(unknown_status),
    "; expected one of ",
    toString(statuses$status)
  )
}
duplicated_edges <- duplicated(edges[c("from", "to")])
if (any(duplicated_edges)) {
  stop(
    "Relation given twice: ",
    toString(paste(
      edges$from[duplicated_edges],
      "->",
      edges$to[duplicated_edges]
    ))
  )
}

all_nodes <- unique(rbind(
  data.frame(label = edges$from, type = edges$from_type),
  data.frame(label = edges$to, type = edges$to_type)
))
two_types <- unique(all_nodes$label[duplicated(all_nodes$label)])
if (length(two_types) > 0) {
  stop("Node given two types: ", toString(two_types))
}
# Labels that differ only by case or spacing are a typo, not two nodes.
squashed <- tolower(gsub("[[:space:]]+", " ", trimws(all_nodes$label)))
near_twins <- all_nodes$label[squashed %in% squashed[duplicated(squashed)]]
if (length(near_twins) > 0) {
  stop("Labels differing only by case or spacing: ", toString(near_twins))
}

# ---- Uses -------------------------------------------------------------------
use_number <- function(label) {
  sub("^Use ([0-9]+):.*$", "\\1", label)
}
use_labels <- all_nodes$label[all_nodes$type == "scientific question"]
not_numbered <- use_labels[use_number(use_labels) == use_labels]
if (length(not_numbered) > 0) {
  stop(
    "Scientific questions must be labelled \"Use <k>: ...\": ",
    toString(not_numbered)
  )
}
names(use_labels) <- use_number(use_labels)
use_labels <- use_labels[order(as.integer(names(use_labels)))]

row_uses <- lapply(strsplit(trimws(edges$uses), "[[:space:]]+"), function(u) {
  if (identical(u, "all")) {
    names(use_labels)
  } else {
    u
  }
})
unknown_uses <- setdiff(unlist(row_uses), names(use_labels))
if (length(unknown_uses) > 0) {
  stop(
    "Unknown values in `uses`: ",
    toString(unknown_uses),
    "; expected \"all\" or numbers among ",
    toString(names(use_labels))
  )
}
# A row of a scientific question is drawn for that question alone.
own_use <- ifelse(
  edges$from_type == "scientific question",
  use_number(edges$from),
  ifelse(
    edges$to_type == "scientific question",
    use_number(edges$to),
    NA_character_
  )
)
wrong_use <- which(
  !is.na(own_use) & !mapply(identical, row_uses, as.list(own_use))
)
if (length(wrong_use) > 0) {
  stop(
    "CSV rows ",
    toString(wrong_use + 1),
    ": a row of \"Use <k>\" must have `uses` = k"
  )
}

# ---- Nodes ------------------------------------------------------------------
# The nodes of a set of relations, in the order of `node_types`, with their
# id, cluster, fill, pinned column and status (the one their relations all
# share, else "built").
build_nodes <- function(edges) {
  nodes <- unique(rbind(
    data.frame(label = edges$from, type = edges$from_type),
    data.frame(label = edges$to, type = edges$to_type)
  ))
  nodes <- nodes[order(match(nodes$type, node_types$type)), ]
  nodes$id <- paste0("n", seq_len(nrow(nodes)))
  nodes$cluster <- node_types$cluster[match(nodes$type, node_types$type)]
  nodes$fill <- node_types$fill[match(nodes$type, node_types$type)]
  nodes$column <- node_types$column[match(nodes$type, node_types$type)]
  nodes$status <- vapply(
    nodes$label,
    function(l) {
      s <- unique(edges$status[edges$from == l | edges$to == l])
      if (length(s) == 1) {
        s
      } else {
        "built"
      }
    },
    character(1),
    USE.NAMES = FALSE
  )
  nodes
}

# A DOT string: quotes escaped, long labels wrapped on several lines.
dot_label <- function(x, width = wrap_width) {
  vapply(
    x,
    function(s) {
      s <- gsub("\"", "\\\\\"", s)
      paste(strwrap(s, width = width), collapse = "\\n")
    },
    character(1),
    USE.NAMES = FALSE
  )
}

node_line <- function(nodes, i) {
  line <- statuses$line[statuses$status == nodes$status[i]]
  style <- if (line == "solid") {
    "filled,rounded"
  } else {
    paste0("filled,rounded,", line)
  }
  sprintf(
    "    %s [label=\"%s\", fillcolor=\"%s\", style=\"%s\"];",
    nodes$id[i],
    dot_label(nodes$label[i]),
    nodes$fill[i],
    style
  )
}

# A cluster holding several types boxes each type of two nodes or more on
# its own, so that methods, databases or datasets are not interleaved.
cluster_block <- function(nodes, cl, titles) {
  rows <- which(nodes$cluster == cl)
  types <- unique(nodes$type[rows])
  body <- unlist(lapply(types, function(ty) {
    type_rows <- rows[nodes$type[rows] == ty]
    lines <- vapply(type_rows, \(i) node_line(nodes, i), character(1))
    if (length(types) == 1 || length(type_rows) < 2) {
      return(lines)
    }
    c(
      sprintf("    subgraph cluster_%s_%s {", cl, gsub("[^a-z]", "_", ty)),
      sprintf(
        "      label=\"%s\"; fontsize=12; style=\"rounded,dashed\";",
        node_types$title[node_types$type == ty]
      ),
      paste0("  ", lines),
      "    }"
    )
  }))
  c(
    sprintf("  subgraph cluster_%s {", cl),
    sprintf("    label=\"%s\";", titles[[cl]]),
    body,
    "  }"
  )
}

# ---- Edges ------------------------------------------------------------------
edge_line <- function(edges, nodes, i) {
  from <- match(edges$from[i], nodes$label)
  to <- match(edges$to[i], nodes$label)
  status <- statuses[statuses$status == edges$status[i], ]
  color <- if (is.na(status$color)) {
    cluster_edge_colors[[nodes$cluster[from]]]
  } else {
    status$color
  }
  label <- if (show_relations) {
    sprintf(
      ", label=\"%s\", fontcolor=\"%s\"",
      dot_label(edges$relation[i], width = 18),
      color
    )
  } else {
    ""
  }
  sprintf(
    "  %s -> %s [style=%s, color=\"%s\"%s];",
    nodes$id[from],
    nodes$id[to],
    status$line,
    color,
    label
  )
}

# ---- Pinned columns ---------------------------------------------------------
# The nodes of one pinned column share a rank; an invisible chain between the
# columns, as long as the gap between their numbers, keeps them in order and
# leaves room for the unpinned clusters in between.
column_lines <- function(nodes) {
  columns <- sort(unique(stats::na.omit(nodes$column)))
  blocks <- vapply(
    columns,
    function(k) {
      sprintf(
        "  { rank=same; %s; }",
        paste(nodes$id[which(nodes$column == k)], collapse = "; ")
      )
    },
    character(1)
  )
  anchors <- nodes$id[match(columns, nodes$column)]
  chain <- if (length(columns) > 1) {
    sprintf(
      "  %s -> %s [style=invis, minlen=%d];",
      anchors[-length(anchors)],
      anchors[-1],
      as.integer(diff(columns))
    )
  } else {
    character(0)
  }
  c(blocks, chain)
}

# ---- Legend -----------------------------------------------------------------
# The node types and statuses drawn, one swatch or one sample edge each.
legend_block <- function(nodes, edges) {
  types <- node_types[node_types$type %in% nodes$type, ]
  shown <- statuses[statuses$status %in% edges$status, ]
  legend_rows <- sprintf(
    "<TR><TD BGCOLOR=\"%s\" WIDTH=\"18\"> </TD><TD ALIGN=\"LEFT\">%s</TD></TR>",
    types$fill,
    types$type
  )
  c(
    "  subgraph cluster_legend {",
    "    label=\"Legend\";",
    "    style=\"rounded\"; color=\"grey60\";",
    paste0(
      "    legend_types [shape=plaintext, style=\"\", label=<",
      "<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLSPACING=\"2\">",
      paste(legend_rows, collapse = ""),
      "</TABLE>>];"
    ),
    unlist(lapply(seq_len(nrow(shown)), function(i) {
      color <- if (is.na(shown$color[i])) "grey35" else shown$color[i]
      c(
        sprintf("    legend_from%d [shape=point, width=0.05, label=\"\"];", i),
        sprintf(
          "    legend_to%d [shape=plaintext, style=\"\", label=\"%s\"];",
          i,
          shown$legend[i]
        ),
        sprintf(
          "    legend_from%d -> legend_to%d [style=%s, color=\"%s\"];",
          i,
          i,
          shown$line[i],
          color
        )
      )
    })),
    "  }"
  )
}

# ---- DOT and SVG ------------------------------------------------------------
# The DOT source of a set of relations. `title`, when given, is written above
# the graph; it is set last so that the clusters do not inherit it.
schema_dot <- function(edges, titles = cluster_titles, title = NULL) {
  nodes <- build_nodes(edges)
  c(
    "digraph approach {",
    "  graph [rankdir=LR, newrank=true, splines=true, nodesep=0.12,",
    "         ranksep=0.45, fontname=\"Helvetica\", fontsize=16,",
    "         labeljust=l, style=\"rounded,filled\", fillcolor=\"#FAFAFA\",",
    "         color=\"grey75\"];",
    "  node [shape=box, fontname=\"Helvetica\", fontsize=11, margin=\"0.08,0.04\",",
    "        color=\"grey35\"];",
    "  edge [fontname=\"Helvetica\", fontsize=8, fontcolor=\"grey35\",",
    "        color=\"grey45\", arrowsize=0.6];",
    unlist(lapply(
      unique(nodes$cluster),
      cluster_block,
      nodes = nodes,
      titles = titles
    )),
    legend_block(nodes, edges),
    column_lines(nodes),
    vapply(
      seq_len(nrow(edges)),
      \(i) edge_line(edges, nodes, i),
      character(1)
    ),
    if (!is.null(title)) {
      sprintf(
        "  label=\"%s\"; labelloc=t; fontsize=22;",
        dot_label(title, width = 80)
      )
    },
    "}"
  )
}

render_svg <- function(dot, path) {
  dot_path <- tempfile(fileext = ".dot")
  writeLines(dot, dot_path, useBytes = TRUE)
  out <- suppressWarnings(system2(
    "dot",
    c("-Tsvg", shQuote(dot_path), "-o", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("dot failed on ", path, ":\n", paste(out, collapse = "\n"))
  }
  invisible(path)
}

if (!nzchar(Sys.which("dot"))) {
  stop("Graphviz `dot` is not on PATH.")
}

render_svg(schema_dot(edges), svg_path)
message(
  "approach_schema: ",
  nrow(all_nodes),
  " nodes, ",
  nrow(edges),
  " relations -> ",
  svg_path
)

drawn <- !(edges$from_type %in% hidden_types | edges$to_type %in% hidden_types)
use_titles <- cluster_titles
use_titles[["uses"]] <- use_cluster_title
for (k in names(use_labels)) {
  rows <- drawn & vapply(row_uses, \(u) k %in% u, logical(1))
  use_edges <- edges[rows, ]
  path <- use_svg_path(k)
  render_svg(
    schema_dot(use_edges, use_titles, title = use_labels[[k]]),
    path
  )
  message(
    "approach_schema: use ",
    k,
    ", ",
    nrow(build_nodes(use_edges)),
    " nodes, ",
    nrow(use_edges),
    " relations -> ",
    path
  )
}
