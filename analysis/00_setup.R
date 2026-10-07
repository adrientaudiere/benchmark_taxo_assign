# Shared setup for every chapter of analysis/. First line of each qmd:
#   source(here::here("analysis/00_setup.R"))
#
# Loads the packages and the pqverse checkouts, defines the project-wide
# constants (rank order, default-settings filter, preference pattern) and the
# helpers used across chapters (save_fig, cache_write / cache_read,
# plot_tc_metrics_mock).

library("conflicted")
library("targets")
library("here")
library("tibble")
library("tidyr")
library("dplyr")
library("ggplot2")
library("patchwork")

conflicted::conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::rename,
  dplyr::count,
  .quiet = TRUE
)

here::i_am("analysis/00_setup.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/values_map.R"))
load_pqverse(c("MiscMetabar", "comparpq", "dbpq", "greenAlgoR", "tidypq"))

# ---- stores and directories ------------------------------------------------

# BENCHMARK_ANALYSIS_MINI=TRUE renders the chapters on the smoke-test stores
# of the *_mini targets projects, with their own cache and figure folders, so
# mini results can never overwrite production figures:
#   BENCHMARK_ANALYSIS_MINI=TRUE quarto render analysis/01_load_and_clean.qmd
analysis_mini <- as.logical(Sys.getenv(
  "BENCHMARK_ANALYSIS_MINI",
  unset = "FALSE"
))
if (is.na(analysis_mini)) {
  stop("BENCHMARK_ANALYSIS_MINI must be TRUE or FALSE.")
}
variant_suffix <- if (analysis_mini) "_mini" else ""

store_assign_taxo <- here(paste0("store_assign_taxo", variant_suffix))
store_cross_val <- here(paste0("store_cross_val", variant_suffix))

# Datasets chapter 01 scores: the Pauvert mock of the production projects and
# every biological dataset of `config.R::bio_datasets` that has a truth table
# in `config.R::mock_sanger_tables`. Only the ones whose assignment store
# exists are kept, so the chapter renders before the other runs exist; `store`
# holds the assignments and `dada2_store` the objects the truth is built from
# (`d_asv` and `d_vs`, never the assignment inputs: the ITSx input carries
# trimmed sequences).
analysis_datasets <- tibble::tibble(
  dataset = mock_sanger_tables$dataset,
  store = ifelse(
    mock_sanger_tables$dataset == "pauvert",
    here(paste0("store_assign_taxo", variant_suffix)),
    here(paste0(
      "store_assign_taxo_",
      mock_sanger_tables$dataset,
      variant_suffix
    ))
  ),
  dada2_store = ifelse(
    mock_sanger_tables$dataset == "pauvert",
    here("store_dada2"),
    here(paste0("store_dada2_", mock_sanger_tables$dataset))
  )
)
# ITSx regions and units of each dataset (OD 34), as config.R derives them for
# the running project: ITS1 with ASV and OTU for the Pauvert mock, the
# `bio_datasets` row otherwise.
dataset_regions <- function(dataset) {
  if (dataset == "pauvert") {
    return("ITS1")
  }
  region <- bio_datasets$itsx_region[bio_datasets$dataset == dataset]
  strsplit(region, ",", fixed = TRUE)[[1]]
}
dataset_units <- function(dataset) {
  if (
    dataset == "pauvert" ||
      bio_datasets$otu_input[bio_datasets$dataset == dataset]
  ) {
    c("asv", "otu")
  } else {
    "asv"
  }
}
analysis_datasets$regions <- lapply(analysis_datasets$dataset, dataset_regions)
analysis_datasets$units <- lapply(analysis_datasets$dataset, dataset_units)
# A store left half-built by an interrupted run (Tedersoo, 2026-10-01) exists
# without its assignment tables, so a dataset is kept only when every
# `d_all_taxo_<region>_<unit>` target is in its store.
assignments_complete <- function(store, regions, units) {
  dir.exists(store) &&
    all(targets::tar_exist_objects(
      paste0("d_all_taxo_", as.vector(outer(regions, units, paste, sep = "_"))),
      store = store
    ))
}
complete <- dir.exists(analysis_datasets$dada2_store) &
  mapply(
    assignments_complete,
    analysis_datasets$store,
    analysis_datasets$regions,
    analysis_datasets$units
  )
if (any(dir.exists(analysis_datasets$store) & !complete)) {
  message(
    "analysis: assignment store incomplete, dataset skipped: ",
    paste(
      analysis_datasets$dataset[
        dir.exists(analysis_datasets$store) & !complete
      ],
      collapse = ", "
    )
  )
}
analysis_datasets <- analysis_datasets[complete, ]
figures_dir <- if (analysis_mini) here("figures/mini") else here("figures")
cache_dir <- if (analysis_mini) {
  here("analysis/_cache/mini")
} else {
  here("analysis/_cache")
}
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
message(
  "analysis: ",
  if (analysis_mini) "MINI" else "production",
  " stores (",
  basename(store_assign_taxo),
  ", ",
  basename(store_cross_val),
  ")"
)

# Same grid as pipelines/assign_taxo.R, for one dataset: only the column-name
# machinery (full_name, method, db, ...) is used here, never db_path.
values_map_for <- function(dataset) {
  build_values_map(
    dbs = db_list,
    regions = dataset_regions(dataset),
    units = dataset_units(dataset)
  )
}
# Every name any dataset can produce, to read the method, database, parameters
# and region of a scored column by joining on its full_name (C22-S3).
name_map <- build_values_map(
  dbs = db_list,
  regions = c("ITS1", "ITS2", "full"),
  units = c("asv", "otu")
)

# ---- constants -------------------------------------------------------------

tax_order <- c(
  "Kingdom",
  "Phylum",
  "Class",
  "Order",
  "Family",
  "Genus",
  "Gen_sp"
)
single_methods <- c("dada2", "sintax", "lca", "blastn")
# Negative controls of every input: the shuffled sequences of
# comparpq::add_shuffle_seq_pq() and the non-fungal ones of
# add_external_seq_pq(). They carry no reads, so an abundance filter must skip
# them (chapter 01, M1.5).
control_pattern <- "^(fake_|external_)"
# One point shape per database: ggplot's default shape palette stops at 6
# values and silently drops the points of the other databases.
db_shapes <- setNames(
  c(16, 17, 15, 18, 1, 2, 0, 5, 6)[seq_along(db_list)],
  db_list
)

# (method, db, bootstrap) column preferred by the `preference` consensus
# strategy (Q3.5); preference_method and preference_db come from config.R.
# Unprefixed: the consensus of an input adds the input, "<input>_" (chapter
# 01, one consensus per region).
preference_pattern <- full_name_for(preference_method, preference_db, 0.5)
stopifnot(paste0("ITS1_asv_", preference_pattern) %in% name_map$full_name)

# ---- helpers ---------------------------------------------------------------

# Rows of res_comp_tax at the "default settings" used by every publication
# figure: the four single methods, bootstrap 0.5 (dada2 / sintax),
# lca_cutoff = 0.51 (lca, the QIIME 2 / Bokulich consensus, OD 33),
# rel_majority vote with min_id = 95 (blastn, the assign_blastn() default used
# before the min_id sweep). Same four rules as the `consensus_names` of chapter
# 01, which votes over exactly these columns (decision 0.4d).
# The itsx_ and otu_ rows are dropped with the consensus ones: strip the prefix
# of `algo` first to keep them (chapter 01, Q4 and Q5 views).
filter_default_settings <- function(res) {
  res |>
    dplyr::filter(algo %in% single_methods) |>
    dplyr::filter(
      (algo %in% c("dada2", "sintax") & bootstrap_num == 0.5) |
        (algo == "lca" & lca_cutoff == 0.51) |
        (algo == "blastn" & vote_algorithm == "rel_majority" & min_id == 95)
    )
}

# Adds the columns every figure groups by, read from `name_map` joined on the
# scored column (`method_db` = its full_name, C22-S3) rather than parsed from
# the name: algo (the method), db, bootstrap_num (dada2 / sintax threshold),
# lca_cutoff, vote_algorithm and min_id (blastn), and the ITSx `region`.
# Consensus columns ("<input>_<strategy>_consensus", chapter 01) are not in
# `name_map`: their region is read from the prefix, their algo is the
# strategy, and they have no database. Adds the TAR / TDR aliases of the
# manuscript (TAR = PPV = TP/(TP+FP); TDR = TPR = TP/(TP+FN)).
enrich_metrics <- function(res) {
  cols <- name_map |>
    dplyr::transmute(
      method_db = full_name,
      algo = method,
      db,
      region,
      bootstrap_num = ifelse(method == "lca", NA_real_, min_bootstrap),
      lca_cutoff,
      vote_algorithm,
      min_id
    )
  res <- res |>
    dplyr::select(-dplyr::any_of(setdiff(names(cols), "method_db"))) |>
    dplyr::left_join(cols, by = "method_db")
  consensus_input <- "^(ITS1|ITS2|full)_(asv|otu)_"
  res <- res |>
    dplyr::mutate(
      is_consensus = is.na(algo),
      region = ifelse(
        is_consensus,
        sub(paste0(consensus_input, ".*$"), "\\1", method_db),
        region
      ),
      algo = ifelse(is_consensus, sub(consensus_input, "", method_db), algo),
      consensus_strict = is_consensus & grepl("_strict_", method_db),
      nb_agree = dplyr::case_when(
        grepl("_nb2_", method_db) ~ 2L,
        grepl("_nb3_", method_db) ~ 3L,
        TRUE ~ 1L
      )
    )
  unknown <- res$method_db[
    res$is_consensus & !grepl(consensus_input, res$method_db)
  ]
  if (length(unknown) > 0) {
    stop("enrich_metrics(): unknown column name ", unique(unknown)[1])
  }

  dplyr::bind_rows(
    res,
    res |> dplyr::filter(metrics == "PPV") |> dplyr::mutate(metrics = "TAR"),
    res |> dplyr::filter(metrics == "TPR") |> dplyr::mutate(metrics = "TDR")
  )
}

# The view of res_comp_tax chapters 02 to 05 read: one dataset, one ITSx
# region, one unit, and one reading of the consensus columns. Chapter 01 scores
# every region of every dataset that has a store (OD 34), three units per
# region (ASV, OTU, OTU_derived), and every consensus column twice (decision
# 2026-09-22: `external_scoring` "matrix" and "aside", the consensus having no
# database of its own). Without this filter a chapter would average regions
# and units and count each consensus twice, silently. Q4 (regions), Q5 (units)
# and the dataset comparison read the cache directly.
default_view <- function(
  res,
  keep_unit = "ASV",
  keep_dataset = NULL,
  keep_region = NULL,
  consensus_scoring = "aside"
) {
  if (is.null(keep_dataset)) {
    keep_dataset <- res$dataset[[1]]
  }
  if (is.null(keep_region)) {
    # The first region of the dataset, in the order ITS1, ITS2, full.
    present <- unique(res$region[res$dataset == keep_dataset])
    keep_region <- intersect(c("ITS1", "ITS2", "full"), present)[1]
  }
  res |>
    dplyr::filter(
      unit == keep_unit,
      dataset == keep_dataset,
      region == keep_region
    ) |>
    dplyr::filter(!is_consensus | external_scoring == consensus_scoring)
}

# The views chapters 02 to 05 show, one tab each (developer, 2026-10-02):
# every dataset x ITSx region of `res`, ASV unit, Pauvert first and the
# regions in the order ITS1, ITS2, full. `id` suffixes the saved figures.
analysis_views <- function(res) {
  res |>
    dplyr::distinct(dataset, region) |>
    dplyr::mutate(
      region = factor(region, levels = c("ITS1", "ITS2", "full")),
      not_pauvert = dataset != "pauvert"
    ) |>
    dplyr::arrange(not_pauvert, dataset, region) |>
    dplyr::mutate(
      region = as.character(region),
      id = paste0(dataset, "_", region),
      label = paste0(dataset, " · ", region)
    ) |>
    dplyr::select(dataset, region, id, label)
}

# Write one panel-tabset with a tab per view, in a chunk with
# `results = "asis"`. `fun(res_view, view)` receives the default_view() of
# that dataset x region and the view row, and prints its figures and tables.
# `fun` must print() a figure BEFORE save_fig(): the ggsave() devices opened
# in between delay the capture of the plot, which then lands in the next tab.
view_tabs <- function(res, views, fun) {
  cat("\n::: {.panel-tabset}\n\n")
  for (i in seq_len(nrow(views))) {
    cat("### ", views$label[i], "\n\n", sep = "")
    fun(
      default_view(
        res,
        keep_dataset = views$dataset[i],
        keep_region = views$region[i]
      ),
      views[i, ]
    )
    cat("\n\n")
  }
  cat(":::\n\n")
}

# The default_view() of every view, stacked, with a `view` column (the tab
# label, as an ordered factor) and its `view_id`: what the chapter summaries
# read, so that they filter exactly as the tabs below them.
view_rows <- function(res, views, keep_unit = "ASV") {
  rows <- lapply(seq_len(nrow(views)), function(i) {
    default_view(
      res,
      keep_unit = keep_unit,
      keep_dataset = views$dataset[i],
      keep_region = views$region[i]
    ) |>
      dplyr::mutate(view = views$label[i], view_id = views$id[i])
  })
  dplyr::bind_rows(rows) |>
    dplyr::mutate(view = factor(view, levels = views$label))
}

# "0.880" or "0.873–0.903": the range of a few values, for summary sentences.
fmt_range <- function(x, digits = 3) {
  x <- round(x, digits)
  if (min(x) == max(x)) {
    return(format(min(x), nsmall = digits))
  }
  # An en dash between two negatives reads as a minus: "to" then.
  sep <- if (min(x) < 0) " to " else "–"
  paste0(format(min(x), nsmall = digits), sep, format(max(x), nsmall = digits))
}

# One phrase naming the best `what` (a column of `tbl`, e.g. "algo" or "db")
# of each view by `value`, for the summaries: "**dada2** in every view
# (0.823–0.928)" when one wins everywhere, else one "**x** on view (value)"
# per view. Computed, so a summary never states a winner the data dropped.
# `show_value = FALSE` names the winners only (when `value` is a ranking key
# rather than a figure to quote).
best_phrase <- function(tbl, what, show_value = TRUE) {
  best <- tbl |>
    dplyr::group_by(view) |>
    dplyr::slice_max(value, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(view)
  if (dplyr::n_distinct(best[[what]]) == 1) {
    value <- if (show_value) paste0(" (", fmt_range(best$value), ")") else ""
    return(paste0("**", best[[what]][1], "** in every view", value))
  }
  value <- if (show_value) {
    paste0(" (", format(round(best$value, 3), nsmall = 3), ")")
  } else {
    ""
  }
  paste(
    paste0("**", best[[what]], "** on ", best$view, value),
    collapse = "; "
  )
}

# A kable printed from inside view_tabs().
cat_kable <- function(x, ...) {
  cat(knitr::kable(x, format = "pipe", ...), sep = "\n")
  cat("\n\n")
}

# The benchmark_costs of one dataset's assignment store, stopped when none of
# its assignment names matches the scored columns (a store built before the
# 2026-09-17 grid would silently join to nothing).
costs_for <- function(dataset, res) {
  store <- analysis_datasets$store[analysis_datasets$dataset == dataset]
  costs <- targets::tar_read(benchmark_costs, store = store)
  if (sum(costs$phase %in% res$method_db) == 0) {
    stop(
      "No assignment of benchmark_costs matches res_comp_tax: ",
      basename(store),
      " predates the current grid. Rerun its assign_taxo project."
    )
  }
  costs
}

# Which real units of an ITSx region input ITSx trimmed (OD 38): TRUE when
# the sequence of `d_<region>_asv_real` differs from the denoised ASV of the
# same name in `d_asv`. itsx_pq(keep_undetected = TRUE) keeps a unit whose
# region it did not find on its whole sequence, and nothing in the tax_table
# says so; the units merged after trimming keep the name of one of them.
# Named logical vector over the real units of the input.
itsx_trimmed <- function(store, dada2_store, region) {
  d_asv <- targets::tar_read(d_asv, store = dada2_store)
  orig <- stats::setNames(as.character(d_asv@refseq), taxa_names(d_asv))
  real <- targets::tar_read_raw(paste0("d_", region, "_asv_real"), store = store)
  now <- stats::setNames(as.character(real@refseq), taxa_names(real))
  missing <- setdiff(names(now), names(orig))
  if (length(missing) > 0) {
    stop(
      length(missing), " units of d_", region, "_asv_real are not in d_asv of ",
      basename(dada2_store)
    )
  }
  now != orig[names(now)]
}

# Write a figure as pdf + png (300 dpi) under figures/ and return it invisibly.
save_fig <- function(plot, name, width, height, dpi = 300) {
  # cairo_pdf: the default pdf device replaces non-ASCII labels such as → and −.
  ggsave(
    file.path(figures_dir, paste0(name, ".pdf")),
    plot,
    width = width,
    height = height,
    device = cairo_pdf
  )
  ggsave(
    file.path(figures_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = dpi
  )
  invisible(plot)
}

# Objects handed from one chapter to the next live in analysis/_cache/.
cache_write <- function(object, name) {
  saveRDS(object, file.path(cache_dir, paste0(name, ".rds")))
  invisible(object)
}
cache_read <- function(name) {
  path <- file.path(cache_dir, paste0(name, ".rds"))
  if (!file.exists(path)) {
    stop("Missing ", path, ": render analysis/01_load_and_clean.qmd first.")
  }
  readRDS(path)
}

# Generic faceted bar plot of tc_metrics_mock() output (exploratory views).
plot_tc_metrics_mock <- function(
  tc_metrics_mock,
  metric = NULL,
  modality_x = NULL,
  modality_line = NULL,
  modality_column = NULL,
  modality_fill = NULL,
  taxonomic_ranks = NULL,
  legend = TRUE
) {
  if (is.null(taxonomic_ranks)) {
    taxonomic_ranks <- unique(tc_metrics_mock$tax_level)
  } else {
    tc_metrics_mock <- tc_metrics_mock |>
      dplyr::filter(tax_level %in% taxonomic_ranks)
  }
  if (is.null(modality_fill)) {
    modality_fill <- "method_db"
  }
  if (is.null(modality_x)) {
    modality_x <- "tax_level"
  }
  tc_metrics_mock$tax_level <- factor(
    tc_metrics_mock$tax_level,
    taxonomic_ranks
  )

  if (!is.null(metric)) {
    tc_metrics_mock <- dplyr::filter(tc_metrics_mock, metrics == metric)
  } else {
    if (
      modality_fill != "metrics" &&
        !identical(modality_column, "metrics") &&
        !identical(modality_line, "metrics") &&
        modality_x != "metrics"
    ) {
      stop(
        "Either choose one metric (parameter `metric`) or map 'metrics' to ",
        "one of modality_x, modality_line, modality_column or modality_fill."
      )
    }
  }

  p <- tc_metrics_mock |>
    ggplot(aes(
      fill = .data[[modality_fill]],
      y = values,
      x = .data[[modality_x]]
    )) +
    geom_bar(stat = "identity", position = position_dodge())

  if (is.null(modality_column) && !is.null(modality_line)) {
    p <- p + facet_grid(.data[[modality_line]] ~ ., scales = "free_y")
  } else if (!is.null(modality_column) && is.null(modality_line)) {
    p <- p + facet_grid(. ~ .data[[modality_column]], scales = "free_x")
  } else if (!is.null(modality_column) && !is.null(modality_line)) {
    p <- p +
      facet_grid(
        vars(.data[[modality_line]]),
        vars(.data[[modality_column]]),
        scales = "free"
      )
  }
  if (!legend) {
    p <- p + no_legend()
  }
  p
}
