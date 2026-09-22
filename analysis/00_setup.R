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
analysis_datasets <- analysis_datasets[
  dir.exists(analysis_datasets$store) &
    dir.exists(analysis_datasets$dada2_store),
]
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

# Same grid as pipelines/assign_taxo.R: only the column-name machinery
# (full_name, method, db, ...) is used here, never db_path.
values_map <- build_values_map(dbs = db_list)

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
preference_pattern <- full_name_for(preference_method, preference_db, 0.5)
stopifnot(preference_pattern %in% values_map$full_name)

# ---- helpers ---------------------------------------------------------------

# Rows of res_comp_tax at the "default settings" used by every publication
# figure: the four single methods, bootstrap 0.5 (dada2 / sintax),
# lca_cutoff = 1 (lca, objectives_design §2),
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
        (algo == "lca" & lca_cutoff == 1) |
        (algo == "blastn" & vote_algorithm == "rel_majority" & min_id == 95)
    )
}

# Parses method_db into the columns every figure groups by: algo, db, the
# parameter of the row (bootstrap for dada2 / sintax, lca_cutoff for lca,
# vote_algorithm and min_id for blastn, ROADMAP 0.4) and the consensus
# metadata. Adds the TAR / TDR aliases of the manuscript (TAR = PPV =
# TP/(TP+FP); TDR = TPR = TP/(TP+FN)).
enrich_metrics <- function(res) {
  res <- res |>
    dplyr::mutate(bootstrap = gsub(".*__", "", method_db)) |>
    dplyr::mutate(algo = gsub("__.*", "", method_db)) |>
    dplyr::mutate(db = gsub(".*__", "", gsub("___.*", "", method_db))) |>
    dplyr::mutate(
      vote_algorithm = ifelse(
        grepl("blastn$", algo),
        stringr::str_extract(method_db, "(?<=\\.\\.\\.)[a-z_]+(?=\\.\\.\\.)"),
        NA_character_
      ),
      bootstrap_num = suppressWarnings(
        as.numeric(stringr::str_extract(bootstrap, "^[0-9.]+"))
      ),
      # lca rows carry their lca_cutoff where the others carry their
      # bootstrap ("lca__DB___0.8", ROADMAP 0.4): moved to its own column.
      lca_cutoff = ifelse(grepl("lca$", algo), bootstrap_num, NA_real_),
      bootstrap_num = ifelse(grepl("lca$", algo), NA_real_, bootstrap_num),
      # blastn rows end with "...<nb_voting>...<min_id>" (ROADMAP 0.4).
      min_id = ifelse(
        grepl("blastn$", algo),
        suppressWarnings(as.numeric(
          stringr::str_extract(
            method_db,
            "(?<=\\.\\.\\.[0-9]{1,6}\\.\\.\\.)[0-9]+$"
          )
        )),
        NA_real_
      )
    )

  res <- dplyr::bind_rows(
    res,
    res |> dplyr::filter(metrics == "PPV") |> dplyr::mutate(metrics = "TAR"),
    res |> dplyr::filter(metrics == "TPR") |> dplyr::mutate(metrics = "TDR")
  )

  res |>
    dplyr::mutate(
      # itsx_<method> rows (Q4) and otu_<method> rows (Q5) are single methods
      # on another input, not consensus columns.
      is_consensus = !sub("^(itsx|otu)_", "", algo) %in% single_methods,
      consensus_strict = is_consensus & grepl("_strict_", method_db),
      nb_agree = dplyr::case_when(
        grepl("_nb2_", method_db) ~ 2L,
        grepl("_nb3_", method_db) ~ 3L,
        TRUE ~ 1L
      ),
      db = ifelse(is_consensus, NA_character_, db)
    )
}

# The view of res_comp_tax chapters 02 to 05 read: one dataset, one input, and
# one reading of the consensus columns. Chapter 01 scores three inputs (ASV,
# OTU, OTU_derived), every dataset that has a store, and every consensus column
# twice (decision 2026-09-22: `external_scoring` "matrix" and "aside", the
# consensus having no database of its own). Without this filter a chapter would
# average the three inputs and count each consensus twice, silently.
# Q4 (ITSx), Q5 (OTU) and the dataset comparison read the cache directly.
default_view <- function(
  res,
  keep_unit = "ASV",
  keep_dataset = NULL,
  consensus_scoring = "aside"
) {
  if (is.null(keep_dataset)) {
    keep_dataset <- res$dataset[[1]]
  }
  res |>
    dplyr::filter(unit == keep_unit, dataset == keep_dataset) |>
    dplyr::filter(!is_consensus | external_scoring == consensus_scoring)
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
