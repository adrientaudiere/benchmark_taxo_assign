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

conflicted::conflicts_prefer(dplyr::filter, dplyr::select, dplyr::rename,
                             dplyr::count, .quiet = TRUE)

here::i_am("analysis/00_setup.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/values_map.R"))
load_pqverse(c("MiscMetabar", "comparpq", "dbpq", "greenAlgoR", "tidypq"))

# ---- stores and directories ------------------------------------------------

store_assign_taxo <- here("store_assign_taxo")
store_cross_val   <- here("store_cross_val")
figures_dir       <- here("figures")
cache_dir         <- here("analysis/_cache")
dir.create(figures_dir, showWarnings = FALSE)
dir.create(cache_dir,   showWarnings = FALSE)

# Same grid as pipelines/assign_taxo.R: only the column-name machinery
# (full_name, method, db, ...) is used here, never db_path.
values_map <- build_values_map(dbs = db_list)

# ---- constants -------------------------------------------------------------

tax_order      <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Gen_sp")
single_methods <- c("dada2", "sintax", "lca", "blastn")

# (method, db, bootstrap) column preferred by the `preference` consensus
# strategy (Q3.5). Must be one of values_map$full_name.
preference_pattern <- "sintax__EUK_ITS_v2___0.5"
stopifnot(preference_pattern %in% values_map$full_name)
preference_method  <- "sintax"
preference_db      <- "EUK_ITS_v2"

# ---- helpers ---------------------------------------------------------------

# Rows of res_comp_tax at the "default settings" used by every publication
# figure: the four single methods, bootstrap 0.5 (dada2 / sintax / lca) or
# rel_majority vote (blastn), production databases only.
filter_default_settings <- function(res) {
  res |>
    dplyr::filter(algo %in% single_methods) |>
    dplyr::filter(!startsWith(db, "mini_")) |>
    dplyr::filter(
      (algo %in% c("dada2", "sintax", "lca") & bootstrap_num == 0.5) |
        (algo == "blastn" & vote_algorithm == "rel_majority")
    )
}

# Write a figure as pdf + png (300 dpi) under figures/ and return it invisibly.
save_fig <- function(plot, name, width, height, dpi = 300) {
  ggsave(file.path(figures_dir, paste0(name, ".pdf")), plot,
         width = width, height = height)
  ggsave(file.path(figures_dir, paste0(name, ".png")), plot,
         width = width, height = height, dpi = dpi)
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
plot_tc_metrics_mock <- function(tc_metrics_mock,
                                 metric = NULL,
                                 modality_x = NULL,
                                 modality_line = NULL,
                                 modality_column = NULL,
                                 modality_fill = NULL,
                                 taxonomic_ranks = NULL,
                                 legend = TRUE) {
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
  tc_metrics_mock$tax_level <- factor(tc_metrics_mock$tax_level, taxonomic_ranks)

  if (!is.null(metric)) {
    tc_metrics_mock <- dplyr::filter(tc_metrics_mock, metrics == metric)
  } else {
    if (modality_fill != "metrics" &&
        !identical(modality_column, "metrics") &&
        !identical(modality_line, "metrics") &&
        modality_x != "metrics") {
      stop("Either choose one metric (parameter `metric`) or map 'metrics' to ",
           "one of modality_x, modality_line, modality_column or modality_fill.")
    }
  }

  p <- tc_metrics_mock |>
    ggplot(aes(fill = .data[[modality_fill]], y = values, x = .data[[modality_x]])) +
    geom_bar(stat = "identity", position = position_dodge())

  if (is.null(modality_column) && !is.null(modality_line)) {
    p <- p + facet_grid(.data[[modality_line]] ~ ., scales = "free_y")
  } else if (!is.null(modality_column) && is.null(modality_line)) {
    p <- p + facet_grid(. ~ .data[[modality_column]], scales = "free_x")
  } else if (!is.null(modality_column) && !is.null(modality_line)) {
    p <- p + facet_grid(vars(.data[[modality_line]]), vars(.data[[modality_column]]),
                        scales = "free")
  }
  if (!legend) {
    p <- p + no_legend()
  }
  p
}
