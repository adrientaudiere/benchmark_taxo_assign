# Preparation of the assignment columns before scoring: cleaning of the labels,
# the Genus + Species columns, and the `ranks_df` that
# comparpq::tc_metrics_unit() takes. Moved out of analysis/01_load_and_clean.qmd
# on 2026-09-23 (first half of C22-S4) so that the in silico chapter scores its
# positive control with exactly the same preparation, without entering the
# res_comp_tax of chapter 01 (OD 29).
#
# Requires comparpq (taxtab_replace_pattern_by_NA), MiscMetabar (simplify_taxo)
# and tidypq (tax_table_to_df, mutate_taxa_pq).

# Same cleaning for every object: "<Genus>_sp" species, the incertae sedis and
# unclassified labels and the empty strings become NA. A truth built from
# reference headers must not hold these labels either, or the scoring can
# never match them (ROADMAP SA4, R/insilico_input.R).
clean_taxo <- function(physeq) {
  species_cols <- grep(
    "^Species_",
    colnames(physeq@tax_table),
    value = TRUE
  )
  physeq |>
    taxtab_replace_pattern_by_NA(
      patterns = ".*_sp$",
      progress_bar = FALSE,
      taxonomic_ranks = species_cols
    ) |>
    taxtab_replace_pattern_by_NA(
      patterns = c(
        ".*_incertae_sedis",
        ".*_Incertae_sedis",
        "unclassified.*",
        "^$"
      ),
      progress_bar = FALSE
    ) |>
    simplify_taxo()
}

# Assignment labels of a tax_table read with tidypq::tax_table_to_df(): the
# column names without their rank prefix, reference ranks, the `Taxa_name_*`
# columns and the ITSx columns of the unit (`ITSx_origin` and, since
# 2026-10-07, `ITSx_detected`, written by MiscMetabar::itsx_pq(), OD 34)
# excluded.
assignment_labels <- function(taxtab_df) {
  labels <- gsub(
    "Kingdom|Phylum|Phyla|Class|Order|Family|Genus|Species",
    NA,
    unique(gsub(
      "Kingdom_|Phyla_|Phylum_|Class_|Order_|Family_|Genus_|Species_",
      "",
      setdiff(names(taxtab_df), "taxon")
    ))
  ) |>
    stats::na.omit() |>
    as.vector()
  labels <- grep("Taxa_name_", labels, invert = TRUE, value = TRUE)
  setdiff(labels, c("ITSx_origin", "ITSx_detected"))
}

# Species when it already carries the genus (Genus_species), else
# Genus_Species; any "NA" left by the paste becomes a real NA.
gen_sp_for <- function(tt, cl) {
  sp <- if (paste0("Species_", cl) %in% names(tt)) {
    tt[[paste0("Species_", cl)]]
  } else {
    rep(NA_character_, nrow(tt))
  }
  gen <- if (paste0("Genus_", cl) %in% names(tt)) {
    tt[[paste0("Genus_", cl)]]
  } else {
    rep(NA_character_, nrow(tt))
  }
  out <- if (sum(grepl("_", sp)) == sum(!is.na(sp))) {
    sp
  } else {
    paste(gen, sp, sep = "_")
  }
  out[grepl("NA", out)] <- NA_character_
  out
}

# Genus + Species concatenated per assignment column, as `Gen_sp_<label>`.
# The tax_table is read with tidypq::tax_table_to_df() and the new columns are
# written back in one tidypq::mutate_taxa_pq() call (decision 13).
add_gen_sp <- function(physeq) {
  taxtab_df <- tax_table_to_df(physeq, convert = FALSE)
  cl_names <- assignment_labels(taxtab_df)
  gen_sp_cols <- stats::setNames(
    lapply(cl_names, \(cl) gen_sp_for(taxtab_df, cl)),
    paste0("Gen_sp_", cl_names)
  )
  do.call(mutate_taxa_pq, c(list(physeq = physeq), gen_sp_cols))
}

# The `ranks_df` that tc_metrics_unit() takes: one row per rank of `ranks`
# (whose last element is `Gen_sp`, because the truth holds the underscored
# binomial and EUKARYOME stores the epithet alone), one column per assignment.
# An assignment with no column at all is dropped only when `allow_empty` is
# TRUE — on the mini stores, where blastn finds no hit in a mini_* database and
# combine_taxo_assignments(allow_empty = TRUE) adds no column; in production an
# empty assignment is a failure and stops.
build_ranks_df <- function(physeq, assignments, ranks, allow_empty = FALSE) {
  ranks_df <- as.data.frame(stats::setNames(
    lapply(assignments, \(a) paste0(ranks, "_", a)),
    assignments
  ))
  present <- names(tax_table_to_df(physeq, convert = FALSE))
  missing_cols <- setdiff(as.vector(unlist(ranks_df)), present)
  empty <- names(ranks_df)[vapply(
    ranks_df,
    \(cols) all(cols %in% missing_cols),
    logical(1)
  )]
  if (allow_empty && all(missing_cols %in% unlist(ranks_df[empty]))) {
    if (length(empty) > 0) {
      message(
        "mini stores: ",
        length(empty),
        " empty assignments dropped: ",
        paste(utils::head(empty, 5), collapse = ", "),
        if (length(empty) > 5) ", ..."
      )
    }
    ranks_df <- ranks_df[setdiff(names(ranks_df), empty)]
  } else if (length(missing_cols) > 0) {
    print(utils::head(missing_cols, 20))
    stop("Columns expected from values_map / consensus suffixes are missing.")
  }
  ranks_df
}
