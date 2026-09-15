# Extract the suffixed columns that each assignment phyloseq added to the
# shared base tax_table, then add them all back onto the base. Replaces the
# sequential `previous_target` chain of the archived script_assign_taxo.R:
# instead of every (method, db) consuming the prior step, each one consumes
# the same base and this function reconstitutes the all-methods phyloseq.
#
# Written with tidypq (decision 13): tax_table_to_df() reads each tax_table as
# a tibble, rows are aligned on the taxon name (not assumed to be in the same
# order), and one mutate_taxa_pq() call writes every new column back.
# Requires tidypq to be loaded (load_pqverse("tidypq")).
#
# An assignment that adds no column is an error by default: MiscMetabar's
# assign_blastn() returns the input phyloseq unchanged (with only a message)
# when no BLAST hit passes the score filters, which targets records as a
# success. Without this guard the missing columns only surface much later, in
# the analysis (ROADMAP S1.3). Pass allow_empty = TRUE to tolerate it.
combine_taxo_assignments <- function(base_pq, ..., allow_empty = FALSE) {
  assignments <- list(...)
  # Label each assignment for error messages: its argument name if given,
  # else the symbol passed (tar_combine passes the target names as symbols),
  # else its position.
  dots <- match.call(expand.dots = FALSE)$...
  fallback <- vapply(seq_along(assignments), function(i) {
    if (is.symbol(dots[[i]])) {
      as.character(dots[[i]])
    } else {
      paste0("assignment ", i)
    }
  }, character(1))
  labels <- names(assignments)
  if (is.null(labels)) {
    labels <- fallback
  } else {
    labels[labels == ""] <- fallback[labels == ""]
  }

  base_df <- tidypq::tax_table_to_df(base_pq, convert = FALSE)
  base_cols <- setdiff(names(base_df), "taxon")

  added <- lapply(seq_along(assignments), function(i) {
    df <- tidypq::tax_table_to_df(assignments[[i]], convert = FALSE)
    if (!setequal(df$taxon, base_df$taxon)) {
      stop("Taxa of ", labels[i], " differ from the taxa of the base phyloseq.")
    }
    new_cols <- setdiff(names(df), c("taxon", base_cols))
    df[match(base_df$taxon, df$taxon), new_cols, drop = FALSE]
  })

  n_added <- vapply(added, ncol, integer(1))
  if (!allow_empty && any(n_added == 0)) {
    stop(
      "These assignments added no taxonomy column (assign_* returned the ",
      "input unchanged, e.g. no BLAST hit passed the score filters): ",
      paste(labels[n_added == 0], collapse = ", ")
    )
  }
  if (sum(n_added) == 0) {
    return(base_pq)
  }

  new_columns <- dplyr::bind_cols(added[n_added > 0], .name_repair = "check_unique")
  do.call(
    tidypq::mutate_taxa_pq,
    c(list(physeq = base_pq), as.list(new_columns))
  )
}
