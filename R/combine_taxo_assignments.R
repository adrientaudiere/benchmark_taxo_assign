# Extract the suffixed columns that each assignment phyloseq added to the
# shared base tax_table, then cbind them all back onto the base. Replaces the
# sequential `previous_target` chain of the archived script_assign_taxo.R:
# instead of every (method, db) consuming the prior step, each one consumes
# the same base and this function reconstitutes the all-methods phyloseq.
#
# An assignment that adds no column is an error by default: MiscMetabar's
# assign_blastn() returns the input phyloseq unchanged (with only a message)
# when no BLAST hit passes the score filters, which targets records as a
# success. Without this guard the missing columns only surface much later, in
# the analysis (ROADMAP S1.3). Pass allow_empty = TRUE to tolerate it.
combine_taxo_assignments <- function(base_pq, ..., allow_empty = FALSE) {
  assignments <- list(...)
  base_cols <- colnames(base_pq@tax_table)

  added <- lapply(assignments, function(pq) {
    tt <- as.matrix(unclass(pq@tax_table))
    tt[, setdiff(colnames(tt), base_cols), drop = FALSE]
  })

  n_added <- vapply(added, ncol, integer(1))
  if (!allow_empty && any(n_added == 0)) {
    labels <- names(assignments)
    if (is.null(labels) || all(labels == "")) {
      labels <- paste0("assignment ", seq_along(assignments))
    }
    stop(
      "These assignments added no taxonomy column (assign_* returned the ",
      "input unchanged, e.g. no BLAST hit passed the score filters): ",
      paste(labels[n_added == 0], collapse = ", ")
    )
  }

  combined_tt <- do.call(
    cbind,
    c(list(as.matrix(unclass(base_pq@tax_table))), added)
  )

  base_pq@tax_table <- phyloseq::tax_table(combined_tt)
  base_pq
}
