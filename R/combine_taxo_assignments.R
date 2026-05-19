# Extract the suffixed columns that each assignment phyloseq added to the
# shared base tax_table, then cbind them all back onto the base. Replaces the
# sequential `previous_target` chain in the original script_assign_taxo.R:
# instead of every (method, db) consuming the prior step, each one consumes
# the same base and this function reconstitutes the all-methods phyloseq.
combine_taxo_assignments <- function(base_pq, ...) {
  assignments <- list(...)
  base_cols <- colnames(base_pq@tax_table)

  added <- lapply(assignments, function(pq) {
    tt <- as.matrix(unclass(pq@tax_table))
    tt[, setdiff(colnames(tt), base_cols), drop = FALSE]
  })

  combined_tt <- do.call(
    cbind,
    c(list(as.matrix(unclass(base_pq@tax_table))), added)
  )

  base_pq@tax_table <- phyloseq::tax_table(combined_tt)
  base_pq
}
