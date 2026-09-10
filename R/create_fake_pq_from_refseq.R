# Helpers used during development of the benchmark pipeline. Sourced by
# script_dada2.R, but the pipeline does not call anything defined here at
# runtime — the helpers are kept available for interactive use.
#
# Database derivation (the fake-reference set, format conversions, Fungi-only
# filtering, mini DBs, vsearch clustering, cutadapt trimming) lives in
# make_databases.R.
#
# Cross-validation tooling (cross_val, cross_val_param) lives in
# R/cross_val.R. Worked examples live in R/examples_cross_val.R. Neither is
# sourced by the pipeline.

# Build a degenerate phyloseq from a reference fasta. If `taxonomy_in_names`
# is TRUE, the tax_table is reconstructed from the sintax-style header
# (`;tax=k:Fungi,p:Ascomycota,...`).
create_fake_pq_from_refseq <- function(references_sequences,
                                       taxonomy_in_names = TRUE,
                                       taxa_ranks = c(
                                         "Kingdom", "Phylum", "Class",
                                         "Order", "Family", "Genus", "Species"
                                       )) {
  if (is.character(references_sequences)) {
    references_sequences <- Biostrings::readDNAStringSet(references_sequences)
  }
  n_taxrank <- length(taxa_ranks)
  if (taxonomy_in_names) {
    taxtab <- stringr::str_split_fixed(names(references_sequences), ";tax=", n = 2)[, 2] |>
      stringr::str_split_fixed(",", n_taxrank)
    colnames(taxtab) <- taxa_ranks
    taxtab <- tax_table(taxtab)
  } else {
    taxtab <- tax_table(matrix(
      data = "FAKE",
      nrow = length(references_sequences),
      ncol = n_taxrank
    ))
  }

  otutab <- otu_table(matrix(
    data = 0,
    nrow = length(references_sequences),
    ncol = 1
  ), taxa_are_rows = TRUE)

  refseq <- refseq(references_sequences)

  taxa_names(otutab) <- taxa_names(refseq)
  taxa_names(taxtab) <- taxa_names(refseq)
  MiscMetabar::simplify_taxo(phyloseq(taxtab, otutab, refseq))
}
