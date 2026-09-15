# Build a degenerate phyloseq (one zero-count sample) from a set of reference
# sequences, so that the assign_*() functions can consume one cross-validation
# fold. Sourced by pipelines/cross_val.R; used by R/cross_val.R.
#
# If `taxonomy_in_names` is TRUE, the tax_table is parsed from the fasta
# headers. Both header styles of data/data_raw/refseq/ are understood:
#   sintax_format: "ID;tax=k:Fungi,p:Ascomycota,..."
#   dada2_format:  "Name|Acc|SH|type|k__Fungi;p__Ascomycota;..." (UNITE)
#                  "Fungi;Ascomycota;...;Species;"                 (EUKARYOME)
# Before 2026-09-10 only the sintax style was parsed, so every dada2-format
# fold got an empty truth table and dada2 scored 0 in cross-validation
# (ROADMAP S1.4). Rank prefixes (k__, k:) are removed by simplify_taxo().
# `headers`, when given (one per sequence, same order), carries the taxonomy
# instead of the names: cross_val() names dada2-format records by their
# sintax-format identifiers and keeps the dada2 headers aside (ROADMAP B22).
create_fake_pq_from_refseq <- function(references_sequences,
                                       taxonomy_in_names = TRUE,
                                       taxa_ranks = c(
                                         "Kingdom", "Phylum", "Class",
                                         "Order", "Family", "Genus", "Species"
                                       ),
                                       headers = NULL) {
  if (is.character(references_sequences)) {
    references_sequences <- Biostrings::readDNAStringSet(references_sequences)
  }
  n_taxrank <- length(taxa_ranks)
  if (taxonomy_in_names) {
    if (is.null(headers)) {
      headers <- names(references_sequences)
    }
    taxtab <- parse_header_taxonomy(headers, n_taxrank)
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

# Character matrix (one row per header, `n_taxrank` columns) of the taxonomy
# carried by fasta headers in sintax or dada2 style. Missing trailing ranks
# are "".
parse_header_taxonomy <- function(headers, n_taxrank) {
  is_sintax <- grepl(";tax=", headers, fixed = TRUE)
  # dada2 style: keep what follows the last "|" (UNITE identifiers), drop a
  # trailing ";" (EUKARYOME), then use "," as separator like sintax.
  dada2_strings <- gsub(";", ",", sub(";$", "", sub("^.*\\|", "", headers)))
  tax_strings <- ifelse(is_sintax, sub("^.*;tax=", "", headers), dada2_strings)
  stringr::str_split_fixed(tax_strings, ",", n_taxrank)
}
