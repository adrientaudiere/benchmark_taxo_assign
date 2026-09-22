# One computation per method × database × input, and the benchmark rows
# derived from it (ROADMAP S8.2, docs/experimental_design.md §8). The rows are
# the ones add_new_taxonomy_pq() used to compute one by one:
#
# - dada2: one assignTaxonomy(minBoot = 0, outputBootstraps = TRUE); a row
#   keeps the ranks whose bootstrap is >= 100 * min_bootstrap, which is what
#   assignTaxonomy(minBoot = ) does after computing the same bootstraps.
# - sintax: one assign_sintax(min_bootstrap = 0, behavior = "return_matrix");
#   a row sets the values whose bootstrap is < min_bootstrap to NA, as
#   assign_sintax() does.
# - lca: one vsearch search whose hits are kept
#   (assign_vsearch_lca(behavior = "return_hits")); a row computes the LCA at
#   its lca_cutoff from these hits (assign_vsearch_lca(hits_table = )), which
#   gives what vsearch --lcaout gives at that cutoff (ROADMAP 0.3).
# - blastn: one raw blast table (no score filter); a row applies the score
#   filters (min_id, min_cover) and its vote with assign_blastn(blast_table = ).
#
# The negative controls (fake_*, external_*) have no reads. assign_sintax()
# and assign_vsearch_lca() drop empty taxa before querying unless
# clean_pq = FALSE, so the controls were never assigned by these two methods
# (critique 07, B19).

compute_assignment <- function(physeq, method, ref_fasta, nproc = 1) {
  switch(
    method,
    dada2 = dada2::assignTaxonomy(
      physeq@refseq,
      refFasta = ref_fasta,
      minBoot = 0,
      outputBootstraps = TRUE,
      multithread = nproc
    ),
    sintax = MiscMetabar::assign_sintax(
      physeq,
      ref_fasta = ref_fasta,
      behavior = "return_matrix",
      min_bootstrap = 0,
      clean_pq = FALSE,
      nproc = nproc
    ),
    lca = MiscMetabar::assign_vsearch_lca(
      physeq,
      ref_fasta = ref_fasta,
      behavior = "return_hits",
      clean_pq = FALSE,
      nproc = nproc,
      verbose = FALSE
    ),
    blastn = MiscMetabar::blast_pq(
      physeq,
      fasta_for_db = ref_fasta,
      unique_per_seq = FALSE,
      score_filter = FALSE,
      nproc = nproc
    ),
    stop("Unknown method: ", method)
  )
}

# Returns `physeq` with the columns of one benchmark row (suffix
# "_<full_name>"), or `physeq` unchanged when the computation found nothing.
derive_assignment <- function(
  physeq,
  computed,
  method,
  suffix,
  min_bootstrap = NA,
  vote_algorithm = NA,
  nb_voting = NA,
  min_cover = NA,
  min_id = NA,
  lca_cutoff = NA
) {
  if (is.null(computed)) {
    return(physeq)
  }
  switch(
    method,
    dada2 = derive_dada2_row(physeq, computed, min_bootstrap, suffix),
    sintax = derive_sintax_row(physeq, computed, min_bootstrap, suffix),
    lca = MiscMetabar::assign_vsearch_lca(
      physeq,
      hits_table = computed,
      behavior = "add_to_phyloseq",
      lca_cutoff = lca_cutoff,
      suffix = suffix,
      verbose = FALSE
    ),
    blastn = MiscMetabar::assign_blastn(
      physeq,
      blast_table = computed,
      behavior = "add_to_phyloseq",
      suffix = suffix,
      vote_algorithm = vote_algorithm,
      nb_voting = nb_voting,
      min_id = min_id,
      min_cover = min_cover
    ),
    stop("Unknown method: ", method)
  )
}

# Same columns as add_new_taxonomy_pq(method = "dada2", min_bootstrap = ).
derive_dada2_row <- function(physeq, computed, min_bootstrap, suffix) {
  full <- computed$tax
  tax_tab <- matrix(
    NA_character_,
    nrow(full),
    ncol(full),
    dimnames = dimnames(full)
  )
  for (i in seq_len(nrow(full))) {
    kept <- full[i, computed$boot[i, ] >= 100 * min_bootstrap]
    if (length(kept) > 0) {
      tax_tab[i, seq_along(kept)] <- kept
    }
  }
  colnames(tax_tab) <- make.unique(paste0(colnames(tax_tab), suffix))
  new_physeq <- physeq
  phyloseq::tax_table(new_physeq) <- phyloseq::tax_table(cbind(
    physeq@tax_table,
    tax_tab
  ))
  new_physeq
}

# Same columns as add_new_taxonomy_pq(method = "sintax", min_bootstrap = ).
derive_sintax_row <- function(physeq, computed, min_bootstrap, suffix) {
  taxo <- computed$taxo_value
  taxo[computed$taxo_bootstrap < min_bootstrap] <- NA
  rank_cols <- setdiff(names(taxo), "taxa_names")
  names(taxo)[match(rank_cols, names(taxo))] <- paste0(rank_cols, suffix)

  tax_tab <- as.data.frame(as.matrix(physeq@tax_table))
  tax_tab$taxa_names <- phyloseq::taxa_names(physeq)
  new_tax_tab <- dplyr::left_join(tax_tab, taxo, by = "taxa_names") |>
    dplyr::select(-"taxa_names") |>
    as.matrix()
  new_physeq <- physeq
  new_physeq@tax_table <- phyloseq::tax_table(new_tax_tab)
  phyloseq::taxa_names(new_physeq@tax_table) <- phyloseq::taxa_names(physeq)
  new_physeq
}
