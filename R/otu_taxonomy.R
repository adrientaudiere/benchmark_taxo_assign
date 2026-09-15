# OTU taxonomy built from the member ASVs (ROADMAP Q5, decided 2026-09-15).
# The ASVs are clustered once at 97 % with vsearch, as
# MiscMetabar::postcluster_pq(method = "vsearch"). Each assignment column of an
# OTU then keeps a value only when all its member ASVs assigned at that rank
# agree (MiscMetabar::resolve_vector_ranks(method = "unanimity"), NA ignored).
# No assignment is recomputed: the OTU results are derived from d_all_taxo
# (analysis/01_load_and_clean.qmd).
#
# Requires phyloseq, dplyr, tibble and the MiscMetabar and tidypq checkouts.

# Cluster membership of the taxa of `physeq`. Returns a tibble with one row per
# taxon: `taxon`, `cluster` (vsearch cluster number) and `archetype` (most
# abundant member, the OTU name kept by merge_taxa_vec()).
otu_membership <- function(physeq, id = 0.97, ...) {
  uc <- MiscMetabar::postcluster_pq(
    dna_seq = stats::setNames(
      as.character(physeq@refseq),
      phyloseq::taxa_names(physeq)
    ),
    method = "vsearch",
    id = id,
    ...
  )
  hits <- uc[uc$type != "C", c("query", "cluster")]
  if (!setequal(hits$query, phyloseq::taxa_names(physeq))) {
    stop("vsearch clusters do not cover exactly the taxa of `physeq`.")
  }
  reads <- phyloseq::taxa_sums(physeq)
  tibble::tibble(taxon = hits$query, cluster = hits$cluster) |>
    dplyr::mutate(reads = reads[taxon]) |>
    dplyr::group_by(cluster) |>
    dplyr::mutate(archetype = taxon[which.max(reads)]) |>
    dplyr::ungroup() |>
    dplyr::select(taxon, cluster, archetype)
}

# Merge the taxa of `physeq` into the OTUs of `membership` (otu_membership()):
# counts are summed and the refseq of the most abundant member is kept
# (merge_taxa_vec(tax_adjust = 0)), then every tax_table column of a
# multi-member OTU is resolved across its members with resolve_vector_ranks().
# Taxa absent from `membership` (the negative controls) are kept unchanged.
resolve_otu_taxonomy <- function(physeq, membership, method = "unanimity") {
  if (!all(membership$taxon %in% phyloseq::taxa_names(physeq))) {
    stop("Some taxa of `membership` are not in `physeq`.")
  }
  group <- stats::setNames(
    phyloseq::taxa_names(physeq),
    phyloseq::taxa_names(physeq)
  )
  group[membership$taxon] <- paste0("otu_", membership$cluster)

  tt <- tidypq::tax_table_to_df(physeq, convert = FALSE)
  tax_cols <- setdiff(names(tt), "taxon")
  multi <- names(which(table(group) > 1))
  resolved <- lapply(multi, function(g) {
    members <- tt[group[tt$taxon] == g, tax_cols, drop = FALSE]
    vapply(
      members,
      \(vec) {
        as.character(MiscMetabar::resolve_vector_ranks(vec, method = method))
      },
      character(1)
    )
  })
  names(resolved) <- multi

  otu <- MiscMetabar::merge_taxa_vec(
    physeq,
    group = unname(group),
    tax_adjust = 0L
  )
  otu_tt <- tidypq::tax_table_to_df(otu, convert = FALSE)
  otu_group <- group[otu_tt$taxon]
  for (g in multi) {
    otu_tt[otu_group == g, tax_cols] <- as.list(resolved[[g]])
  }
  do.call(
    tidypq::mutate_taxa_pq,
    c(list(physeq = otu), as.list(otu_tt[tax_cols]))
  )
}
