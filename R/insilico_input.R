# Input records of the in silico set (D1b, OD 29): which UNITE records
# InSilicoSeq generates the reads from.
#
# The set is a positive control, so its truth must be fully nameable and every
# record must survive the pipeline. Three filters, each fixing a defect measured
# on the first generation of 2026-09-16 (100 records drawn at random, archived
# in Archives/insilico_gen1_20260916/):
#
#   - no `Incertae_sedis` / `unclassified` label at any rank: 25 of its 48
#     recovered ASVs carried one, which clean_taxo() turns into NA on the
#     assignment side, so a method giving the exact UNITE label was scored
#     as missing it;
#   - a real species epithet: 70 of its 100 records were "<Genus>_sp", which
#     species_name() refuses as a species truth, so only 13 of the 48 ASVs
#     could be scored at Species;
#   - an amplicon length the 2 x 300 bp MiSeq reads can merge: 24 abundant
#     records exceeded the ~588 bp merge limit and vanished, taking 22 genera
#     with them. `max_len` leaves a 40 bp overlap; `min_len` keeps templates
#     at least one read long.
#
# The cost, to state wherever the set is reported: the hard, poorly known
# lineages are excluded by construction (Rozellomycota and the placeholder
# phyla have no record without an incertae sedis rank).
#
# Requires species_name() from R/mock_truth.R, and Biostrings.

# Records of `seqs` (a named DNAStringSet in the UNITE general-release header
# format, "Name|acc|SH|type|k__...;p__...;...;s__Genus_species") that pass the
# three filters, one per species, drawn round-robin over phyla so that small
# phyla are fully represented before large ones get a second record. Returns
# the selected DNAStringSet; the draw depends on `seed` only.
select_insilico_records <- function(
  seqs,
  n = 100,
  min_len = 300,
  max_len = 560,
  seed = 22
) {
  headers <- names(seqs)
  species <- species_name(sub("^.*s__([^;|]*).*$", "\\1", headers))
  keep <- grepl("k__Fungi", headers, fixed = TRUE) &
    !grepl("ncertae_sedis|unclassified", headers, ignore.case = TRUE) &
    !is.na(species) &
    Biostrings::width(seqs) >= min_len &
    Biostrings::width(seqs) <= max_len
  if (!any(keep)) {
    stop("select_insilico_records(): no record passes the filters.")
  }
  pool <- data.frame(
    index = which(keep),
    species = species[keep],
    phylum = sub("^.*p__([^;|]*).*$", "\\1", headers[keep])
  )

  set.seed(seed)
  pool <- pool[sample(nrow(pool)), ]
  pool <- pool[!duplicated(pool$species), ]

  # Round-robin: the k-th record of every phylum comes before the (k+1)-th of
  # any, so the phyla are balanced as far as their species allow.
  pool$turn <- stats::ave(seq_len(nrow(pool)), pool$phylum, FUN = seq_along)
  pool <- pool[order(pool$turn, pool$phylum), ]
  picked <- utils::head(pool$index, n)
  if (length(picked) < n) {
    message(
      "select_insilico_records(): only ",
      length(picked),
      " species pass the filters; returning them all."
    )
  }
  seqs[picked]
}
