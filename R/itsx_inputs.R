# The ITSx inputs of the assignment grid (OD 34): one input per ITSx region of
# the dataset and per unit (the ASVs, and their 97 % clusters). ITSx itself is
# MiscMetabar::itsx_pq(); this file only builds the negative controls of an
# input in the order decided on 2026-09-27 (ROADMAP Q4):
#
#   1. itsx_pq() on the real units (origin recorded in `ITSx_origin`, and
#      since 2026-10-07 whether ITSx found the region in `ITSx_detected`,
#      code review §1.1);
#   2. itsx_pq() on the external controls alone, which gives their trimmed
#      sequences, their ITSx origin (needed by ROADMAP SA7) and detection;
#   3. `fake_` controls shuffled from the trimmed real units (both columns NA:
#      ITSx never sees them);
#   4. the trimmed external controls added, and their origin and detection
#      written back, since comparpq::add_external_seq_pq() blanks the
#      tax_table rows it adds.
#
# Requires MiscMetabar, comparpq, phyloseq and Biostrings.

# ITSx on the external controls (`ext_seqs`, a named DNAStringSet) alone.
# Returns the trimmed sequences (those without `region` keep their full
# sequence, like the real units), the ITSx origin of each and whether ITSx
# found `region` in it ("TRUE" / "FALSE"), all named like `ext_seqs`. Controls that become identical are merged by itsx_pq(), so the
# result can be shorter than `ext_seqs`.
itsx_controls <- function(
  ext_seqs,
  region,
  origin_col = "ITSx_origin",
  detected_col = "ITSx_detected",
  ...
) {
  otu <- phyloseq::otu_table(
    matrix(
      1,
      nrow = 1,
      ncol = length(ext_seqs),
      dimnames = list("controls", names(ext_seqs))
    ),
    taxa_are_rows = FALSE
  )
  pq <- phyloseq::phyloseq(otu, phyloseq::refseq(ext_seqs))
  trimmed <- MiscMetabar::itsx_pq(
    pq,
    region = region,
    keep_undetected = TRUE,
    add_origin = TRUE,
    origin_col = origin_col,
    add_detected = TRUE,
    detected_col = detected_col,
    duplicated_seqs = "merge",
    verbose = FALSE,
    ...
  )
  list(
    seqs = trimmed@refseq,
    origin = stats::setNames(
      as.character(trimmed@tax_table[, origin_col]),
      phyloseq::taxa_names(trimmed)
    ),
    detected = stats::setNames(
      as.character(trimmed@tax_table[, detected_col]),
      phyloseq::taxa_names(trimmed)
    )
  )
}

# The assignment input: `real` (the ITSx-trimmed real units, with their
# `origin_col`) plus `prop_fake` shuffled copies of them (`fake_`) and the
# trimmed external controls of itsx_controls() (`external_`), whose origin is
# written back into `origin_col`, and their detection into `detected_col` when
# `real` has that column. The shuffle is seeded so that every run
# draws the same fakes.
with_controls <- function(
  real,
  controls,
  prop_fake,
  seed,
  origin_col = "ITSx_origin",
  detected_col = "ITSx_detected"
) {
  shuffled <- withr::with_seed(
    seed,
    comparpq::add_shuffle_seq_pq(real, prop_fake = prop_fake)
  )
  out <- comparpq::add_external_seq_pq(shuffled, controls$seqs)
  tax <- as(out@tax_table, "matrix")
  ext_names <- paste0("external_", names(controls$origin))
  tax[ext_names, origin_col] <- controls$origin
  if (detected_col %in% colnames(tax) && !is.null(controls$detected)) {
    tax[ext_names, detected_col] <- controls$detected
  }
  out@tax_table <- phyloseq::tax_table(tax)
  out
}
