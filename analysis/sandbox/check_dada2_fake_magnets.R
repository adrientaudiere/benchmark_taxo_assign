suppressMessages({
  library(targets)
  library(phyloseq)
})

# Since the dada2 overflow fix (2026-10-01), dada2 on the full EUK_ITS_v2.1
# gives a genus (bootstrap >= 50) to 55 of the 97 `fake_` controls of Pauvert
# ITS1 ASV, against 0 on its `_Fungi` variants, whose records are capped at
# 5 000 bp (OD 32). Are the winning taxonomies those of long records, the
# "magnets" of sintax and the Laccaria clone (SA5)? Patched stored output.

full <- tar_read(
  ITS1_asv_compute_dada2__EUK_ITS_v2.1,
  store = "store_assign_taxo"
)
d <- tar_read(d_all_taxo_ITS1_asv, store = "store_assign_taxo")
unit <- setNames(names(d@refseq), as.character(d@refseq))[rownames(full$tax)]
stopifnot(!anyNA(unit))
class <- ifelse(
  grepl("^external_", unit),
  "external",
  ifelse(grepl("^fake_", unit), "fake", "real")
)

len <- Biostrings::fasta.seqlengths(
  "data/data_raw/refseq/dada2_format/EUK_ITS_v2.1.fasta"
)
tax <- gsub("^\\s+|\\s+$", "", names(len))
# Longest record carrying each taxonomy string (dada2's unit of assignment).
max_len <- tapply(len, tax, max)
n_rec <- tapply(len, tax, length)

winner <- paste0(apply(full$tax, 1, paste, collapse = ";"), ";")
stopifnot(all(winner %in% names(max_len)))
named <- full$boot[, "Genus"] >= 50

summary_of <- function(keep, label) {
  w <- winner[keep]
  data.frame(
    group = label,
    n = sum(keep),
    median_longest_record = median(max_len[w]),
    share_over_5kb = round(mean(max_len[w] > 5000), 2),
    median_records_in_taxonomy = median(n_rec[w]),
    kingdom_fungi = round(mean(sub(";.*", "", w) == "Fungi"), 2)
  )
}
print(rbind(
  summary_of(class == "fake" & named, "fake_ given a genus"),
  summary_of(class == "fake" & !named, "fake_ left NA (winner below 50)"),
  summary_of(class == "real" & named, "real given a genus")
))
message("\nwinning taxonomies of the named fake_, most frequent:")
print(head(sort(table(winner[class == "fake" & named]), decreasing = TRUE), 5))
message(
  "\nreference: ",
  round(100 * mean(len > 5000), 2),
  " % of the full EUK_ITS_v2.1 records are over 5 000 bp"
)
