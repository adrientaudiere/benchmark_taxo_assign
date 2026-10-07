suppressMessages(library(targets))

# dada2 1.40.0, src/taxonomy.cpp, get_best_genus(): the k-mer profile of
# taxonomy g is read at lgk_probability[g * n_kmers] with `int g` and
# `unsigned int n_kmers` (4^8 = 65 536), so the product is computed in 32 bits
# and wraps once g >= 65 536: taxonomies g, g + 65 536, g + 131 072 ... are all
# scored with the profile of taxonomy g mod 65 536, tie exactly, and the tie is
# broken at random. The matrix itself is built with size_t, so it is correct.
#
# Zero-compute test (2026-09-30), Pauvert ITS1 ASV, stored unthresholded
# dada2 outputs (minBoot = 0): on the full EUK_ITS_v2.1 (206 399 taxonomies),
# the taxonomy whose profile dada2 actually read for its winner is
# genus.unq[(index - 1) %% 65536 + 1]. If the collapse is this bug, that
# "profile owner" is the taxonomy the same query wins on EUK_ITS_v2.1_Fungi,
# whatever lineage the reported winner carries.

refseq <- "data/data_raw/refseq/dada2_format"
headers <- names(Biostrings::fasta.seqlengths(file.path(
  refseq,
  "EUK_ITS_v2.1.fasta"
)))
# As dada2::assignTaxonomy() 1.40.0 does (lines 20-21 and 49 of its body).
headers <- gsub("^\\s+|\\s+$", "", headers)
genus_unq <- unique(headers)
message(
  "full EUK_ITS_v2.1: ",
  length(genus_unq),
  " taxonomies, ",
  round(length(genus_unq) / 65536, 2),
  " x 65 536"
)
first_kingdoms <- sub(";.*", "", genus_unq[1:65536])
message(
  "kingdoms of the first 65 536 taxonomies (the only profiles ever read):"
)
print(round(
  100 * sort(table(first_kingdoms), decreasing = TRUE)[1:5] / 65536,
  1
))

as_string <- function(tax) paste0(apply(tax, 1, paste, collapse = ";"), ";")
full <- tar_read(
  ITS1_asv_compute_dada2__EUK_ITS_v2.1,
  store = "store_assign_taxo"
)
fungi <- tar_read(
  ITS1_asv_compute_dada2__EUK_ITS_v2.1_Fungi,
  store = "store_assign_taxo"
)
stopifnot(identical(rownames(full$tax), rownames(fungi$tax)))

winner <- as_string(full$tax)
idx <- match(winner, genus_unq)
message(
  "winners found among the taxonomies: ",
  sum(!is.na(idx)),
  " of ",
  length(idx)
)
owner <- genus_unq[(idx - 1) %% 65536 + 1]
fungi_winner <- as_string(fungi$tax)

res <- data.frame(
  above = idx > 65536,
  winner_kingdom = sub(";.*", "", winner),
  owner_is_fungi_winner = owner == fungi_winner,
  winner_is_fungi_winner = winner == fungi_winner,
  boot_kingdom = full$boot[, 1],
  boot_genus = full$boot[, 6]
)
message(
  "\nqueries whose winner index is above 65 536: ",
  sum(res$above, na.rm = TRUE),
  " of ",
  nrow(res)
)
message(
  "reported winner = the query's winner on EUK_ITS_v2.1_Fungi: ",
  sum(res$winner_is_fungi_winner, na.rm = TRUE)
)
message(
  "profile owner (index mod 65 536) = the query's winner on EUK_ITS_v2.1_Fungi: ",
  sum(res$owner_is_fungi_winner, na.rm = TRUE)
)
message("\nreported kingdom of the winners:")
print(table(res$winner_kingdom))
message(
  "\nmedian bootstrap on the full release, Kingdom / Genus: ",
  median(res$boot_kingdom),
  " / ",
  median(res$boot_genus)
)

# Number of taxonomies sharing each profile read (the aliases), and the
# bootstrap the bug predicts: about 1 / aliases at a rank where they differ.
n_alias <- tabulate((seq_along(genus_unq) - 1) %% 65536 + 1, nbins = 65536)
owner_idx <- (idx - 1) %% 65536 + 1
res$aliases <- n_alias[owner_idx]
message("\naliases per profile read, and median Kingdom bootstrap:")
print(aggregate(boot_kingdom ~ aliases, data = res, FUN = median))
print(table(aliases = res$aliases))

# The same counts per unit class (the rows of the compute are sequences).
d <- tar_read(d_all_taxo_ITS1_asv, store = "store_assign_taxo")
seq_to_name <- setNames(names(d@refseq), as.character(d@refseq))
unit <- seq_to_name[rownames(full$tax)]
stopifnot(!anyNA(unit))
res$class <- ifelse(
  grepl("^external_", unit),
  "external",
  ifelse(grepl("^fake_", unit), "fake", "real")
)
message(
  "\nper class: n, reported winner = Fungi winner, profile owner = Fungi winner"
)
print(aggregate(
  cbind(n = 1, winner_is_fungi_winner, owner_is_fungi_winner) ~ class,
  data = res,
  FUN = sum
))
