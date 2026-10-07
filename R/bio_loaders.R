# Loaders of the D1c datasets that publish their OTUs with representative
# sequences, so they skip dada2 (OD 30): the published OTUs play the role of
# `d_asv` in store_dada2_<dataset>, and assign_taxo_<dataset> assigns them.
#
# The published taxonomy goes to the tax_table under plain rank names, like the
# seed taxonomy of pipelines/dada2_bio.R: add_shuffle_seq_pq() needs a
# tax_table, and assignment_labels() (R/assignment_helpers.R) ignores
# rank-only column names, so it is never scored as an assignment.
#
# Requires phyloseq, Biostrings and data.table.

# Units kept from a published table (OD 30b, revised 2026-09-26): the common
# ones, those occurring in at least `min_occur` samples of `otu` (samples in
# rows), returned sorted. The benchmark measures assignment on the OTUs a study
# actually relies on, not on a random draw that is half singletons.
common_units <- function(otu, min_occur = 1) {
  sort(colnames(otu)[colSums(otu > 0) >= min_occur])
}

# The published taxonomy of a deposit as a character matrix, one row per OTU
# and one column per rank, with the "k__" / "p__" prefixes some deposits keep
# stripped so both deposits read the same way. Empty strings become NA.
published_taxonomy <- function(taxo, id_col, rank_cols) {
  tax <- as.matrix(taxo[, rank_cols])
  tax <- sub("^[kpcofgs]__", "", tax)
  tax[tax == ""] <- NA_character_
  dimnames(tax) <- list(
    taxo[[id_col]],
    c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  )
  tax
}

# Korhonen et al. 2024, two Dryad deposits from the same lab (ITS2, Illumina,
# presence / absence): 10.5061/dryad.3j9kd51s7 (deciduous fine woody debris,
# 100 replicates, plus 20 birch logs) and 10.5061/dryad.44j0zpcng (spruce
# trunks, tops and branches, 93 samples). `dir` holds one sub-folder per
# deposit, named after its DOI as Dryad serves it. Pooled as one dataset (OD
# 30): the OTUs are keyed by the authors' hash `OTU_ID`, which both deposits
# share, so an OTU found in both is one taxon. The OTUs occurring in at least
# `min_occur` samples of the pooled table are kept (the default 1 keeps every
# OTU present somewhere: the deciduous tables also list 155 OTUs present
# nowhere, 4 of them missing from its taxonomy); every published sample is
# kept, the replicates the authors flag as outliers with `outlier = TRUE`.
# 14 RepSeq carry lowercase (soft-masked) stretches; DNAStringSet() stores
# them in upper case, the only difference from the deposit.
load_korhonen_wood <- function(dir, min_occur = 1) {
  deposit_dir <- function(doi) {
    found <- list.files(dir, pattern = doi, full.names = TRUE)
    if (length(found) != 1) {
      stop("load_korhonen_wood(): expected one folder matching ", doi, ".")
    }
    found
  }
  read_semicolon <- function(folder, file) {
    data.table::fread(
      file.path(folder, paste0("Korhonen_et_al_2024_", file, ".csv")),
      sep = ";",
      dec = ",",
      encoding = "Latin-1",
      data.table = FALSE
    )
  }
  presence <- function(tab, id_col, sample_ids, otu_ids) {
    mat <- as.matrix(tab[, setdiff(names(tab), id_col), drop = FALSE])
    dimnames(mat) <- list(sample_ids, otu_ids[colnames(mat)])
    mat[, colSums(mat) > 0, drop = FALSE]
  }
  ranks <- c("k", "p", "c", "o", "f", "g", "s")

  # Deciduous deposit: the tables name OTUs by `OTU_name`, the taxonomy maps
  # it to the hash (written there with a leading ">").
  dec_dir <- deposit_dir("3j9kd51s7")
  dec_taxo <- read_semicolon(dec_dir, "Fungi_ITS2_OTU_taxonomy")
  dec_taxo$OTU_ID <- sub("^>", "", dec_taxo$OTU_ID)
  dec_hash <- stats::setNames(dec_taxo$OTU_ID, dec_taxo$OTU_name)
  dfwd <- read_semicolon(dec_dir, "Fungi_presabs_DFWD_replicates")
  dfwd_meta <- read_semicolon(dec_dir, "Fungi_DFWD_SampleNo_SiteNo_table")
  birch <- read_semicolon(dec_dir, "Fungi_presabs_BirchLog_samples")
  birch_meta <- read_semicolon(dec_dir, "BirchLog_sample_metadata")

  # Spruce deposit: the tables name OTUs by the hash directly.
  spr_dir <- deposit_dir("44j0zpcng")
  spr_taxo <- read_semicolon(spr_dir, "Fungal_OTU_taxonomy_and_metadata")
  spruce <- read_semicolon(spr_dir, "Fungal_OTU_presabs")

  dfwd_ids <- sprintf("DFWD_%03d", dfwd$SampleNo)
  birch_ids <- sprintf("BirchLog_%02d", birch$BirchLogNo)
  tables <- list(
    presence(dfwd, "SampleNo", dfwd_ids, dec_hash),
    presence(birch, "BirchLogNo", birch_ids, dec_hash),
    presence(
      spruce,
      "SampleID",
      spruce$SampleID,
      stats::setNames(spr_taxo$OTU_ID, spr_taxo$OTU_ID)
    )
  )

  keep <- sort(unique(unlist(lapply(tables, colnames))))
  otu <- matrix(
    0L,
    nrow = sum(vapply(tables, nrow, integer(1))),
    ncol = length(keep),
    dimnames = list(unlist(lapply(tables, rownames)), keep)
  )
  for (tab in tables) {
    cols <- intersect(colnames(tab), keep)
    otu[rownames(tab), cols] <- tab[, cols]
  }
  keep <- common_units(otu, min_occur)
  otu <- otu[, keep, drop = FALSE]
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]

  # Taxonomy and sequence of each kept OTU, from the deciduous deposit when it
  # holds the OTU (the 84 shared OTUs carry the same sequence in both).
  tax_rows <- rbind(
    dec_taxo[, c("OTU_ID", "RepSeq", ranks)],
    spr_taxo[
      !spr_taxo$OTU_ID %in% dec_taxo$OTU_ID,
      c("OTU_ID", "RepSeq", ranks)
    ]
  )
  tax_rows <- tax_rows[match(keep, tax_rows$OTU_ID), ]

  outliers <- dfwd_meta$SampleNo[dfwd_meta[["Valid/Outlier"]] == "outlier"]
  meta <- rbind(
    data.frame(
      deposit = "3j9kd51s7",
      substrate = "deciduous_fwd",
      site = paste0(
        "site_",
        dfwd_meta$SiteNo[match(dfwd$SampleNo, dfwd_meta$SampleNo)]
      ),
      outlier = dfwd$SampleNo %in% outliers,
      row.names = dfwd_ids
    ),
    data.frame(
      deposit = "3j9kd51s7",
      substrate = "birch_log",
      site = paste0(
        "site_",
        birch_meta$SiteNo[match(birch$BirchLogNo, birch_meta$BirchLogNo)]
      ),
      outlier = FALSE,
      row.names = birch_ids
    ),
    data.frame(
      deposit = "44j0zpcng",
      substrate = paste0(
        "spruce_",
        tolower(sub("^Fungi_([A-Za-z]+)_.*$", "\\1", spruce$SampleID))
      ),
      site = paste0("tree_", as.integer(sub("^.*_", "", spruce$SampleID))),
      outlier = FALSE,
      row.names = spruce$SampleID
    )
  )
  meta <- meta[rownames(otu), ]

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = FALSE),
    phyloseq::sample_data(meta),
    phyloseq::tax_table(published_taxonomy(tax_rows, "OTU_ID", ranks)),
    Biostrings::DNAStringSet(stats::setNames(tax_rows$RepSeq, tax_rows$OTU_ID))
  )
}

# Global Spore Sampling Project (Ovaskainen et al. 2024, Zenodo 11125610): air,
# ITS2, 27 954 species-level OTUs, all published as Fungi. `dir` holds
# taxonomy.csv (OTU, taxonomy and sequence), otu.table.csv (samples in rows)
# and metadata.csv, all semicolon-separated with decimal commas. The OTUs
# occurring in at least `min_occur` samples are kept (10 of 2 768 in the
# benchmark: 2 590 OTUs); samples left empty are dropped.
load_gssp_air <- function(dir, min_occur = 1) {
  taxo <- data.table::fread(
    file.path(dir, "taxonomy.csv"),
    sep = ";",
    data.table = FALSE
  )
  counts <- data.table::fread(
    file.path(dir, "otu.table.csv"),
    sep = ";",
    data.table = FALSE
  )
  otu <- as.matrix(counts[, -1])
  rownames(otu) <- counts$V1
  keep <- common_units(otu, min_occur)
  otu <- otu[, keep, drop = FALSE]
  otu <- otu[rowSums(otu) > 0, , drop = FALSE]
  taxo <- taxo[match(keep, taxo$OTU), ]

  meta <- data.table::fread(
    file.path(dir, "metadata.csv"),
    sep = ";",
    dec = ",",
    data.table = FALSE
  )
  rownames(meta) <- meta$sample.id
  meta <- meta[rownames(otu), ]

  ranks <- c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
  tax <- as.matrix(taxo[, ranks])
  dimnames(tax) <- list(taxo$OTU, tools::toTitleCase(ranks))
  refseq <- Biostrings::DNAStringSet(stats::setNames(taxo$sequence, taxo$OTU))

  phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = FALSE),
    phyloseq::sample_data(meta),
    phyloseq::tax_table(tax),
    refseq
  )
}
