# ITS extraction of the ASV sequences with ITSx (ROADMAP Q4). The ITS1F-ITS2
# ASVs start with ~45 bp of the 18S end and stop in the 5.8S; ITS references
# start at ITS1, so the flank lowers the BLAST query cover (critique 07, B18)
# and may bias the other methods. No pqverse wrapper exists yet: candidate for
# MiscMetabar, next to cutadapt_remove_primers().

# Run ITSx on `seqs` (a named DNAStringSet) and return the extracted `region`
# (named like `seqs`) and the ITSx positions table. `prelude` activates the
# conda env holding ITSx (config.R::itsx_conda_prelude).
run_itsx <- function(seqs,
                     region = "ITS1",
                     organism_groups = "F",
                     cpu = 1,
                     prelude = "") {
  work_dir <- tempfile("itsx_")
  dir.create(work_dir)
  on.exit(unlink(work_dir, recursive = TRUE), add = TRUE)
  input <- file.path(work_dir, "input.fasta")
  prefix <- file.path(work_dir, "itsx")
  Biostrings::writeXStringSet(seqs, input, width = 20000)

  cmd <- sprintf(
    "%sITSx -i %s -o %s -t %s --cpu %d --preserve T --save_regions %s --graphical F --silent T",
    prelude, shQuote(input), shQuote(prefix), organism_groups, as.integer(cpu), region
  )
  status <- system2("bash", c("-c", shQuote(cmd)), stdout = FALSE, stderr = FALSE)
  if (status != 0) {
    stop("ITSx failed (exit ", status, "): ", cmd)
  }

  region_file <- paste0(prefix, ".", region, ".fasta")
  extracted <- if (file.exists(region_file) && file.size(region_file) > 0) {
    Biostrings::readDNAStringSet(region_file)
  } else {
    Biostrings::DNAStringSet()
  }
  names(extracted) <- sub("\\s.*$", "", names(extracted))

  positions_file <- paste0(prefix, ".positions.txt")
  positions <- if (file.exists(positions_file) && file.size(positions_file) > 0) {
    utils::read.delim(
      positions_file,
      header = FALSE,
      quote = "",
      fill = TRUE,
      col.names = c("taxon", "length", "SSU", "ITS1", "5.8S", "ITS2", "LSU", "comment"),
      check.names = FALSE
    ) |>
      tibble::as_tibble()
  } else {
    tibble::tibble(taxon = character())
  }

  list(sequences = extracted, positions = positions, region = region)
}

# Replace the refseq of the taxa that ITSx detected by the extracted region.
# Undetected taxa keep their full sequence (keep_undetected = TRUE, so the
# taxa set, and therefore the negative controls and the metrics denominators,
# stay those of `physeq`) or are removed.
itsx_replace_refseq <- function(physeq, itsx, keep_undetected = TRUE) {
  seqs <- physeq@refseq
  detected <- names(seqs) %in% names(itsx$sequences)
  if (any(!detected)) {
    message(
      sum(!detected), " of ", length(seqs), " taxa without ", itsx$region,
      " detected by ITSx: ",
      if (keep_undetected) "kept with their full sequence" else "removed"
    )
  }
  if (!keep_undetected && any(!detected)) {
    physeq <- phyloseq::prune_taxa(names(seqs)[detected], physeq)
    seqs <- physeq@refseq
    detected <- rep(TRUE, length(seqs))
  }
  new_seqs <- as.character(seqs)
  new_seqs[detected] <- as.character(itsx$sequences[names(seqs)[detected]])
  physeq@refseq <- Biostrings::DNAStringSet(stats::setNames(new_seqs, names(seqs)))

  # Two ASVs differing only outside `itsx$region` (e.g. in the 18S/5.8S flank)
  # become identical once trimmed, and verify_pq() (run by every assign_*())
  # rejects duplicated refseq entries. They must be removed beforehand, from
  # both inputs, with itsx_duplicated_taxa() (ROADMAP decision 20).
  if (anyDuplicated(new_seqs)) {
    stop(
      sum(duplicated(new_seqs)), " taxa share their ", itsx$region, " sequence ",
      "with another taxon: drop itsx_duplicated_taxa() from `physeq` first."
    )
  }
  physeq
}

# Taxa to remove so that no two taxa share a sequence once trimmed to
# `itsx$region` (undetected taxa count with their full sequence). In each
# group of identical trimmed sequences the most abundant taxon is kept (the
# first one on ties) and the others are returned. Removing them from the raw
# ASVs too keeps the same taxa in both inputs (ROADMAP decision 20); the rule
# applies to any dataset (mock, in silico, ...).
itsx_duplicated_taxa <- function(physeq, itsx) {
  seqs <- as.character(physeq@refseq)
  trimmed <- seqs
  hit <- names(seqs) %in% names(itsx$sequences)
  trimmed[hit] <- as.character(itsx$sequences[names(seqs)[hit]])
  abundance <- phyloseq::taxa_sums(physeq)[names(seqs)]
  by_abundance <- order(-abundance, seq_along(seqs))
  dropped <- names(seqs)[by_abundance][duplicated(trimmed[by_abundance])]
  if (length(dropped) > 0) {
    message(
      length(dropped), " taxa removed from both inputs (same ", itsx$region,
      " as a more abundant taxon): ", paste(dropped, collapse = ", ")
    )
  }
  dropped
}
