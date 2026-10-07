source("analysis/00_setup.R")
source(here("R/mock_truth.R"))

# Compare the outputs of analysis/sandbox/run_dada2_patched.R (installed dada2
# 1.40.0 against the build patched at src/taxonomy.cpp line 85) on the Pauvert
# ITS1 ASV input. Rates on the real units at minBoot = 50 (the pipeline's
# default 0.5): share given a genus, share given the true genus (Genus or
# Genus_accepted of the mock truth, among the units with a genus truth), and
# the median bootstrap at Kingdom and Genus; plus the controls given a genus.
#   Rscript analysis/sandbox/compare_dada2_patched.R <label=file.rds> ...

args <- commandArgs(trailingOnly = TRUE)
runs <- setNames(sub("^[^=]+=", "", args), sub("=.*$", "", args))

d <- tar_read(d_all_taxo_ITS1_asv, store = store_assign_taxo)
seq_to_name <- setNames(names(d@refseq), as.character(d@refseq))
cfg <- mock_sanger_tables[mock_sanger_tables$dataset == "pauvert", ]
truth <- build_mock_truth(tar_read(d_asv, store = here("store_dada2")), cfg)

summarise_run <- function(label, path) {
  res <- readRDS(path)
  unit <- seq_to_name[rownames(res$tax)]
  stopifnot(!anyNA(unit))
  class <- ifelse(
    grepl("^external_", unit),
    "external",
    ifelse(grepl("^fake_", unit), "fake", "real")
  )
  genus <- res$tax[, 6]
  genus[res$boot[, 6] < 50] <- NA
  genus[grepl("unclassified|Incertae_sedis|_DADA2_UNSPECIFIED", genus)] <- NA
  tr <- truth[match(unit, truth$taxon), ]
  has_truth <- class == "real" & !is.na(tr$Genus)
  correct <- !is.na(genus) &
    (genus == tr$Genus |
      genus %in% tr$Genus_accepted & genus == tr$Genus_accepted)
  tibble(
    run = label,
    real_genus = mean(!is.na(genus[class == "real"])),
    real_correct = sum(correct & has_truth) / sum(has_truth),
    boot_kingdom = median(res$boot[class == "real", 1]),
    boot_genus = median(res$boot[class == "real", 6]),
    fake_genus = mean(!is.na(genus[class == "fake"])),
    ext_genus = mean(!is.na(genus[class == "external"])),
    kingdom_fungi = mean(res$tax[class == "real", 1] %in% "Fungi")
  )
}

out <- bind_rows(Map(summarise_run, names(runs), runs))
out |>
  mutate(across(where(is.numeric), \(v) round(v, 3))) |>
  as.data.frame() |>
  print(row.names = FALSE)
