suppressMessages({
  library(targets)
  library(dplyr)
  library(phyloseq)
})
source("config.R")

# Q4 follow-up (ROADMAP §2, 2026-10-05): can blastn return to min_cover = 95
# on the ITSx inputs? Decision 18 lowered it to 80 because the raw ITS1F-ITS2
# ASVs kept 30-45 bp of 18S that ITS references do not contain; ITSx removes
# those flanks since OD 34. No BLAST is rerun: each computation stores its hit
# table (`<input>_compute_blastn__<db>`, with the query cover), and the default
# row (rel_majority, nb_voting 100, min_id 95) is derived again at
# min_cover 80 and 95 with assign_blastn(blast_table = ), as R/assign_compute.R
# does. Real units only (the controls are left out), every mock view.
#
#   Rscript analysis/sandbox/check_blastn_min_cover_itsx.R

views <- tibble::tribble(
  ~dataset            , ~project                        , ~region ,
  "pauvert"           , "assign_taxo"                   , "ITS1"  ,
  "tedersoo_illumina" , "assign_taxo_tedersoo_illumina" , "ITS1"  ,
  "tedersoo_illumina" , "assign_taxo_tedersoo_illumina" , "ITS2"  ,
  "tedersoo_illumina" , "assign_taxo_tedersoo_illumina" , "full"
)
covers <- c(80, 95)
dada2_stores <- c(
  pauvert = "store_dada2",
  tedersoo_illumina = "store_dada2_tedersoo_illumina"
)
# Per-unit truth as in chapter 01 (Genus, or its GBIF accepted name).
source("R/mock_truth.R")
truths <- lapply(names(dada2_stores), function(ds) {
  cfg <- mock_sanger_tables[mock_sanger_tables$dataset == ds, ]
  build_mock_truth(tar_read(d_asv, store = dada2_stores[[ds]]), cfg)
})
names(truths) <- names(dada2_stores)
control_pattern <- "^(fake_|external_)"

genus_at <- function(physeq, hits, cover) {
  out <- MiscMetabar::assign_blastn(
    physeq,
    blast_table = hits,
    behavior = "add_to_phyloseq",
    suffix = "_check",
    vote_algorithm = "rel_majority",
    nb_voting = 100,
    min_id = 95,
    min_cover = cover
  )
  tt <- as.data.frame(tax_table(out), stringsAsFactors = FALSE)
  setNames(tt[["Genus_check"]], rownames(tt))
}

check_view <- function(v) {
  Sys.setenv(TAR_PROJECT = v$project)
  input <- paste0(v$region, "_asv")
  physeq <- tar_read_raw(paste0("d_", input, "_for_assignation"))
  real <- taxa_names(physeq)[!grepl(control_pattern, taxa_names(physeq))]
  bind_rows(lapply(db_list, function(db) {
    hits <- tar_read_raw(paste0(input, "_compute_blastn__", db))
    real_hits <- hits[hits[["Query name"]] %in% real, ]
    # Query cover of the best hit (highest bit score) of each real unit.
    best_cover <- real_hits |>
      group_by(`Query name`) |>
      slice_max(`bit score`, n = 1, with_ties = FALSE) |>
      pull(`Query cover`)
    g <- lapply(covers, \(cv) genus_at(physeq, hits, cv)[real])
    g80 <- g[[1]]
    g95 <- g[[2]]
    truth <- truths[[v$dataset]]
    tr <- truth[match(real, truth$taxon), ]
    has_truth <- !is.na(tr$Genus) & !tr$foreign
    right_80 <- !is.na(g80) &
      has_truth &
      (g80 == tr$Genus | (!is.na(tr$Genus_accepted) & g80 == tr$Genus_accepted))
    is_lost <- !is.na(g80) & is.na(g95)
    tibble(
      dataset = v$dataset,
      region = v$region,
      db = db,
      n_real = length(real),
      best_cover_ge95 = mean(best_cover >= 95),
      best_cover_median = median(best_cover),
      genus_80 = sum(!is.na(g80)),
      genus_95 = sum(!is.na(g95)),
      lost = sum(is_lost),
      lost_right = sum(is_lost & right_80),
      lost_wrong = sum(is_lost & has_truth & !right_80),
      lost_no_truth = sum(is_lost & !has_truth),
      gained = sum(is.na(g80) & !is.na(g95)),
      changed = sum(!is.na(g80) & !is.na(g95) & g80 != g95)
    )
  }))
}

res <- bind_rows(lapply(seq_len(nrow(views)), \(i) check_view(views[i, ])))
print(as.data.frame(res), digits = 3)

message("Summary per view (sums over the nine databases):")
print(
  res |>
    group_by(dataset, region) |>
    summarise(
      best_cover_ge95 = paste0(
        round(min(best_cover_ge95), 3),
        "-",
        round(max(best_cover_ge95), 3)
      ),
      genus_80 = sum(genus_80),
      lost = sum(lost),
      lost_right = sum(lost_right),
      lost_wrong = sum(lost_wrong),
      lost_no_truth = sum(lost_no_truth),
      gained = sum(gained),
      changed = sum(changed),
      .groups = "drop"
    ) |>
    as.data.frame()
)
