# stampa LCA (current) against the common LCA of QIIME 2 / Bokulich et al. 2018
# (developer, 2026-09-26). The pipeline runs assign_vsearch_lca() at its
# MiscMetabar defaults, the stampa algorithm (top_hits_only = TRUE, id = 0.5,
# maxaccepts = 0), which timing_by_query_class.R found near-exhaustive on real
# ITS queries (~4.7 s per wood OTU against 102 137 records). The common
# configuration is the "balanced" one Bokulich et al. 2018 optimised on 16S
# and fungal ITS and set as the defaults of QIIME 2's
# classify-consensus-vsearch: top_hits_only = FALSE, id (perc_identity) = 0.8,
# maxaccepts = 10, consensus (lca_cutoff) = 0.51; it is also run at vsearch's
# own strict default lca_cutoff = 1. maxrejects stays at vsearch's 32.
#
# Same queries as timing_by_query_class.R (up to 100 per class and dataset),
# against Unite_all_20250219_Fungi, one computation at a time. Per
# configuration: time per query, share of queries given a genus, and for the
# real units of the mock and the in silico set, the share given the true
# genus (submitted or GBIF accepted name) or a wrong one.
#
# Run from the project root (a few minutes, plus ~20 min for stampa):
#   Rscript analysis/sandbox/compare_lca_stampa_vs_common.R

here::i_am("analysis/sandbox/compare_lca_stampa_vs_common.R")
source(here::here("config.R"))
source(here::here("R/load_pqverse.R"))
source(here::here("R/mock_truth.R"))
suppressMessages(load_pqverse(c("MiscMetabar", "comparpq")))
suppressMessages(library(phyloseq))

db <- "Unite_all_20250219_Fungi"
ref <- here::here("data/data_raw/refseq/sintax_format", paste0(db, ".fasta"))
out_csv <- here::here("data/data_final/timing/compare_lca_stampa_vs_common.csv")
stores <- c(
  mock = "store_assign_taxo",
  insilico = "store_assign_taxo_insilico",
  wood = "store_assign_taxo_korhonen_wood"
)
configs <- list(
  stampa = list(top_hits_only = TRUE, id = 0.5, maxaccepts = 0, lca_cutoff = 1),
  common_0.51 = list(
    top_hits_only = FALSE,
    id = 0.8,
    maxaccepts = 10,
    lca_cutoff = 0.51
  ),
  common_1 = list(
    top_hits_only = FALSE,
    id = 0.8,
    maxaccepts = 10,
    lca_cutoff = 1
  )
)
not_a_genus <- "ncertae_sedis|^unclassified|^$"

# Genus truth of the real units, NULL when the dataset has none.
truth_mock <- build_mock_truth(
  targets::tar_read(d_asv, store = here::here("store_dada2")),
  mock_sanger_tables[mock_sanger_tables$dataset == "pauvert", ]
)
input_seqs <- Biostrings::readDNAStringSet(
  here::here("data/data_raw/insilico/iss_input_fungi.fasta")
)
records <- data.frame(
  strain = sprintf("rec_%03d", seq_along(input_seqs)),
  name = names(input_seqs),
  sequence = as.character(input_seqs)
)
tax_str <- sub(".*\\|", "", records$name)
lineage <- data.frame(strain = records$strain)
for (i in seq_along(truth_ranks)) {
  prefix <- c("k__", "p__", "c__", "o__", "f__", "g__", "s__")[i]
  lineage[[truth_ranks[i]]] <- sub(
    paste0("^.*", prefix, "([^;]*).*$"),
    "\\1",
    tax_str
  )
}
lineage$Species <- species_name(lineage$Species)
lineage$Genus_accepted <- NA_character_
lineage$Species_accepted <- NA_character_
truth_insilico <- mock_truth_table(
  targets::tar_read(d_asv, store = here::here("store_dada2_insilico"))@refseq,
  records,
  lineage
)
truths <- list(mock = truth_mock, insilico = truth_insilico, wood = NULL)

rows <- list()
for (s in names(stores)) {
  pq <- targets::tar_read(
    d_asv_for_assignation,
    store = here::here(stores[[s]])
  )
  cls <- ifelse(
    grepl("^fake_", taxa_names(pq)),
    "fake",
    ifelse(grepl("^external_", taxa_names(pq)), "external", "real")
  )
  for (k in c("real", "fake", "external")) {
    set.seed(targets_seed)
    pick <- sample(taxa_names(pq)[cls == k], min(100, sum(cls == k)))
    sub <- prune_taxa(pick, pq)
    for (cfg in names(configs)) {
      args <- c(
        list(
          sub,
          ref_fasta = ref,
          behavior = "return_matrix",
          clean_pq = FALSE,
          nproc = assign_threads_fast,
          verbose = FALSE
        ),
        configs[[cfg]]
      )
      elapsed <- system.time(
        res <- do.call(MiscMetabar::assign_vsearch_lca, args)
      )[["elapsed"]]
      genus <- res$Genus_sintax[match(pick, res$taxa_names)]
      genus[!is.na(genus) & grepl(not_a_genus, genus)] <- NA
      correct <- wrong <- NA_real_
      truth <- truths[[s]]
      if (k == "real" && !is.null(truth)) {
        tr <- truth[match(pick, truth$taxon), ]
        scored <- !is.na(tr$Genus) &
          tr$truth_depth %in% c("Genus", "Species")
        hit <- !is.na(genus) &
          (genus == tr$Genus |
            (!is.na(tr$Genus_accepted) & genus == tr$Genus_accepted))
        correct <- mean(hit[scored])
        wrong <- mean((!is.na(genus) & !hit)[scored])
      }
      rows[[length(rows) + 1]] <- data.frame(
        dataset = s,
        class = k,
        config = cfg,
        n_query = length(pick),
        s_per_query = elapsed / length(pick),
        genus_assigned = mean(!is.na(genus)),
        genus_correct = correct,
        genus_wrong = wrong
      )
      message(sprintf(
        "%-8s %-8s %-11s %6.2f s/q  genus %4.2f  correct %4.2f  wrong %4.2f",
        s,
        k,
        cfg,
        elapsed / length(pick),
        mean(!is.na(genus)),
        correct,
        wrong
      ))
      utils::write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
    }
  }
}
