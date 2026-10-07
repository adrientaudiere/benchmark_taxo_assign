# Why is lca 50 times slower in the pipeline than in timing_query_x_refsize.R
# (developer, 2026-09-26)? On Unite_all_20250219_Fungi, lca took 80 s for the
# mock's 392 units but 366 s for the in silico set's 206 and 77 min for wood's
# 1412, while 100 real mock ASVs alone take 27 s. The per-query cost must
# depend on what the queries are. This times lca and blastn alone, on up to 100
# queries of each class (real units, fake_ shuffled controls, external_
# non-fungal controls) of each assigned dataset, with the pipeline's own
# compute_assignment() and threads, and records the query lengths.
#
# Run from the project root (a few minutes to an hour, depending on the answer):
#   Rscript analysis/sandbox/timing_by_query_class.R

here::i_am("analysis/sandbox/timing_by_query_class.R")
source(here::here("config.R"))
source(here::here("R/load_pqverse.R"))
source(here::here("R/assign_compute.R"))
suppressMessages(load_pqverse(c("MiscMetabar", "comparpq")))
suppressMessages(library(phyloseq))

db <- "Unite_all_20250219_Fungi"
ref <- here::here("data/data_raw/refseq/sintax_format", paste0(db, ".fasta"))
out_csv <- here::here("data/data_final/timing/timing_by_query_class.csv")
stores <- c(
  mock = "store_assign_taxo",
  insilico = "store_assign_taxo_insilico",
  wood = "store_assign_taxo_korhonen_wood"
)

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
    for (method in c("lca", "blastn")) {
      elapsed <- system.time(
        compute_assignment(sub, method, ref, nproc = assign_threads_fast)
      )[["elapsed"]]
      rows[[length(rows) + 1]] <- data.frame(
        dataset = s,
        class = k,
        method = method,
        n_query = length(pick),
        median_len = median(Biostrings::width(sub@refseq)),
        elapsed_s = elapsed,
        s_per_query = elapsed / length(pick)
      )
      message(sprintf(
        "%-8s %-8s %-6s n=%3d len %4.0f  %7.1f s  (%.2f s/query)",
        s,
        k,
        method,
        length(pick),
        median(Biostrings::width(sub@refseq)),
        elapsed,
        elapsed / length(pick)
      ))
      utils::write.csv(do.call(rbind, rows), out_csv, row.names = FALSE)
    }
  }
}
