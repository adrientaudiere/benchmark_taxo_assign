# Which vsearch command finds every reference record that holds a query at
# 100 % (species assignment by exact match, developer 2026-10-06)?
#
# The queries are ITSx-cut units (ITS1, ITS2) while the records are whole ITS
# regions, so an exact match is the query contained in a longer record.
# Ground truth: every record holding the query as an exact substring
# (Biostrings::vcountPDict(), as dada2::assignSpecies() does). Compared on
# Pauvert ITS1 ASV against two databases, single thread:
#   search_exact     vsearch --search_exact (full-length identity only)
#   global_default   --usearch_global --id 1 --query_cov 1, default heuristics
#                    (--maxaccepts 1 --maxrejects 32)
#   global_all       --usearch_global --id 1 --query_cov 1 --maxaccepts 0
#                    --maxrejects 0 (exhaustive: every hit)
#   global_top       --usearch_global --id 1 --query_cov 1 --maxaccepts 0
#                    --maxrejects 32 (all hits among the k-mer candidates)
# For each: time, queries matched, and the share of the true (query, record)
# pairs found. Reads store_assign_taxo; writes
# data/snapshots/species_timing/vsearch_exact_check.csv.
#   Rscript analysis/sandbox/check_vsearch_exact_species.R

suppressMessages({
  library(here)
  library(targets)
  library(dplyr)
  library(Biostrings)
  library(phyloseq)
  source(here("R/load_pqverse.R"))
  load_pqverse("dbpq")
})
out_dir <- here("data/snapshots/species_timing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dbs <- c("Unite_all_20250219_Fungi", "EUK_ITS_v2.1_Fungi")

q <- tar_read(
  d_ITS1_asv_for_assignation,
  store = here("store_assign_taxo")
)@refseq
q <- q[!grepl("[^ACGT]", as.character(q))]
q_file <- tempfile(fileext = ".fasta")
writeXStringSet(q, q_file, width = 20000)

truth_pairs <- function(db_file) {
  ref <- readDNAStringSet(db_file)
  ref_id <- sub(";.*$", "", names(ref))
  bind_rows(lapply(unique(width(q)), function(len) {
    sub_q <- q[width(q) == len]
    hit <- vcountPDict(PDict(sub_q), ref) > 0
    idx <- which(hit, arr.ind = TRUE)
    tibble(query = names(sub_q)[idx[, 1]], target = ref_id[idx[, 2]])
  }))
}

vsearch_pairs <- function(db_file, args) {
  o <- tempfile(fileext = ".tsv")
  on.exit(unlink(o))
  elapsed <- system.time(system2(
    dbpq::find_vsearch(),
    c(
      args[1],
      q_file,
      "--db",
      db_file,
      args[-1],
      "--strand",
      "plus",
      "--threads",
      "1",
      "--userout",
      o,
      "--userfields",
      "query+target",
      "--quiet"
    )
  ))[["elapsed"]]
  hits <- if (file.size(o) > 0) {
    read.delim(o, header = FALSE, col.names = c("query", "target"))
  } else {
    tibble(query = character(), target = character())
  }
  hits$target <- sub(";.*$", "", hits$target)
  structure(hits, elapsed = elapsed)
}

global <- c("--usearch_global", "--id", "1", "--query_cov", "1")
modes <- list(
  search_exact = "--search_exact",
  global_default = global,
  global_all = c(global, "--maxaccepts", "0", "--maxrejects", "0"),
  global_top = c(global, "--maxaccepts", "0", "--maxrejects", "32")
)

res <- bind_rows(lapply(dbs, function(db) {
  db_file <- here("data/data_raw/refseq/sintax_format", paste0(db, ".fasta"))
  t_truth <- system.time(truth <- truth_pairs(db_file))[["elapsed"]]
  key_truth <- paste(truth$query, truth$target)
  bind_rows(
    tibble(
      db = db,
      mode = "truth (substring)",
      seconds = t_truth,
      queries_matched = n_distinct(truth$query),
      pairs = nrow(truth),
      pairs_found = 1
    ),
    bind_rows(lapply(names(modes), function(m) {
      message(format(Sys.time(), "%H:%M"), " ", db, " ", m)
      hits <- vsearch_pairs(db_file, modes[[m]])
      tibble(
        db = db,
        mode = m,
        seconds = attr(hits, "elapsed"),
        queries_matched = n_distinct(hits$query),
        pairs = nrow(hits),
        pairs_found = mean(key_truth %in% paste(hits$query, hits$target))
      )
    }))
  )
}))
write.csv(res, file.path(out_dir, "vsearch_exact_check.csv"), row.names = FALSE)
cat("\nPauvert ITS1 ASV,", length(q), "queries\n")
print(
  mutate(res, seconds = round(seconds, 1), pairs_found = round(pairs_found, 3)),
  n = Inf
)
