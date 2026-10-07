# Timing of the four methods against query count and reference size
# (developer, 2026-09-26). On the wood and air runs lca took ~4 h and blastn
# ~2 h per computation, far beyond what the mock suggested, so this separates
# the fixed cost of a computation (reading or indexing the reference) from its
# cost per query, and shows how both scale with the reference.
#
# Grid: 10, 25, 50 and 100 ASVs of the Pauvert mock (one seeded permutation,
# nested) x Unite_all_20250219_Fungi whole and as random subsets 2, 5, 100 and
# 1000 times smaller (the same records in the dada2 and sintax files) x the 4
# methods. Each cell runs ALONE, through the pipeline's own
# compute_assignment() (R/assign_compute.R), with the pipeline's threads
# (assign_threads_dada2 for dada2, assign_threads_fast otherwise, config.R),
# so a time here is what one computation costs without contention (OD 31).
# Rows are appended to the csv as they are measured and a rerun skips them.
#
# The larger points (mock, in silico, wood: 206 to 1477 units, controls
# included) come from the benchmark_costs of the stores, same database, but
# measured with 2 workers + the dada2 worker running at once.
#
# Run from the project root (about 2-3 h):
#   Rscript analysis/sandbox/timing_query_x_refsize.R

here::i_am("analysis/sandbox/timing_query_x_refsize.R")
source(here::here("config.R"))
source(here::here("R/load_pqverse.R"))
source(here::here("R/assign_compute.R"))
suppressMessages(load_pqverse(c("MiscMetabar", "comparpq")))
suppressMessages(library(phyloseq))

db <- "Unite_all_20250219_Fungi"
n_queries <- c(10, 25, 50, 100)
divisors <- c(1, 2, 5, 100, 1000)
methods <- c("sintax", "lca", "blastn", "dada2")
out_csv <- here::here("data/data_final/timing/timing_query_x_refsize.csv")
sub_dir <- here::here("data/data_intermediate/timing_refs")
dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
dir.create(sub_dir, showWarnings = FALSE, recursive = TRUE)

# Queries: one permutation of the Pauvert ASVs, the first n of it for each n.
asv <- targets::tar_read(d_asv_common, store = here::here("store_assign_taxo"))
set.seed(targets_seed)
order_asv <- sample(taxa_names(asv))

# References: one permutation of the records, the first n / divisor of it,
# written in both formats from the same positions.
ref_path <- function(fmt) {
  here::here(
    "data/data_raw/refseq",
    paste0(fmt, "_format"),
    paste0(db, ".fasta")
  )
}
index <- lapply(
  c(dada2 = "dada2", sintax = "sintax"),
  \(fmt) Biostrings::fasta.index(ref_path(fmt))
)
stopifnot(identical(index$dada2$seqlength, index$sintax$seqlength))
set.seed(targets_seed)
order_ref <- sample(nrow(index$sintax))
subset_ref <- function(divisor, fmt) {
  if (divisor == 1) {
    return(ref_path(fmt))
  }
  path <- file.path(sub_dir, sprintf("%s_div%d_%s.fasta", db, divisor, fmt))
  if (!file.exists(path)) {
    keep <- sort(order_ref[seq_len(ceiling(length(order_ref) / divisor))])
    Biostrings::writeXStringSet(
      Biostrings::readBStringSet(index[[fmt]][keep, ]),
      path,
      width = 20000
    )
  }
  path
}

done <- if (file.exists(out_csv)) {
  utils::read.csv(out_csv)
} else {
  data.frame(method = character(), divisor = numeric(), n_query = numeric())
}

for (divisor in divisors) {
  for (method in methods) {
    fmt <- if (method == "dada2") "dada2" else "sintax"
    ref <- subset_ref(divisor, fmt)
    for (n in n_queries) {
      if (
        any(
          done$method == method & done$divisor == divisor & done$n_query == n
        )
      ) {
        next
      }
      physeq <- prune_taxa(order_asv[seq_len(n)], asv)
      threads <- if (method == "dada2") {
        assign_threads_dada2
      } else {
        assign_threads_fast
      }
      elapsed <- system.time(
        compute_assignment(physeq, method, ref, nproc = threads)
      )[["elapsed"]]
      row <- data.frame(
        method = method,
        db = db,
        divisor = divisor,
        n_ref = ceiling(nrow(index$sintax) / divisor),
        n_query = n,
        threads = threads,
        elapsed_s = elapsed,
        source = "timing run, alone",
        date = format(Sys.time(), "%Y-%m-%d %H:%M")
      )
      utils::write.table(
        row,
        out_csv,
        sep = ",",
        row.names = FALSE,
        col.names = !file.exists(out_csv),
        append = file.exists(out_csv)
      )
      message(sprintf(
        "%-6s 1/%-4d n=%3d  %8.1f s",
        method,
        divisor,
        n,
        elapsed
      ))
    }
  }
}

# The larger points already measured on the same database, per store.
stores <- c(
  mock = "store_assign_taxo",
  insilico = "store_assign_taxo_insilico",
  wood = "store_assign_taxo_korhonen_wood"
)
pipeline_points <- do.call(
  rbind,
  lapply(names(stores), function(s) {
    costs <- targets::tar_read(benchmark_costs, store = here::here(stores[[s]]))
    units <- ntaxa(targets::tar_read(
      d_asv_for_assignation,
      store = here::here(stores[[s]])
    ))
    costs <- unique(costs[
      costs$db == db & costs$preprocess == "none",
      c("method", "wall_time_s")
    ])
    data.frame(
      method = costs$method,
      db = db,
      divisor = 1,
      n_ref = nrow(index$sintax),
      n_query = units,
      threads = ifelse(
        costs$method == "dada2",
        assign_threads_dada2,
        assign_threads_fast
      ),
      elapsed_s = costs$wall_time_s,
      source = paste0("pipeline run (", s, "), concurrent"),
      date = NA
    )
  })
)
utils::write.csv(
  pipeline_points,
  sub("\\.csv$", "_pipeline_points.csv", out_csv),
  row.names = FALSE
)
message("done: ", out_csv)
