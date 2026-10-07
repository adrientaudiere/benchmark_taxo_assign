# Is each method still reproducible when a computation uses every core
# (developer, 2026-09-27)? The final timing run is to run one computation at a
# time but give it all the threads of the machine (OD 31, revised). vsearch
# documents a seeded sintax as reproducible only with one thread (ROADMAP
# SA6), so each method is run twice with nproc = all cores, through the
# pipeline's compute_assignment(), and the two outputs compared. sintax runs
# against EUK_ITS_v2.1_Fungi, where the unseeded runs differed; the others
# against Unite_all_20250219_Fungi. Queries: the 206 units of the in silico
# assignment input.
#   Rscript analysis/sandbox/check_threads_determinism.R

here::i_am("analysis/sandbox/check_threads_determinism.R")
source(here::here("config.R"))
source(here::here("R/load_pqverse.R"))
source(here::here("R/assign_compute.R"))
suppressMessages(load_pqverse(c("MiscMetabar", "comparpq")))
suppressMessages(library(phyloseq))

pq <- targets::tar_read(
  d_asv_for_assignation,
  store = here::here("store_assign_taxo_insilico")
)
threads <- parallel::detectCores()
ref_of <- function(method) {
  db <- if (method == "sintax") {
    "EUK_ITS_v2.1_Fungi"
  } else {
    "Unite_all_20250219_Fungi"
  }
  fmt <- if (method == "dada2") "dada2" else "sintax"
  here::here(
    "data/data_raw/refseq",
    paste0(fmt, "_format"),
    paste0(db, ".fasta")
  )
}
for (method in c("sintax", "lca", "blastn", "dada2")) {
  runs <- lapply(1:2, function(i) {
    t <- system.time(
      out <- compute_assignment(pq, method, ref_of(method), nproc = threads)
    )[["elapsed"]]
    list(out = out, t = t)
  })
  a <- runs[[1]]$out
  b <- runs[[2]]$out
  if (method == "blastn") {
    key <- function(x) x[order(x[[1]], x[[3]], x[[6]]), ]
    a <- key(a)
    b <- key(b)
    rownames(a) <- rownames(b) <- NULL
  }
  if (method == "lca") {
    a <- a[order(a$query, a$target), ]
    b <- b[order(b$query, b$target), ]
    rownames(a) <- rownames(b) <- NULL
  }
  message(sprintf(
    "%-6s %d threads: identical %s | %.1f s and %.1f s",
    method,
    threads,
    identical(a, b),
    runs[[1]]$t,
    runs[[2]]$t
  ))
}
