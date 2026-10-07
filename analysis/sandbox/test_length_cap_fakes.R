suppressMessages({
  library(targets)
  library(phyloseq)
  library(dplyr)
  library(here)
})
source("config.R")
source("R/load_pqverse.R")
source("R/assign_compute.R")
source("make_databases.R")
load_pqverse(c("MiscMetabar", "dbpq"))

# Dedicated test of OD 35 for the article (developer, 2026-10-01): the
# shuffled `fake_` controls only, assigned against the whole EUK_ITS_v2.1
# release with and without its records over 5 000 bp, by the four methods at
# their default setting. Self-contained: both references are rebuilt from the
# source release, the queries are read from the stores once and saved.
#
#   queries: the `fake_` controls of Pauvert ITS1 ASV (97) and of Tedersoo
#            full-ITS ASV (118), seeded shuffles of the real units (their only
#            correct answer is NA)
#   references: uncapped = dada2 / sintax files rebuilt from
#            data/data_raw/refseq/sources/EUK_ITS_v2.1.fasta into
#            data/data_raw/refseq/test_uncapped/; capped = the project files
#            (make_databases.R step 2a-bis)
#   methods: compute_assignment() then derive_assignment() of the pipeline
#            (R/assign_compute.R), same threads, sintax seed and parameters
#
# One raw result per method x reference in data/snapshots/length_cap_test/,
# so the methods can run at different times (dada2 needs ~54 GB of RAM per
# reference: run it with the desktop closed). The summary is printed from
# whatever results exist.
#   Rscript analysis/sandbox/test_length_cap_fakes.R sintax lca blastn
#   Rscript analysis/sandbox/test_length_cap_fakes.R dada2
#   Rscript analysis/sandbox/test_length_cap_fakes.R          (summary only)

methods <- commandArgs(trailingOnly = TRUE)
stopifnot(all(methods %in% c("dada2", "sintax", "lca", "blastn")))
out_dir <- "data/snapshots/length_cap_test"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Queries ---------------------------------------------------------------------

query_file <- file.path(out_dir, "fake_queries.rds")
if (!file.exists(query_file)) {
  read_fakes <- function(target, store, dataset) {
    d <- tar_read_raw(target, store = store)
    fakes <- d@refseq[grepl("^fake_", names(d@refseq))]
    names(fakes) <- paste0(dataset, "__", names(fakes))
    fakes
  }
  seqs <- c(
    read_fakes("d_ITS1_asv_for_assignation", "store_assign_taxo", "pauvert"),
    read_fakes(
      "d_full_asv_for_assignation",
      "store_assign_taxo_tedersoo_illumina",
      "tedersoo"
    )
  )
  saveRDS(seqs, query_file)
}
seqs <- readRDS(query_file)
# A one-sample phyloseq holding the fakes, as compute_assignment() expects.
queries <- phyloseq(
  otu_table(
    matrix(1, nrow = 1, ncol = length(seqs), dimnames = list("s", names(seqs))),
    taxa_are_rows = FALSE
  ),
  # The derivations add their columns to an existing tax_table.
  tax_table(matrix(
    sub("__.*$", "", names(seqs)),
    ncol = 1,
    dimnames = list(names(seqs), "Dataset")
  )),
  refseq(seqs)
)
message(
  "queries: ",
  sum(startsWith(names(seqs), "pauvert__")),
  " Pauvert ITS1 ASV fakes, ",
  sum(startsWith(names(seqs), "tedersoo__")),
  " Tedersoo full-ITS ASV fakes"
)

# References ------------------------------------------------------------------

source_fasta <- refseq_path("EUK_ITS_v2.1", "general")
uncapped_dir <- "data/data_raw/refseq/test_uncapped"
reference <- function(db, method) {
  fmt <- if (method == "dada2") "dada2" else "sintax"
  if (db == "capped") {
    path <- refseq_path("EUK_ITS_v2.1", fmt)
    stopifnot(max(Biostrings::fasta.seqlengths(path)) <= 5000)
    return(path)
  }
  path <- file.path(uncapped_dir, fmt, "EUK_ITS_v2.1.fasta")
  if (!file.exists(path)) {
    derive <- if (fmt == "dada2") derive_dada2 else derive_sintax
    derive(source_fasta, path)
  }
  path
}

# Runs ------------------------------------------------------------------------

for (method in methods) {
  for (db in c("uncapped", "capped")) {
    ref <- reference(db, method)
    t0 <- Sys.time()
    res <- compute_assignment(
      queries,
      method,
      ref,
      nproc = if (method == "dada2") {
        assign_threads_dada2
      } else {
        assign_threads_fast
      }
    )
    minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    saveRDS(
      list(
        result = res,
        minutes = minutes,
        records = length(Biostrings::fasta.seqlengths(ref))
      ),
      file.path(out_dir, paste0(method, "__", db, ".rds"))
    )
    message(method, " x ", db, ": ", round(minutes, 1), " min")
  }
}

# Summary ---------------------------------------------------------------------

defaults <- list(
  dada2 = list(min_bootstrap = 0.5),
  sintax = list(min_bootstrap = 0.5),
  lca = list(lca_cutoff = 0.51),
  blastn = list(
    vote_algorithm = "rel_majority",
    nb_voting = 100,
    min_id = 95,
    min_cover = blastn_min_cover
  )
)
placeholder <- "[Ii]ncertae_sedis|^unclassified|^$"

# Records removed by the cap, from the uncapped source release.
src <- Biostrings::fasta.index(source_fasta)
removed <- src[src$seqlength > max_record_length[["EUK_ITS_v2.1"]], ]
removed_id <- sub(";.*$", "", removed$desc)
removed_tax <- vapply(
  strsplit(sub("^[^;]*;", "", removed$desc), ";", fixed = TRUE),
  \(r) paste(sub("^[a-z]__", "", r), collapse = ";"),
  character(1)
)

# Does the answer to each query rest on a removed record?
on_removed <- function(method, computed, units) {
  switch(
    method,
    dada2 = apply(computed$tax, 1, paste, collapse = ";") %in% removed_tax,
    sintax = {
      v <- computed$taxo_value
      tax <- apply(
        v[, setdiff(names(v), "taxa_names")],
        1,
        paste,
        collapse = ";"
      )
      tax[match(units, v$taxa_names)] %in% removed_tax
    },
    lca = units %in%
      computed$query[sub(";.*$", "", computed$target) %in% removed_id],
    blastn = units %in%
      computed$`Query name`[
        sub(";.*$", "", computed$`Taxa name`) %in% removed_id
      ]
  )
}

rows <- list()
for (f in list.files(
  out_dir,
  pattern = "^(dada2|sintax|lca|blastn)__(uncapped|capped)\\.rds$"
)) {
  method <- sub("__.*$", "", f)
  db <- sub("^.*__(.*)\\.rds$", "\\1", f)
  run <- readRDS(file.path(out_dir, f))
  derived <- do.call(
    derive_assignment,
    c(
      list(
        physeq = queries,
        computed = run$result,
        method = method,
        suffix = "_test"
      ),
      defaults[[method]]
    )
  )
  # A computation with no hit at all (blastn on random sequences) returns
  # NULL, and derive_assignment() then adds no column: no fake is named.
  if (is.null(run$result)) {
    named <- rep(FALSE, ntaxa(queries))
    removed_hit <- rep(FALSE, ntaxa(queries))
  } else {
    genus <- as.character(derived@tax_table[, "Genus_test"])
    named <- !is.na(genus) & !grepl(placeholder, genus)
    removed_hit <- on_removed(method, run$result, taxa_names(queries))
  }
  for (ds in c("pauvert", "tedersoo")) {
    q <- startsWith(taxa_names(queries), paste0(ds, "__"))
    rows[[length(rows) + 1]] <- tibble(
      queries = ds,
      method = method,
      reference = db,
      fakes = sum(q),
      named_at_genus = sum(named & q),
      share = round(mean(named[q]), 3),
      named_via_removed_record = sum(named & q & removed_hit),
      minutes = round(run$minutes, 1),
      records = run$records
    )
  }
}
res <- bind_rows(rows) |>
  arrange(queries, method, desc(reference))
print(as.data.frame(res), row.names = FALSE)
write.csv(res, file.path(out_dir, "summary.csv"), row.names = FALSE)
