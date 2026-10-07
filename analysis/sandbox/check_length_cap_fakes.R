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
load_pqverse("MiscMetabar")

# Dedicated test of OD 35 (developer, 2026-10-01): the 5 000 bp record cap
# applied to the whole EUK_ITS_v2.1 release, read on the `fake_` controls
# only. The fakes are seeded shuffles of the real units, so their only correct
# answer is NA, and the same fakes are queried before and after the cap: each
# control is its own pair.
#   uncapped: data/snapshots/uncapped_EUK_ITS_v2.1_2026-10-01/ (the raw
#             computations of the 2026-09-28 to 2026-10-01 runs, dada2 patched)
#   capped:   the same computation targets in the stores, once rebuilt on the
#             capped file (built after it, else reported as pending)
# Default settings (dada2 / sintax bootstrap 0.5, lca 0.51, blastn
# rel_majority 95 %), a genus counted only when clean_taxo() would keep it.
# For the uncapped fakes given a genus, the test also says whether the answer
# rests on one of the 7 871 removed records (> 5 000 bp): the dada2 winner's
# taxonomy, the sintax taxonomy, or any lca / blastn hit.
#   Rscript analysis/sandbox/check_length_cap_fakes.R

snap_dir <- "data/snapshots/uncapped_EUK_ITS_v2.1_2026-10-01"
index <- read.csv(file.path(snap_dir, "index.csv"))
stores <- c(
  pauvert = "store_assign_taxo",
  korhonen_wood = "store_assign_taxo_korhonen_wood",
  tedersoo_illumina = "store_assign_taxo_tedersoo_illumina"
)
capped_file <- refseq_path("EUK_ITS_v2.1", "dada2")
cap_time <- file.mtime(capped_file)
stopifnot(max(Biostrings::fasta.seqlengths(capped_file)) <= 5000)

# The removed records, read from the uncapped source release.
src <- Biostrings::fasta.index(refseq_path("EUK_ITS_v2.1", "general"))
removed <- src[src$seqlength > max_record_length[["EUK_ITS_v2.1"]], ]
removed_id <- sub(";.*$", "", removed$desc)
removed_ranks <- strsplit(sub("^[^;]*;", "", removed$desc), ";", fixed = TRUE)
removed_ranks <- vapply(
  removed_ranks,
  \(r) paste(sub("^[a-z]__", "", r), collapse = ";"),
  character(1)
)
message(length(removed_id), " removed records")

placeholder <- "[Ii]ncertae_sedis|^unclassified|^$"
genus_named <- function(physeq_new, suffix) {
  col <- paste0("Genus", suffix)
  v <- as.character(physeq_new@tax_table[, col])
  !is.na(v) & !grepl(placeholder, v)
}

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

# Does the answer to each unit rest on a removed record?
rests_on_removed <- function(method, computed, units) {
  switch(
    method,
    dada2 = {
      tax <- apply(computed$tax, 1, \(r) paste(r, collapse = ";"))
      tax %in% removed_ranks
    },
    sintax = {
      v <- computed$taxo_value
      ranks <- setdiff(names(v), "taxa_names")
      tax <- apply(v[, ranks], 1, \(r) paste(r, collapse = ";"))
      tax[match(units, v$taxa_names)] %in% removed_ranks
    },
    lca = {
      hit <- sub(";.*$", "", computed$target) %in% removed_id
      units %in% computed$query[hit]
    },
    blastn = {
      hit <- sub(";.*$", "", computed$`Taxa name`) %in% removed_id
      units %in%
        sub("\\|.*$", "", computed$`Query name`[hit]) |
        units %in% computed$`Query name`[hit]
    }
  )
}

rows <- lapply(seq_len(nrow(index)), function(i) {
  ds <- index$dataset[i]
  target <- index$target[i]
  input <- sub("_compute_.*$", "", target)
  method <- sub("^.*_compute_([a-z0-9]+)__.*$", "\\1", target)
  physeq <- tar_read_raw(
    paste0("d_", input, "_for_assignation"),
    store = stores[[ds]]
  )
  units <- taxa_names(physeq)
  fake <- grepl("^fake_", units)
  derive <- function(computed) {
    do.call(
      derive_assignment,
      c(
        list(
          physeq = physeq,
          computed = computed,
          method = method,
          suffix = "_cap_test"
        ),
        defaults[[method]]
      )
    )
  }
  uncapped <- readRDS(file.path(snap_dir, index$file[i]))
  named_un <- genus_named(derive(uncapped), "_cap_test")
  on_removed <- rests_on_removed(method, uncapped, units)

  meta <- tar_meta(names = all_of(target), store = stores[[ds]], fields = time)
  rebuilt <- nrow(meta) == 1 && !is.na(meta$time) && meta$time > cap_time
  named_cap <- if (rebuilt) {
    genus_named(derive(tar_read_raw(target, store = stores[[ds]])), "_cap_test")
  } else {
    rep(NA, length(units))
  }
  tibble(
    dataset = ds,
    input = input,
    method = method,
    fakes = sum(fake),
    named_uncapped = sum(named_un & fake),
    of_which_on_removed = sum(named_un & fake & on_removed),
    named_capped = if (rebuilt) sum(named_cap & fake) else NA_integer_,
    status = if (rebuilt) "paired" else "capped run pending"
  )
})
res <- bind_rows(rows) |> arrange(dataset, input, method)
print(as.data.frame(res), row.names = FALSE)
