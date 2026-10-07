# Cost of a species-assignment layer on the whole grid (developer,
# 2026-10-06): every current input of every dataset x the nine databases, with
# the two exact-match methods asked for.
#
# dada2::assignSpecies(), against data/data_raw/refseq/species_format/
# (make_databases.R::derive_species()): a query is assigned when it is an
# exact substring of a record. Its cost is one full scan of the reference per
# distinct query length, plus reading the reference. So:
#   1. calibration: Pauvert ITS1 ASV against the nine databases, the
#      reading and the per-length scan of each database measured apart;
#   2. validation: the model against measured runs on the inputs with the
#      most distinct lengths (air ITS2, Tedersoo full);
#   3. prediction: every input x database.
#
# vsearch --search_exact, against the sintax files: full-length identity only
# (query and record identical end to end), so an ITSx-cut query can match only
# a record cut to the same region. Measured directly on Pauvert ITS1 ASV and
# Tedersoo full ASV against the nine databases: time and queries matched.
#
# Queries with a non-ACGT base are left out (assignSpecies() refuses them).
# Single thread for both methods. Reads the stores; writes
# data/snapshots/species_timing/. Run from the project root:
#   Rscript analysis/sandbox/timing_species_assignment.R

suppressMessages({
  library(here)
  library(targets)
  library(dplyr)
  library(Biostrings)
  library(phyloseq)
  source(here("config.R"))
  source(here("R/load_pqverse.R"))
  load_pqverse("dbpq")
})

out_dir <- here("data/snapshots/species_timing")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stores <- c(
  pauvert = "store_assign_taxo",
  tedersoo = "store_assign_taxo_tedersoo_illumina",
  insilico = "store_assign_taxo_insilico",
  wood = "store_assign_taxo_korhonen_wood",
  air = "store_assign_taxo_gssp_air"
)
species_file <- function(db) {
  here("data/data_raw/refseq/species_format", paste0(db, ".fasta"))
}
sintax_file <- function(db) {
  here("data/data_raw/refseq/sintax_format", paste0(db, ".fasta"))
}

# ---- the queries of every input ----------------------------------------------
inputs <- bind_rows(lapply(names(stores), function(ds) {
  names_in <- grep(
    "^d_(ITS1|ITS2|full)_(asv|otu)_for_assignation$",
    tar_meta(store = here(stores[[ds]]), fields = "name")$name,
    value = TRUE
  )
  tibble(dataset = ds, target = names_in)
}))
queries <- lapply(seq_len(nrow(inputs)), function(i) {
  s <- tar_read_raw(
    inputs$target[i],
    store = here(stores[[inputs$dataset[i]]])
  )@refseq
  s[!grepl("[^ACGT]", as.character(s))]
})
inputs$input <- sub("^d_(.*)_for_assignation$", "\\1", inputs$target)
inputs$queries <- lengths(queries)
inputs$lengths <- vapply(queries, \(s) length(unique(width(s))), integer(1))
print(select(inputs, dataset, input, queries, lengths), n = Inf)

assign_species_time <- function(seqs, db) {
  system.time(
    out <- dada2::assignSpecies(
      as.character(seqs),
      species_file(db),
      allowMultiple = TRUE,
      tryRC = FALSE
    )
  )[["elapsed"]] |>
    structure(assigned = sum(!is.na(out[, "Species"])))
}

# ---- 1. calibration ------------------------------------------------------------
pauvert <- which(inputs$dataset == "pauvert" & inputs$input == "ITS1_asv")
calib_file <- file.path(out_dir, "assign_species_calibration.csv")
if (file.exists(calib_file)) {
  calib <- read.csv(calib_file)
} else {
  calib <- bind_rows(lapply(db_list, function(db) {
    message(format(Sys.time(), "%H:%M"), " calibration ", db)
    read_s <- system.time(ShortRead::readFasta(species_file(db)))[["elapsed"]]
    total <- assign_species_time(queries[[pauvert]], db)
    tibble(
      db = db,
      read_s = read_s,
      total_s = as.numeric(total),
      scan_s = (as.numeric(total) - read_s) / inputs$lengths[pauvert],
      assigned = attr(total, "assigned")
    )
  }))
  write.csv(calib, calib_file, row.names = FALSE)
}
cat(
  "\nassignSpecies calibration on Pauvert ITS1 ASV (",
  inputs$queries[pauvert],
  " queries, ",
  inputs$lengths[pauvert],
  " distinct lengths)\n",
  sep = ""
)
print(calib)

predict_s <- function(lengths, db) {
  k <- match(db, calib$db)
  calib$read_s[k] + lengths * calib$scan_s[k]
}

# ---- 2. validation ---------------------------------------------------------------
valid_file <- file.path(out_dir, "assign_species_validation.csv")
checks <- tibble::tribble(
  ~dataset   , ~input     , ~db                        ,
  "air"      , "ITS2_asv" , "Unite_all_20250219_Fungi" ,
  "air"      , "ITS2_asv" , "EUK_ITS_v2.1_Fungi"       ,
  "tedersoo" , "full_asv" , "EUK_ITS_v2.1"
)
if (file.exists(valid_file)) {
  valid <- read.csv(valid_file)
} else {
  valid <- bind_rows(lapply(seq_len(nrow(checks)), function(j) {
    i <- which(
      inputs$dataset == checks$dataset[j] & inputs$input == checks$input[j]
    )
    message(
      format(Sys.time(), "%H:%M"),
      " validation ",
      checks$dataset[j],
      " ",
      checks$db[j]
    )
    measured <- assign_species_time(queries[[i]], checks$db[j])
    tibble(
      checks[j, ],
      lengths = inputs$lengths[i],
      predicted_s = predict_s(inputs$lengths[i], checks$db[j]),
      measured_s = as.numeric(measured),
      assigned = attr(measured, "assigned")
    )
  }))
  write.csv(valid, valid_file, row.names = FALSE)
}
cat("\nassignSpecies validation, predicted against measured\n")
print(valid)

# ---- 3. prediction on the whole grid ---------------------------------------------
grid <- tidyr::expand_grid(
  select(inputs, dataset, input, lengths),
  db = db_list
) |>
  mutate(predicted_s = predict_s(lengths, db))
write.csv(
  grid,
  file.path(out_dir, "assign_species_prediction.csv"),
  row.names = FALSE
)
cat(sprintf(
  "\nassignSpecies, %d inputs x %d databases: %.1f h on one core\n",
  nrow(inputs),
  length(db_list),
  sum(grid$predicted_s) / 3600
))
print(
  grid |>
    group_by(dataset) |>
    summarise(hours = round(sum(predicted_s) / 3600, 2)),
  n = Inf
)

# ---- vsearch --search_exact ---------------------------------------------------------
exact_file <- file.path(out_dir, "search_exact.csv")
if (file.exists(exact_file)) {
  exact <- read.csv(exact_file)
} else {
  sample_inputs <- which(
    (inputs$dataset == "pauvert" & inputs$input == "ITS1_asv") |
      (inputs$dataset == "tedersoo" & inputs$input == "full_asv")
  )
  exact <- bind_rows(lapply(sample_inputs, function(i) {
    bind_rows(lapply(db_list, function(db) {
      message(
        format(Sys.time(), "%H:%M"),
        " search_exact ",
        inputs$dataset[i],
        " ",
        db
      )
      q <- tempfile(fileext = ".fasta")
      o <- tempfile(fileext = ".tsv")
      on.exit(unlink(c(q, o)))
      writeXStringSet(queries[[i]], q, width = 20000)
      elapsed <- system.time(system2(
        dbpq::find_vsearch(),
        c(
          "--search_exact",
          q,
          "--db",
          sintax_file(db),
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
      hits <- if (file.exists(o) && file.size(o) > 0) {
        read.delim(o, header = FALSE)
      } else {
        data.frame(V1 = character())
      }
      tibble(
        dataset = inputs$dataset[i],
        input = inputs$input[i],
        db = db,
        queries = inputs$queries[i],
        matched = length(unique(hits$V1)),
        seconds = elapsed
      )
    }))
  }))
  write.csv(exact, exact_file, row.names = FALSE)
}
cat("\nvsearch --search_exact (full-length identity), single thread\n")
print(exact, n = Inf)
