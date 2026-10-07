# OD 39 (developer, 2026-10-06): redraw the external controls with the fixed
# derive_fake_ref() (code review §1.2), then the `_Fungi_rep` databases, whose
# representatives must stay away from the controls, then their species_format
# files (OD 43). The previous draw is kept in
# data/snapshots/controls_before_OD39_2026-10-06/.
#
# Every assignment store tracks fake_ref_asv_100.fasta and the reference
# files, so this outdates every input and computation: the grid reruns.
#   Rscript analysis/sandbox/redraw_controls_od39.R

suppressMessages(source("make_databases.R"))
before_dir <- here("data/snapshots/controls_before_OD39_2026-10-06")

phyla_of <- function(path) {
  h <- names(Biostrings::readDNAStringSet(path))
  stringr::str_match(h, "p__\\s*(.*?)\\s*;c__")[, 2]
}

# ---- 1. controls -------------------------------------------------------------
pool <- file.path(sources_dir, paste0(fake_ref_source, "_wo_Fungi.fasta"))
kingdoms <- external_control_kingdoms()
derive_fake_ref(
  pool,
  here(fake_ref_fasta),
  seed = targets_seed,
  kingdoms = kingdoms,
  force = TRUE
)
pool_seqs <- Biostrings::readDNAStringSet(pool)
pool_kingdom <- stringr::str_match(names(pool_seqs), "k__([^;]+);")[, 2]
pool_phyla <- unique(stringr::str_match(
  names(pool_seqs)[pool_kingdom %in% kingdoms],
  "p__\\s*(.*?)\\s*;c__"
)[, 2])
new_phyla <- phyla_of(here(fake_ref_fasta))
old_phyla <- phyla_of(file.path(before_dir, "fake_ref_asv_100.fasta"))
new_names <- names(Biostrings::readDNAStringSet(here(fake_ref_fasta)))
old_names <- names(Biostrings::readDNAStringSet(
  file.path(before_dir, "fake_ref_asv_100.fasta")
))
message(sprintf(
  "controls: %d records, %d of %d phyla (before: %d records, %d phyla); %d records kept from the previous draw",
  length(new_names),
  length(unique(new_phyla)),
  length(pool_phyla),
  length(old_names),
  length(unique(old_phyla)),
  length(intersect(new_names, old_names))
))
stopifnot(
  length(new_names) == 100,
  setequal(unique(new_phyla), pool_phyla)
)

# ---- 2. fungi + rep, then their species files ---------------------------------
for (s in unique(benchmark_dbs$source[
  benchmark_dbs$simplification == "Fungi+rep"
])) {
  message(format(Sys.time(), "%H:%M"), " fungi_rep ", s)
  derive_fungi_rep(s, force = TRUE)
  db <- paste0(s, db_suffix[["Fungi+rep"]])
  old <- read.csv(file.path(before_dir, paste0(db, "_reps.csv")))
  new <- read.csv(file.path(sources_dir, paste0(db, "_reps.csv")))
  message(sprintf(
    "  %s: %d representatives (before %d), %d in common",
    db,
    nrow(new),
    nrow(old),
    length(intersect(new$record, old$record))
  ))
  derive_species(
    refseq_path(db, "sintax"),
    refseq_path(db, "species"),
    force = TRUE
  )
}
message(format(Sys.time(), "%H:%M"), " done")
