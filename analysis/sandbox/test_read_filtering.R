# Read filtering of Rolling et al. 2022 on Pauvert (ROADMAP §1, code review
# §11 Q4, Q9, Q10; developer 2026-10-06). Decides, before M2.4, whether the
# main grid denoises with the settings recommended for fungal ITS: maximum
# expected errors = 8 and reads truncated at the first quality score <= 8.
# No decision rule was fixed in advance: the script reports, the developer
# decides.
#
# Two arms, everything else equal:
#   current: filterAndTrim(trimLeft = 1, trimRight = 1), defaults otherwise
#            (maxEE = Inf, truncQ = 2), as pipelines/dada2.R;
#   rolling: the same plus maxEE = 8, truncQ = 8 (trimLeft / trimRight kept).
# Each arm runs the rest of pipelines/dada2.R from the cutadapt output stored in
# store_dada2 (derep, learnErrors, dada, mergePairs, vsearch chimera removal,
# >= 200 bp, the seed taxonomy), then ITSx ITS1 and the controls as
# pipelines/assign_taxo.R (the stored ITSx external controls, the same shuffle
# seed), then the four methods on the two zoom databases (EUK_ITS_v2.1_Fungi_rep,
# Unite_all_20250219) with the pipeline's threads and seeds.
#
# The current arm is a replica: it must reproduce the stored d_asv and the
# stored d_ITS1_asv_for_assignation exactly, and is then scored from the
# stored computations (test_seed_spread.R checks that path against
# res_comp_tax). If the replica fails, the script stops before any
# comparison.
#
# The filter changes the units, not only their calls, so the report gives the
# reads and units of each arm, the scores on all units (as the pipeline would
# score them) and on the real units common to both arms (matched by ITS1
# sequence, shuffled controls left out), and the seed range of
# test_seed_spread.R beside each difference.
#
# Writes data/data_intermediate/read_filtering_test/ (filtered reads) and
# data/snapshots/read_filtering/ (objects, computations, summaries); reads the
# stores, never writes to them, and logs no autometric file. Run from the
# project root, after test_seed_spread.R (no two dada2 at once):
#   Rscript analysis/sandbox/test_read_filtering.R

suppressMessages({
  source(here::here("analysis/00_setup.R"))
  source(here("analysis/sandbox/seed_filter_helpers.R"))
  source(here("R/itsx_inputs.R"))
})
assert_dada2_patched()

out_dir <- here("data/snapshots/read_filtering")
reads_dir <- here("data/data_intermediate/read_filtering_test")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dstore <- here("store_dada2")
store <- here("store_assign_taxo")
config_sintax_seed <- sintax_seed
methods <- c("dada2", "sintax", "lca", "blastn")

arms <- tibble::tribble(
  ~arm      , ~max_ee , ~trunc_q ,
  "current" , Inf     ,        2 ,
  "rolling" ,       8 ,        8
)

# Saved once in out_dir, read back on a rerun.
cached_rds <- function(name, expr) {
  path <- file.path(out_dir, paste0(name, ".rds"))
  if (file.exists(path)) {
    return(readRDS(path))
  }
  value <- expr
  saveRDS(value, path)
  value
}

# ---- denoising: pipelines/dada2.R from the stored cutadapt output ------------
denoise <- function(arm, max_ee, trunc_q) {
  message(format(Sys.time(), "%H:%M"), " denoising, arm ", arm)
  fnfs <- tar_read(data_fnfs, store = dstore)
  fnrs <- tar_read(data_fnrs, store = dstore)
  filtered <- filter_trim(
    output_fw = file.path(reads_dir, arm, "fwd"),
    output_rev = file.path(reads_dir, arm, "rev"),
    fw = fnfs,
    rev = fnrs,
    multithread = n_threads,
    compress = TRUE,
    trimLeft = 1,
    trimRight = 1,
    maxEE = max_ee,
    truncQ = trunc_q
  )
  derep_fs <- dada2::derepFastq(filtered[[1]])
  derep_rs <- dada2::derepFastq(filtered[[2]])
  set.seed(target_seed("err_fs"))
  err_fs <- dada2::learnErrors(derep_fs, multithread = n_threads)
  set.seed(target_seed("err_rs"))
  err_rs <- dada2::learnErrors(derep_rs, multithread = n_threads)
  set.seed(target_seed("ddF"))
  dd_f <- dada2::dada(derep_fs, err_fs, multithread = n_threads)
  set.seed(target_seed("ddR"))
  dd_r <- dada2::dada(derep_rs, err_rs, multithread = n_threads)
  merged <- dada2::mergePairs(
    dadaF = dd_f,
    dadaR = dd_r,
    derepF = derep_fs,
    derepR = derep_rs
  )
  seq_tab <- dada2::makeSequenceTable(merged)
  wo_chimera <- chimera_removal_vs(seq_tab)
  seqtab <- wo_chimera[, nchar(colnames(wo_chimera)) >= seq_len_min]

  s_d <- tar_read(s_d, store = dstore)
  sam_tab <- tar_read(sam_tab, store = dstore)
  samp_n <- s_d$sam_names_matching$common_names[
    match(rownames(seqtab), s_d$sam_names_matching$raw_fastq)
  ]
  single <- !(duplicated(samp_n) | duplicated(samp_n, fromLast = TRUE))
  asv_tab <- otu_table(
    rename_samples(
      otu_table(seqtab[single, ], taxa_are_rows = FALSE),
      names_of_samples = samp_n[single]
    ),
    taxa_are_rows = FALSE
  )
  set.seed(target_seed("tax_tab"))
  tax_tab <- dada2::assignTaxonomy(
    seqtab,
    refFasta = here(
      "data/data_raw/refseq/dada2_format",
      paste0(seed_taxonomy_db, ".fasta")
    ),
    taxLevels = c(
      "Kingdom",
      "Phylum",
      "Class",
      "Order",
      "Family",
      "Genus",
      "Species"
    ),
    multithread = n_threads
  )
  d_asv <- add_dna_to_phyloseq(
    phyloseq(
      asv_tab,
      sam_tab,
      tax_table(as.matrix(tax_tab, dimnames = rownames(tax_tab)))
    )
  )
  attr(d_asv, "reads") <- tibble(
    arm = arm,
    fastq = basename(fnfs),
    raw = ShortRead::countFastq(fnfs)$records,
    filtered = vapply(derep_fs, \(x) sum(dada2::getUniques(x)), numeric(1)),
    merged = vapply(merged, \(x) sum(x$abundance), numeric(1)),
    final = rowSums(seqtab)
  )
  d_asv
}

# ---- ITSx ITS1 and the controls: pipelines/assign_taxo.R ----------------------
assignment_input <- function(d_asv) {
  real <- MiscMetabar::itsx_pq(
    d_asv,
    region = "ITS1",
    organism_groups = itsx_organism_groups,
    keep_undetected = TRUE,
    add_origin = TRUE,
    duplicated_seqs = "merge",
    nproc = n_threads,
    args_before_itsx = itsx_conda_prelude,
    verbose = FALSE
  )
  with_controls(
    real,
    tar_read(itsx_controls_ITS1, store = store),
    prop_fake = prop_fake,
    seed = targets_seed
  )
}

d_asv <- list(
  current = cached_rds("d_asv_current", denoise("current", Inf, 2)),
  rolling = cached_rds("d_asv_rolling", denoise("rolling", 8, 8))
)
input <- list(
  current = cached_rds("input_current", assignment_input(d_asv$current)),
  rolling = cached_rds("input_rolling", assignment_input(d_asv$rolling))
)

# ---- the replica must be exact ------------------------------------------------
stored_asv <- tar_read(d_asv, store = dstore)
stored_input <- tar_read(d_ITS1_asv_for_assignation, store = store)
seqs_of <- function(pq) stats::setNames(as.character(pq@refseq), taxa_names(pq))
stopifnot(
  "replica: refseq differs from store_dada2 d_asv" = identical(
    seqs_of(d_asv$current),
    seqs_of(stored_asv)
  ),
  "replica: otu_table differs from store_dada2 d_asv" = isTRUE(all.equal(
    as(otu_table(d_asv$current), "matrix"),
    as(otu_table(stored_asv), "matrix")
  )),
  "replica: assignment input differs from d_ITS1_asv_for_assignation" = identical(
    seqs_of(input$current),
    seqs_of(stored_input)
  )
)
message("✓ the current arm reproduces d_asv and d_ITS1_asv_for_assignation")

# ---- computations -------------------------------------------------------------
# Current arm: the stored computations (same input, same seeds). Rolling arm:
# the same computations on its own input, with the seed targets would give.
grid <- tidyr::expand_grid(
  arm = arms$arm,
  method = methods,
  db = zoom_test_dbs
) |>
  mutate(
    compute_name = compute_name_for(method, db, "ITS1_asv"),
    file = file.path(out_dir, paste0(arm, "__", compute_name, ".rds"))
  )
for (i in seq_len(nrow(grid))) {
  g <- grid[i, ]
  if (file.exists(g$file)) {
    next
  }
  computed <- if (g$arm == "current") {
    tar_read_raw(g$compute_name, store = store)
  } else {
    message(format(Sys.time(), "%H:%M"), " rolling ", g$compute_name)
    seed <- if (g$method == "dada2") {
      target_seed(g$compute_name)
    } else {
      config_sintax_seed
    }
    run_computation(input$rolling, g$method, g$db, seed)
  }
  saveRDS(computed, g$file)
}

# ---- units --------------------------------------------------------------------
truth <- list(
  current = build_mock_truth(d_asv$current, pauvert_cfg),
  rolling = build_mock_truth(d_asv$rolling, pauvert_cfg)
)
real_seqs <- lapply(input, function(pq) {
  s <- seqs_of(pq)
  s[!grepl(control_pattern, names(s))]
})
common <- intersect(real_seqs$current, real_seqs$rolling)
reads <- bind_rows(lapply(d_asv, attr, "reads"))
units <- tibble(
  arm = arms$arm,
  reads_final = vapply(d_asv, \(x) sum(sample_sums(x)), numeric(1)),
  asv = vapply(d_asv, ntaxa, integer(1)),
  asv_with_truth = vapply(truth, \(t) sum(!is.na(t$truth_depth)), integer(1)),
  asv_foreign = vapply(truth, \(t) sum(t$foreign), integer(1)),
  its1_units = lengths(real_seqs),
  its1_common = length(common)
)
write.csv(reads, file.path(out_dir, "reads.csv"), row.names = FALSE)
write.csv(units, file.path(out_dir, "units.csv"), row.names = FALSE)
cat(
  "\nReads per sample (raw -> filtered -> merged -> after chimeras and length)\n"
)
print(reads, width = 200)
cat("\nUnits per arm\n")
print(units, width = 200)

# ---- scores -------------------------------------------------------------------
# "all": every unit and control of the arm, as the pipeline scores it.
# "common": the real units present in both arms (same ITS1 sequence) and the
# external controls (the same in both arms); the shuffled controls differ
# between arms and are left out.
score_arm <- function(arm, reading) {
  pq_base <- input[[arm]]
  if (reading == "common") {
    keep <- names(real_seqs[[arm]])[real_seqs[[arm]] %in% common]
    keep <- c(keep, grep("^external_", taxa_names(pq_base), value = TRUE))
    pq_base <- prune_taxa(taxa_names(pq_base) %in% keep, pq_base)
  }
  rows <- which(grid$arm == arm)
  bind_rows(lapply(rows, function(i) {
    g <- grid[i, ]
    computed <- readRDS(g$file)
    pq <- derive_columns(input[[arm]], computed, g$method, g$db)
    pq <- prune_taxa(taxa_names(pq) %in% taxa_names(pq_base), pq)
    cols <- pauvert_vm$full_name[pauvert_vm$compute_name == g$compute_name]
    score_columns(pq, cols, truth[[arm]])
  })) |>
    mutate(arm = arm, reading = reading)
}
scores <- bind_rows(
  score_arm("current", "all"),
  score_arm("rolling", "all"),
  score_arm("current", "common"),
  score_arm("rolling", "common")
)
write.csv(scores, file.path(out_dir, "scores.csv"), row.names = FALSE)

seed_range_file <- here("data/snapshots/seed_spread/spread_default.csv")
seed_range <- if (file.exists(seed_range_file)) {
  read.csv(seed_range_file) |>
    select(method, db, tax_level, metrics, seed_range = range)
} else {
  tibble(
    method = character(),
    db = character(),
    tax_level = character(),
    metrics = character(),
    seed_range = numeric()
  )
}
comparison <- scores |>
  default_rows() |>
  endpoint_rows() |>
  select(reading, arm, method, db, tax_level, metrics, values) |>
  tidyr::pivot_wider(names_from = arm, values_from = values) |>
  mutate(delta = rolling - current) |>
  left_join(seed_range, by = c("method", "db", "tax_level", "metrics")) |>
  mutate(
    seed_range = ifelse(method %in% c("lca", "blastn"), 0, seed_range),
    beyond_seed = abs(delta) > seed_range
  ) |>
  arrange(reading, tax_level, metrics, method, db)
write.csv(
  comparison,
  file.path(out_dir, "comparison_default.csv"),
  row.names = FALSE
)
cat("\nOD 37 endpoints at default settings, rolling minus current\n")
print(
  mutate(
    comparison,
    across(c(current, rolling, delta, seed_range), \(x) round(x, 4))
  ),
  n = Inf,
  width = 200
)
cat("\nMean difference per endpoint and reading (8 columns each)\n")
print(
  comparison |>
    group_by(reading, tax_level, metrics) |>
    summarise(
      mean_delta = round(mean(delta), 4),
      better = sum(delta > 0),
      worse = sum(delta < 0),
      beyond_seed = sum(beyond_seed, na.rm = TRUE),
      .groups = "drop"
    ),
  n = Inf
)
