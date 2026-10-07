suppressMessages({
  library(targets)
  library(phyloseq)
  library(dplyr)
  library(here)
})
source("config.R")
source("R/load_pqverse.R")
source("R/mock_truth.R")
load_pqverse(c("MiscMetabar", "dbpq"))

# Why is the Tedersoo Genus F1 flat at 0.52-0.59 for every method and
# database (chapter 01, 2026-10-02)? comparpq::tc_metrics_unit() counts a real
# unit without truth as FP when a method names it and FN when it leaves it NA,
# whatever the method. This script measures, on the denoised ASVs the truth is
# built from:
#   1. how many units and reads get a truth, Tedersoo against Pauvert;
#   2. why the others match no Sanger read: best identity and coverage to the
#      Sanger reads, length, N content;
#   3. what the methods call the units without truth (sintax x
#      Unite_all_20250219 at bootstrap 0.5, full-ITS ASV), and whether those
#      genera belong to the mock;
#   4. on the units with a truth, the share each method gets right at Genus.
#   Rscript analysis/sandbox/check_tedersoo_truth.R

cfgs <- split(mock_sanger_tables, mock_sanger_tables$dataset)
stores <- c(
  pauvert = "store_dada2",
  tedersoo_illumina = "store_dada2_tedersoo_illumina"
)

best_hits <- function(seqs, sanger) {
  q <- tempfile(fileext = ".fasta")
  d <- tempfile(fileext = ".fasta")
  o <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(q, d, o)))
  Biostrings::writeXStringSet(seqs, q, width = 20000)
  Biostrings::writeXStringSet(
    Biostrings::DNAStringSet(setNames(sanger$sequence, sanger$strain)),
    d,
    width = 20000
  )
  system2(
    dbpq::find_vsearch(),
    c(
      "--usearch_global",
      q,
      "--db",
      d,
      "--id",
      "0.5",
      "--maxaccepts",
      "0",
      "--maxrejects",
      "0",
      "--userfields",
      "query+target+id+qcov+tcov",
      "--userout",
      o,
      "--threads",
      "4",
      "--quiet"
    ),
    stdout = FALSE,
    stderr = FALSE
  )
  read.delim(
    o,
    header = FALSE,
    col.names = c("taxon", "strain", "id", "qcov", "tcov")
  ) |>
    group_by(taxon) |>
    slice_max(id, n = 1, with_ties = FALSE) |>
    ungroup()
}

units <- list()
for (ds in names(stores)) {
  d_asv <- tar_read(d_asv, store = stores[[ds]])
  truth <- build_mock_truth(d_asv, cfgs[[ds]])
  sanger <- read_mock_sanger(cfgs[[ds]])
  hits <- best_hits(d_asv@refseq, sanger)
  seqs <- d_asv@refseq
  units[[ds]] <- tibble(
    dataset = ds,
    taxon = taxa_names(d_asv),
    reads = taxa_sums(d_asv),
    length = Biostrings::width(seqs),
    n_N = as.vector(Biostrings::letterFrequency(seqs, "N")),
    truth_depth = truth$truth_depth[match(taxa_names(d_asv), truth$taxon)],
    truth_genus = truth$Genus[match(taxa_names(d_asv), truth$taxon)]
  ) |>
    left_join(hits, by = "taxon") |>
    mutate(
      class = case_when(
        !is.na(truth_depth) ~ "truth",
        is.na(id) ~ "no hit >= 50 %",
        id >= 99.5 ~ "id >= 99.5, cover < 99.5",
        id >= 97 ~ "id 97-99.5",
        id >= 90 ~ "id 90-97",
        TRUE ~ "id < 90"
      )
    )
}
u <- bind_rows(units)

message("\n1. Units and reads with a truth (denoised ASVs)")
u |>
  group_by(dataset) |>
  summarise(
    asv = n(),
    with_truth = sum(class == "truth"),
    share_units = round(mean(class == "truth"), 3),
    share_reads = round(sum(reads[class == "truth"]) / sum(reads), 3),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()

message("\n2. Units without truth, by best identity to a Sanger read")
u |>
  group_by(dataset, class) |>
  summarise(
    units = n(),
    median_reads = median(reads),
    reads = sum(reads),
    median_length = median(length),
    with_N = sum(n_N > 0),
    median_qcov = median(qcov),
    median_tcov = median(tcov),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()

ted <- u |> filter(dataset == "tedersoo_illumina")
message("\nTedersoo length of the ASVs (truth / no truth)")
print(tapply(ted$length, ted$class == "truth", summary))

message(
  "\n3. What sintax x Unite_all (bootstrap 0.5, full ASV) calls the Tedersoo units"
)
d_full <- tar_read_raw(
  "d_all_taxo_full_asv",
  store = "store_assign_taxo_tedersoo_illumina"
)
col <- grep(
  "^Genus.*sintax.*Unite_all_20250219___0\\.5",
  colnames(d_full@tax_table),
  value = TRUE
)
if (length(col) != 1) {
  col <- grep(
    "^Genus.*sintax.*Unite_all_20250219",
    colnames(d_full@tax_table),
    value = TRUE
  )
  message("candidate columns: ", paste(col, collapse = ", "))
  col <- col[1]
}
message("column: ", col)
sintax_genus <- setNames(
  as.character(d_full@tax_table[, col]),
  taxa_names(d_full)
)
mock_genera <- unique(genus_name(read_mock_sanger(cfgs$tedersoo_illumina)$name))
ted <- ted |>
  mutate(
    in_full_input = taxon %in% names(sintax_genus),
    sintax = sintax_genus[taxon],
    sintax_in_mock = sintax %in% mock_genera
  )
message(
  "ASVs in the full-ITS ASV input: ",
  sum(ted$in_full_input),
  " of ",
  nrow(ted)
)
ted |>
  filter(in_full_input) |>
  group_by(truth = class == "truth") |>
  summarise(
    units = n(),
    named = sum(!is.na(sintax)),
    named_mock_genus = sum(sintax_in_mock, na.rm = TRUE),
    named_other_genus = sum(!is.na(sintax) & !sintax_in_mock),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()
message(
  "\nGenera given to units without truth that are not mock genera (top 25)"
)
ted |>
  filter(in_full_input, class != "truth", !is.na(sintax), !sintax_in_mock) |>
  count(sintax, sort = TRUE) |>
  head(25) |>
  as.data.frame() |>
  print()
message("\nUnits without truth named with a mock genus: best Sanger identity")
ted |>
  filter(in_full_input, class != "truth", sintax_in_mock) |>
  count(class) |>
  as.data.frame() |>
  print()

message(
  "\n4. Units with a Genus truth, full ASV: right / wrong / NA, sintax x Unite_all"
)
ted |>
  filter(in_full_input, !is.na(truth_genus)) |>
  mutate(
    call = case_when(
      is.na(sintax) ~ "NA",
      sintax == truth_genus ~ "right",
      TRUE ~ "wrong"
    )
  ) |>
  count(call) |>
  as.data.frame() |>
  print()
message("Wrong pairs (truth -> sintax)")
ted |>
  filter(
    in_full_input,
    !is.na(truth_genus),
    !is.na(sintax),
    sintax != truth_genus
  ) |>
  count(truth_genus, sintax, sort = TRUE) |>
  as.data.frame() |>
  print()
message("(the scoring also accepts the GBIF accepted name: see 7)")

# 5. Genus F1 of the 36 default full-ITS ASV columns, controls left aside, on
# all real units (the chapter figure), without the non-mock units (best
# Sanger identity < 90 % or no hit), and on the units with a Genus truth.
truth <- build_mock_truth(
  tar_read(d_asv, store = stores[["tedersoo_illumina"]]),
  cfgs$tedersoo_illumina
)
tt <- as.data.frame(as(d_full@tax_table, "matrix"))
default_cols <- grep(
  paste0(
    "^Genus_full_asv_(dada2|sintax)__.*___0\\.5$|",
    "^Genus_full_asv_lca__.*___0\\.51$|",
    "^Genus_full_asv_blastn__.*rel_majority\\.\\.\\.100\\.\\.\\.95$"
  ),
  colnames(tt),
  value = TRUE
)
real <- ted$taxon[ted$in_full_input]
tr <- truth[match(real, truth$taxon), ]
non_mock <- real %in% ted$taxon[ted$class %in% c("id < 90", "no hit >= 50 %")]
has_truth <- !is.na(tr$Genus)
f1 <- function(keep, v) {
  ok <- has_truth &
    !is.na(v) &
    (v == tr$Genus | v %in% tr$Genus_accepted & v == tr$Genus_accepted)
  ok[is.na(ok)] <- FALSE
  2 *
    sum(keep & ok) /
    (2 * sum(keep & ok) + sum(keep & !is.na(v) & !ok) + sum(keep & is.na(v)))
}
cols_f1 <- lapply(default_cols, function(cl) {
  v <- tt[real, cl]
  v[grepl("[Ii]ncertae_sedis|^unclassified|^$", v)] <- NA
  tibble(
    algo = sub("^Genus_full_asv_([a-z0-9]+)__.*$", "\\1", cl),
    all_real = f1(rep(TRUE, length(real)), v),
    without_non_mock = f1(!non_mock, v),
    truth_only = f1(has_truth, v),
    named_non_mock = mean(!is.na(v[non_mock]))
  )
}) |>
  bind_rows()
message(
  "\n5. Genus F1, 36 default full-ASV columns (",
  sum(non_mock),
  " non-mock units)"
)
cols_f1 |>
  group_by(algo) |>
  summarise(
    across(everything(), \(x) paste(round(range(x), 3), collapse = "-")),
    .groups = "drop"
  ) |>
  as.data.frame() |>
  print()

# 6. Identity histogram of the units without truth: where to cut non-mock.
message("\n6. Best Sanger identity of the units without truth")
u |>
  filter(class != "truth") |>
  mutate(
    bin = cut(id, c(50, 70, 80, 90, 95, 97, 99, 99.5, 100.01), right = FALSE)
  ) |>
  count(dataset, bin) |>
  tidyr::pivot_wider(names_from = dataset, values_from = n, values_fill = 0) |>
  as.data.frame() |>
  print()

# 7. Accepted names of the split genera, and the units at >= 99.5 % identity
# whose cover fails: where they sit on their Sanger read.
message("\n7. GBIF accepted genus of the Cortinarius / Calocybe strains")
truth |>
  filter(Genus %in% c("Cortinarius", "Calocybe")) |>
  distinct(Genus, Species, Genus_accepted) |>
  count(Genus, Genus_accepted) |>
  as.data.frame() |>
  print()
part <- ted |> filter(class == "id >= 99.5, cover < 99.5")
sanger_t <- read_mock_sanger(cfgs$tedersoo_illumina)
part$start_on_sanger <- mapply(
  function(taxon, strain) {
    s <- sanger_t$sequence[sanger_t$strain == strain]
    regexpr(
      substr(as.character(d_asv@refseq[[taxon]]), 1, 20),
      s,
      fixed = TRUE
    )[[1]]
  },
  part$taxon,
  part$strain
)
part$sanger_length <- nchar(sanger_t$sequence[match(
  part$strain,
  sanger_t$strain
)])
message(
  "units at >= 99.5 % with partial cover (start_on_sanger -1: the ASV starts before the read)"
)
print(as.data.frame(part[, c(
  "taxon",
  "strain",
  "reads",
  "length",
  "id",
  "qcov",
  "tcov",
  "start_on_sanger",
  "sanger_length"
)]))

saveRDS(
  u,
  file.path(Sys.getenv("SCRATCH", tempdir()), "tedersoo_truth_units.rds")
)
