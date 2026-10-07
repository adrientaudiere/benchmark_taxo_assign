source("analysis/00_setup.R")
source(here("R/assignment_helpers.R"))
source(here("R/score_assignments.R"))

# Decision 26 (ROADMAP §2, "Wood, then air ... against the mocks"): the
# three rates a dataset without truth can give, on the biological datasets
# (wood and air, ITS2 ASVs) against the mocks (Pauvert ITS1, Tedersoo ITS1 /
# ITS2 / full ASVs), per method x database at the default settings, labels
# cleaned with clean_taxo() as the scoring does. First look 2026-09-30 (wood
# against Pauvert), extended 2026-10-05 to Tedersoo and air; Tedersoo ITS2 is
# the mock with the region of wood and air.
#   NA_real   share of the real units left NA at Genus
#   NA_fake   share of the `fake_` controls left NA at Genus (1 is right)
#   ext_fungi share of the `external_` controls named Fungi at Kingdom
#             (always an error)
# The question of decision 26: do the mocks behave like real communities in
# their NA rates, i.e. are they representative for that part?

tedersoo_store <- here("store_assign_taxo_tedersoo_illumina")
inputs <- tibble::tribble(
  ~dataset        , ~target               , ~store                                  ,
  "pauvert ITS1"  , "d_all_taxo_ITS1_asv" , store_assign_taxo                       ,
  "tedersoo ITS1" , "d_all_taxo_ITS1_asv" , tedersoo_store                          ,
  "tedersoo ITS2" , "d_all_taxo_ITS2_asv" , tedersoo_store                          ,
  "tedersoo full" , "d_all_taxo_full_asv" , tedersoo_store                          ,
  "wood ITS2"     , "d_all_taxo_ITS2_asv" , here("store_assign_taxo_korhonen_wood") ,
  "air ITS2"      , "d_all_taxo_ITS2_asv" , here("store_assign_taxo_gssp_air")
)
mocks <- c("pauvert ITS1", "tedersoo ITS1", "tedersoo ITS2", "tedersoo full")
biological <- c("wood ITS2", "air ITS2")

is_default <- function(algo, method_db) {
  case_when(
    algo %in% c("dada2", "sintax") ~ grepl("___0\\.5$", method_db),
    algo == "lca" ~ grepl("___0\\.51$", method_db),
    algo == "blastn" ~
      grepl("\\.\\.\\.rel_majority\\.\\.\\.100\\.\\.\\.95$", method_db),
    .default = FALSE
  )
}

rates_of <- function(dataset, target, store) {
  d <- tar_read_raw(target, store = store) |> clean_taxo()
  tt <- as.data.frame(as(tax_table(d), "matrix"), stringsAsFactors = FALSE)
  class <- ifelse(
    grepl("^external_", rownames(tt)),
    "external",
    ifelse(grepl("^fake_", rownames(tt)), "fake", "real")
  )
  genus <- grep("^Genus_.*__", colnames(tt), value = TRUE)
  tibble(
    dataset = dataset,
    method_db = sub("^Genus_", "", genus),
    algo = sub("__.*$", "", sub("^Genus_[^_]+_[a-z]+_", "", genus)),
    db = assignment_db(sub("^Genus_", "", genus)),
    NA_real = vapply(
      genus,
      \(x) mean(is.na(tt[[x]][class == "real"])),
      numeric(1)
    ),
    NA_fake = vapply(
      genus,
      \(x) mean(is.na(tt[[x]][class == "fake"])),
      numeric(1)
    ),
    ext_fungi = vapply(
      sub("^Genus_", "Kingdom_", genus),
      \(x) mean(tt[[x]][class == "external"] %in% "Fungi"),
      numeric(1)
    )
  ) |>
    filter(is_default(algo, method_db)) |>
    left_join(select(db_meta, db, simplification), by = "db")
}

rates <- bind_rows(Map(rates_of, inputs$dataset, inputs$target, inputs$store))
stopifnot(all(count(rates, dataset)$n == 36))

# The Pauvert rates must be chapter 01's.
ch01 <- readRDS("analysis/_cache/res_comp_tax.rds") |>
  filter(
    dataset == "pauvert",
    unit == "ASV",
    metrics %in% c("NA_real", "NA_fake")
  ) |>
  filter(tax_level == "Genus") |>
  distinct(method_db, metrics, values) |>
  pivot_wider(
    names_from = metrics,
    values_from = values,
    names_prefix = "ch01_"
  )
check <- inner_join(
  filter(rates, dataset == "pauvert ITS1"),
  ch01,
  by = "method_db"
)
stopifnot(nrow(check) == 36)
message(
  "NA_fake matches chapter 01 on Pauvert ITS1 ASV: ",
  all(abs(check$NA_fake - check$ch01_NA_fake) < 1e-9),
  "; NA_real differs only by the units left out below their truth depth: max gap ",
  round(max(abs(check$NA_real - check$ch01_NA_real)), 3)
)

show <- function(x) print(as.data.frame(x), row.names = FALSE)

for (metric in c("NA_real", "NA_fake", "ext_fungi")) {
  message(
    "\n== ",
    metric,
    ", mean over the three databases of each simplification"
  )
  rates |>
    group_by(simplification, algo, dataset) |>
    summarise(v = round(mean(.data[[metric]]), 2), .groups = "drop") |>
    pivot_wider(names_from = dataset, values_from = v) |>
    arrange(simplification, algo) |>
    show()
}

message("\n== Biological dataset against each mock, over the 36 columns")
message("gap = mean(biological - mock); rho = Spearman rank correlation")
wide <- rates |>
  mutate(col = paste(algo, db)) |>
  select(dataset, col, NA_real, NA_fake, ext_fungi)
pairs <- tidyr::expand_grid(
  biological = biological,
  mock = mocks,
  metric = c("NA_real", "NA_fake", "ext_fungi")
)
pairs$gap <- NA_real_
pairs$rho <- NA_real_
for (i in seq_len(nrow(pairs))) {
  a <- filter(wide, dataset == pairs$biological[i])
  b <- filter(wide, dataset == pairs$mock[i])
  ab <- inner_join(a, b, by = "col", suffix = c("_bio", "_mock"))
  stopifnot(nrow(ab) == 36)
  va <- ab[[paste0(pairs$metric[i], "_bio")]]
  vb <- ab[[paste0(pairs$metric[i], "_mock")]]
  pairs$gap[i] <- mean(va - vb)
  pairs$rho[i] <- suppressWarnings(cor(va, vb, method = "spearman"))
}
pairs |>
  mutate(gap = round(gap, 2), rho = round(rho, 2)) |>
  pivot_wider(names_from = metric, values_from = c(gap, rho)) |>
  show()

message("\n== NA_real per method, mean over the nine databases")
rates |>
  group_by(algo, dataset) |>
  summarise(v = round(mean(NA_real), 2), .groups = "drop") |>
  pivot_wider(names_from = dataset, values_from = v) |>
  show()
