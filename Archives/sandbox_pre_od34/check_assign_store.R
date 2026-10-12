# Read-only checks of an assign_taxo store after tar_make(): run status, ITSx
# duplicates and the two inputs, completeness of d_all_taxo, assignments of
# the negative controls, benchmark_costs. Run from the project root:
#   TAR_PROJECT=assign_taxo_mini Rscript analysis/sandbox/check_assign_store.R
# or, in R: Sys.setenv(TAR_PROJECT = "assign_taxo"); source("analysis/sandbox/check_assign_store.R")
# Unset TAR_PROJECT defaults to assign_taxo_mini.

if (Sys.getenv("TAR_PROJECT") == "") {
  Sys.setenv(TAR_PROJECT = "assign_taxo_mini")
}
suppressMessages({
  library(targets)
  library(phyloseq)
  library(dplyr)
})
source("config.R")
source("R/values_map.R")

cat("== 1. Run status (", Sys.getenv("TAR_PROJECT"), ")\n")
m <- tar_meta(fields = c(error, warnings, time), targets_only = TRUE)
cat("targets:", nrow(m), "| errors:", sum(!is.na(m$error)), "| warnings:", sum(!is.na(m$warnings)), "\n")
for (k in which(!is.na(m$error) | !is.na(m$warnings))) {
  cat("  ", m$name[k], ":", substr(dplyr::coalesce(m$error[k], m$warnings[k]), 1, 120), "\n")
}
cat("outdated:", length(tar_outdated(callr_function = NULL)), "\n")

cat("\n== 2. ITSx duplicates and inputs\n")
d_asv <- tar_read(d_asv)
itsx <- tar_read(itsx_asv)
dropped <- tar_read(itsx_dropped_taxa)
its1 <- as.character(itsx$sequences)
cat("ASVs with ITS1 detected:", length(its1), "of", ntaxa(d_asv), "\n")
for (t in dropped) {
  twin <- setdiff(names(its1)[its1 == its1[[t]]], t)
  cat("dropped", t, "(", taxa_sums(d_asv)[[t]], "reads ), same ITS1 as", twin,
      "(", taxa_sums(d_asv)[twin], "reads )\n")
}
raw_in <- tar_read(d_asv_for_assignation)
itsx_in <- tar_read(d_asv_itsx_for_assignation)
is_asv <- \(pq) !grepl("^(fake|external)_", taxa_names(pq))
cat("d_asv_common:", ntaxa(tar_read(d_asv_common)), "| raw input:", ntaxa(raw_in),
    "| ITSx input:", ntaxa(itsx_in),
    "| same taxa:", setequal(taxa_names(raw_in), taxa_names(itsx_in)), "\n")
cat("fake_ / external_:", sum(grepl("^fake_", taxa_names(raw_in))), "/",
    sum(grepl("^external_", taxa_names(raw_in))), "\n")
cat("median ASV length, raw:", median(Biostrings::width(raw_in@refseq[is_asv(raw_in)])),
    "| ITSx:", median(Biostrings::width(itsx_in@refseq[is_asv(itsx_in)])),
    "| duplicated ITSx sequences:", anyDuplicated(as.character(itsx_in@refseq)) > 0, "\n")

cat("\n== 3. d_all_taxo\n")
d <- tar_read(d_all_taxo)
tt <- as.data.frame(unclass(tax_table(d)))
# lca writes "" for unassigned ranks; count them as NA like 01_load_and_clean.qmd.
assigned <- \(x) !is.na(x) & x != ""
vm <- build_values_map(dbs = db_list, mini_db = grepl("_mini$", Sys.getenv("TAR_PROJECT")))
has_cols <- vapply(vm$full_name, \(n) any(names(tt) %in% paste0(c("Kingdom", "Genus"), "_", n)), logical(1))
cat("taxa:", ntaxa(d), "| columns:", ncol(tt), "| rows:", nrow(vm),
    "| rows with no column (empty assignment):", sum(!has_cols), "\n")
for (n in vm$full_name[!has_cols]) cat("  empty:", n, "\n")
cat("rank values with a parenthesis:",
    sum(vapply(tt, \(x) any(grepl("(", x, fixed = TRUE)), logical(1))), "\n")
kind <- ifelse(grepl("^fake_", rownames(tt)), "fake",
               ifelse(grepl("^external_", rownames(tt)), "external", "asv"))
share <- \(col, who) if (col %in% names(tt)) mean(assigned(tt[[col]][kind == who])) else NA_real_
count <- \(col, who) if (col %in% names(tt)) sum(assigned(tt[[col]][kind == who])) else NA_integer_
print(as.data.frame(
  vm[has_cols, ] |>
    mutate(
      asv_genus = vapply(paste0("Genus_", full_name), count, integer(1), who = "asv"),
      fake_kingdom = vapply(paste0("Kingdom_", full_name), share, numeric(1), who = "fake"),
      external_kingdom = vapply(paste0("Kingdom_", full_name), share, numeric(1), who = "external")
    ) |>
    group_by(preprocess, method) |>
    summarise(
      rows = n(),
      median_asv_with_genus = median(asv_genus, na.rm = TRUE),
      fake_with_kingdom = round(mean(fake_kingdom, na.rm = TRUE), 3),
      external_with_kingdom = round(mean(external_kingdom, na.rm = TRUE), 3),
      .groups = "drop"
    )
))

cat("\n== 4. benchmark_costs\n")
bc <- tar_read(benchmark_costs)
cat("rows:", nrow(bc), "| computations:", n_distinct(bc$compute_name),
    "| rows with a measured wall time:", sum(!is.na(bc$wall_time_s)), "\n")
print(as.data.frame(bc |> distinct(compute_name, wall_time_s) |> arrange(desc(wall_time_s)) |> head(5)))
