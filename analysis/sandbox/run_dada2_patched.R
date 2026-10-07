# One dada2::assignTaxonomy() call, with either the installed dada2 or the
# patched build of 1.40.0 (src/taxonomy.cpp line 85:
# lgk_probability[(size_t) g * n_kmers], see check_dada2_alias.R), on the
# pipeline's own input (Pauvert ITS1 ASV with its controls), minBoot = 0 and
# outputBootstraps = TRUE as R/assign_compute.R. Run one call per process:
#   Rscript analysis/sandbox/run_dada2_patched.R <build> <db> <out.rds>
# <build> is "installed" or the path of the library holding the patched build;
# <db> is a database of dada2_format/ or the path of a dada2-format fasta.

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3)
build <- args[[1]]
db <- args[[2]]
out <- args[[3]]
if (build != "installed") {
  .libPaths(c(build, .libPaths()))
}
suppressMessages(library(dada2))
message("dada2 ", packageVersion("dada2"), " from ", find.package("dada2"))

d <- targets::tar_read(
  ITS1_asv_compute_dada2__EUK_ITS_v2.1,
  store = "store_assign_taxo"
)
queries <- rownames(d$tax)
ref <- if (file.exists(db)) {
  db
} else {
  file.path("data/data_raw/refseq/dada2_format", paste0(db, ".fasta"))
}
t0 <- Sys.time()
res <- assignTaxonomy(
  queries,
  refFasta = ref,
  minBoot = 0,
  outputBootstraps = TRUE,
  multithread = 3
)
message(
  db,
  ": ",
  round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1),
  " min"
)
saveRDS(res, out)
