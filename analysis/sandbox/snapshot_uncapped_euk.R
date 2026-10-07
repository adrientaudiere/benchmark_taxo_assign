# Snapshot of the raw computations against the full EUK_ITS_v2.1 *before* its
# 5 000 bp record cap (developer, 2026-10-01): the next run overwrites them,
# and the fake_-only test of the cap (check_length_cap_fakes.R) compares each
# control's uncapped answer with its capped one. dada2 results are those of
# the patched build. One rds per target, plus an index.
#   Rscript analysis/sandbox/snapshot_uncapped_euk.R

out_dir <- "data/snapshots/uncapped_EUK_ITS_v2.1_2026-10-01"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stores <- c(
  pauvert = "store_assign_taxo",
  korhonen_wood = "store_assign_taxo_korhonen_wood",
  tedersoo_illumina = "store_assign_taxo_tedersoo_illumina"
)
pattern <- "^(ITS1|ITS2|full)_(asv|otu)_compute_(dada2|sintax|lca|blastn)__EUK_ITS_v2\\.1$"

index <- do.call(
  rbind,
  lapply(names(stores), function(ds) {
    meta <- targets::tar_meta(
      store = stores[[ds]],
      fields = c(name, time, error)
    )
    hit <- meta$name[grepl(pattern, meta$name) & is.na(meta$error)]
    hit <- hit[targets::tar_exist_objects(hit, store = stores[[ds]])]
    do.call(
      rbind,
      lapply(hit, function(target) {
        file <- file.path(out_dir, paste0(ds, "__", target, ".rds"))
        saveRDS(targets::tar_read_raw(target, store = stores[[ds]]), file)
        data.frame(
          dataset = ds,
          target = target,
          built = format(meta$time[meta$name == target]),
          file = basename(file)
        )
      })
    )
  })
)
write.csv(index, file.path(out_dir, "index.csv"), row.names = FALSE)
print(table(
  index$dataset,
  sub(".*_compute_([a-z0-9]+)__.*", "\\1", index$target)
))
