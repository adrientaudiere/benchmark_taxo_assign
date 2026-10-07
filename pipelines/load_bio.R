# Published OTUs of a biological-community dataset -> d_asv, without dada2
# (OD 30). For the D1c datasets that publish their representative sequences:
# the loader of R/bio_loaders.R turns the published tables into the phyloseq
# that pipelines/dada2_bio.R would otherwise denoise from fastq, and writes it
# to the same store (store_dada2_<dataset>) under the same name, so
# assign_taxo_<dataset> reads it unchanged. The project keeps the `dada2_`
# prefix for that reason: config.R derives the dataset and the store from it.
#
# No `d_vs`: the published OTUs are already species-level clusters, so the OTU
# input of assign_taxo is not run on these datasets (OD 30), and mumu is
# meaningless on them.
#
# Run with:
#   Sys.setenv(TAR_PROJECT = "dada2_gssp_air"); targets::tar_make()

library("targets")
library("here")

here::i_am("pipelines/load_bio.R")
source(here("config.R"))
source(here("R/bio_loaders.R"))
source(here("R/load_pqverse.R"))

bio_loaders <- list(
  gssp_air = load_gssp_air,
  korhonen_wood = load_korhonen_wood
)

if (is.null(bio_cfg) || !bio_dataset %in% names(bio_loaders)) {
  stop(
    "pipelines/load_bio.R needs a published-OTU dataset project. Set ",
    "TAR_PROJECT to one of: ",
    paste(paste0("dada2_", names(bio_loaders)), collapse = ", "),
    "."
  )
}

tar_option_set(seed = targets_seed)

list(
  tar_target(
    published_files,
    list.files(
      here(bio_cfg$raw_dir),
      pattern = "\\.csv$",
      full.names = TRUE,
      recursive = TRUE
    ),
    format = "file"
  ),
  tar_target(
    d_asv,
    {
      published_files
      bio_loaders[[bio_dataset]](
        here(bio_cfg$raw_dir),
        min_occur = bio_cfg$min_occur
      )
    }
  ),
  # The selection, with the published taxonomy of each OTU, so it can be read
  # without opening the store. A pipeline output, so not under data_raw/
  # (code review §3.4; there before 2026-10-07).
  tar_target(
    selected_units_csv,
    {
      out <- here(
        "data/data_final",
        paste0("selected_units_", bio_dataset, ".csv")
      )
      tax <- as.data.frame(as.matrix(d_asv@tax_table))
      utils::write.csv(
        data.frame(unit = rownames(tax), tax, row.names = NULL),
        out,
        row.names = FALSE
      )
      out
    },
    format = "file"
  ),
  tar_target(
    session_info,
    {
      selected_units_csv
      list(
        pqverse = pqverse_versions(),
        project = project_version(),
        session = sessioninfo::session_info()
      )
    }
  )
)
