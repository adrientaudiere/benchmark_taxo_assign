library("conflicted")
library("MiscMetabar")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")


#tar_manifest(script="../arround_MiscMetabar/benchmarking_taxonomy.R")

values_map <- tidyr::expand_grid(
  method = c("dada2", "sintax"),
  #, "lca", "idtaxa"
  # min_boostrap = c(0.3, 0.5, 0.7, 0.8, 0.9),
  db = c("Maarjam", "Eukaryome"),
  db_filter = c("Global", "Fungi", "Glomeromycota"),
  # "Fungi_sparce", "Fungi_with_other", "Glomero_with_other"
  cutadapted_db = c("cut", "")
) |>
  mutate()

data_pq <- data_fungi_AM_mini

# Primer AM
fw_primer_AM <- "AAGCTCGTAGTTGAATTTCG" # AMV4.5NF Sato et al. (2005)
rev_primer_AM <- "CCCAACTATCCCTATTAATCAT" # AMDGR Sato et al. (2005)

# Primer ITS fungi
fw_primer_ITS_ <- "GCATCGATGAAGAACGCAGC"
rev_primer_ITS <- "TCCTCCGCTTATTGATATGC"

# Primer ITS fungi Pauvert et al. 2018
fw_primer_ITS_ <- "CTTGGTCATTTAGAGGAAGTAA" # ITS-1F Gardes and Bruns, 1993
rev_primer_ITS <- "GCTGCGTTCTTCATCGATGC" #ITS2 White 1990

tarchetypes::tar_map(
  unlist = FALSE,
  # Return a nested list from tar_map()
  values = values_map,
  names = "data_source",
  # Select columns from `values` for target names.
  tar_target(ref_fasta_db, if (cutadapted_db == "cut") {
    cutadapt_rm_primers_db(
      db_filter,
      output = paste0("data_intermediate/", db_filter, "cut"),
      primer_fw = fw_primer_ITS,
      primer_rev = rev_primer_ITS,
      return_file_path = TRUE
    )
  } else {
    db_filter
  }, format = file),
  tar_target(
    d_pq,
    add_new_taxonomy_pq(
      physeq = data_pq,
      method = method,
      ref_fasta = db_filter
    )
  )
)






