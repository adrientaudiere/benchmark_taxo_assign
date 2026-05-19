library("conflicted")
library("MiscMetabar")
devtools::load_all("~/Nextcloud/IdEst/Projets/MiscMetabar/R/")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
library("autometric")


# unlink("data/data_final/autometric_log_assign_taxo.txt")
if (tar_active()) {
  log_start(path = "data/data_final/autometric_log_assign_taxo.txt", seconds = 1)
}

# tar_manifest(script="script_assign_taxo.R")
# tar_visnetwork(script = "script_assign_taxo.R", targets_only = T)
# tar_make(script = "script_assign_taxo.R", store = "store_assign_taxo")
# tar_destroy("all", store = "store_assign_taxo/")

## In a second terminal at the root of the R project
# library(targets)
# tar_poll(store="store_assign_taxo")


here::i_am("script_assign_taxo.R")
source("~/Nextcloud/IdEst/Projets/arround_MiscMetabar/comparpq/R/compare_taxo.R")
source("~/Nextcloud/IdEst/Projets/arround_MiscMetabar/comparpq/R/fake_creation.R")



# Primer AM
fw_primer_AM <- "AAGCTCGTAGTTGAATTTCG" # AMV4.5NF Sato et al. (2005)
rev_primer_AM <- "CCCAACTATCCCTATTAATCAT" # AMDGR Sato et al. (2005)

##  Primer ITS fungi
# fw_primer_ITS <- "GCATCGATGAAGAACGCAGC"
# rev_primer_ITS <- "TCCTCCGCTTATTGATATGC"

# Primer ITS fungi Pauvert et al. 2018
fw_primer_sequences <- "CTTGGTCATTTAGAGGAAGTAA" # ITS-1F Gardes and Bruns, 1993
rev_primer_sequences <- "GCTGCGTTCTTCATCGATGC" #ITS2 White 1990
n_threads <- 4
tar_option_set(seed = 22)



if (FALSE) {
  source(here::here("R/cutadapt_rm_primers_db.R"))
  database_cut <- values_map[values_map$cutadapted_db != "cut" &
                               grepl("EUK_", values_map$db) &
                               values_map$method != "lca", ]
  for (db in here::here(database_cut$db_path)) {
    cutadapt_rm_primers_db(
      db,
      output = paste0(dirname(db), "/", add_suffix_before_ext(db, "_cut")),
      primer_fw = fw_primer_sequences,
      primer_rev = rev_primer_sequences,
      return_file_path = TRUE
    )
  }
}

if (FALSE) {
  system("cd data/data_raw/refseq/sintax_format/;
for f in *; do
  head -n 10000 $f > mini_$f
done")

  system("cd data/data_raw/refseq/dada_format/;
for f in *; do
  head -n 10000 $f > mini_$f
done")
}


methods <- tidyr::expand_grid(
  method = c("dada2", "sintax", "lca"),
  min_bootstrap = c(0.4, 0.5, 0.6) # for idtaxa it defined the threshold
) |> full_join(tidyr::expand_grid(
  method = c("blastn"),
  vote_algorithm = c("rel_majority", "abs_majority", "unanimity"),
  nb_voting = 100, # nb_voting = c(NULL, 5, 20, 100),
  min_bootstrap = 0.5 # useless, only to make the put default value
 )) # |>
#   full_join(tidyr::expand_grid(
#     method = "idtaxa",
#     min_bootstrap = c(0.6)) # for idtaxa it defined the threshold
#   )

values_map <-
  tidyr::expand_grid(
    methods,
    db = c(
      # UNITE ——————————————————————————————
      "Unite",
      "Unite_Fungi",
      # Eukaryome ITS ——————————————————————————————
      "EUK_ITS_v1_9_3",
      #    "EUK_ITS_v1_9_3_Fungi",
      #"EUK_ITS_99_v1_9_3",
      # "EUK_ITS_99_v1_9_3_Fungi",
      #"EUK_ITS_v1_9_3_cut",
      #"EUK_ITS_v1_9_3_Fungi_cut",
      #    "EUK_ITS_99_v1_9_3_cut",
      #    "EUK_ITS_99_v1_9_3_Fungi_cut",
      # Eukaryome SSU ——————————————————————————————
      "EUK_SSU_v1_9_3",
      "EUK_SSU_v1_9_3_Fungi",
      #    "EUK_SSU_99_v1_9_3",
      #    "EUK_SSU_99_v1_9_3_Fungi",
      "EUK_SSU_v1_9_3_cut",
      "EUK_SSU_v1_9_3_Fungi_cut"#,
      #    "EUK_SSU_99_v1_9_3_cut",
      #    "EUK_SSU_99_v1_9_3_Fungi_cut",
    ),
  ) |>
  # mutate(db = paste0("mini_", db)) |> # comment to run the TRUE analysis
  mutate(cutadapted_db = ifelse(grepl("cut", db), "cut", ""))  |>
  mutate(db_filter = ifelse(grepl("Fungi", db), "Fungi", ""))  |>
  mutate(db_name = db) |>
  mutate(do_clean_pq = ifelse(method == "dada2", TRUE, FALSE)) |>
  mutate(db_path = paste0(ifelse(
    method == "dada2",
    paste0("data/data_raw/refseq/dada2_format/", db_name),
    paste0("data/data_raw/refseq/sintax_format/", db_name)
  ), ".fasta")) |>
  mutate(full_name = gsub("...NA...NA", "", paste0(method, "__", db, "___",
                                           min_bootstrap, "...",
                                           vote_algorithm, "...",
                                           nb_voting))) |>
  mutate(full_name = ifelse(row_number() == n(), paste0(full_name, "_all_taxo"), full_name)) |>
  mutate(previous_target = dplyr::lag(full_name)) |>
  mutate(previous_target = ifelse(row_number() == 1, "d_asv_for_assignation", previous_target)) |>
  mutate(previous_target = lapply(previous_target, as.symbol)) #|>
  #mutate(method = ifelse(method == "dada2", "dada2_2steps", method)) # uncomment to use both assignTaxonomy and assignSpecies



tar_plan(
  tar_target(d_asv, tar_read(d_asv, store = here::here("store_dada2"))),
  # tar_target(
  #   d_asv,
  #   subset_taxa(tar_read_raw("d_asv", store=here::here("store_dada2")),
  #               taxa_sums(tar_read_raw("d_asv", store=here::here("store_dada2")))>4000)
  # ), # case with less ASV to test for long-time new algo/db
  tar_target(
    file_taxo_mock,
    here("data/data_raw/metadata/taxo_mock.csv"),
    format = "file"
  ),
  tar_target(
    taxo_mock,
    read.csv(file_taxo_mock)  |>
      select(any_of(
        c(
          "Kingdom",
          "Phylum",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        )
      )) |>
      magrittr::set_rownames(read.csv(file_taxo_mock)$MockStrain)
  ),
  tar_target(
    d_asv_for_assignation_fake,
    add_shuffle_seq_pq(d_asv, prop_fake = 0.5)
  ),
  tar_target(
    d_asv_for_assignation,
    add_external_seq_pq(
      d_asv_for_assignation_fake,
      Biostrings::readDNAStringSet("data/data_raw/fake_ref/fake_ref_asv_100.fasta") #  Biostrings::readDNAStringSet("data/data_raw/fake_ref/fake_ref_asv_20.fasta")
    )
  ),
  tarchetypes::tar_eval(tar_target(
    full_name,
    add_new_taxonomy_pq(
      previous_target,
      method = method,
      ref_fasta = db_path,
      suffix = paste0("_", full_name),
      min_bootstrap = min_bootstrap,
      vote_algorithm = vote_algorithm,
      nb_voting = nb_voting
    )
  ), values_map)
)
