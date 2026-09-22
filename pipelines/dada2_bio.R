# DADA2 denoising of a biological-community dataset -> d_asv (plus the vsearch
# OTU variant). Same flow as pipelines/dada2.R, with three differences
# (developer decisions 2026-09-16, ROADMAP D1c):
#   - the raw folder, metadata file, primers and fastq name patterns come from
#     config.R::bio_datasets, selected by TAR_PROJECT;
#   - the `clusterize` post-clustering branch is dropped (it was named
#     `d_idtaxa` in pipelines/dada2.R although it calls method = "clusterize");
#     `d_vs` stays, since Q5 rests on the vsearch OTUs;
#   - no seed taxonomy is needed here: the benchmark assigns in
#     store_assign_taxo_<dataset>, and these datasets have no mock truth to seed
#     from. An empty tax_table would break add_shuffle_seq_pq(), so the seed
#     assignment is kept (see the tax_tab target).
#
# Run with:
#   Sys.setenv(TAR_PROJECT = "dada2_tedersoo_illumina"); targets::tar_make()
#
# Requires cutadapt in the `cutadaptenv` conda env (config.R) and vsearch on
# PATH. Single-process pipeline (no crew controller).

library("conflicted")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
# dada2 is an Import (not a Depends) of MiscMetabar; the bare derepFastq(),
# learnErrors(), dada(), mergePairs(), makeSequenceTable(), assignTaxonomy()
# calls below need it attached.
library("dada2")

here::i_am("pipelines/dada2_bio.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/create_fake_pq_from_refseq.R"))
source(here("R/autometric_helpers.R"))
load_pqverse("MiscMetabar")

if (is.null(bio_cfg)) {
  stop(
    "pipelines/dada2_bio.R needs a dataset project. Set TAR_PROJECT to one of: ",
    paste(paste0("dada2_", bio_datasets$dataset), collapse = ", "),
    "."
  )
}

tar_option_set(seed = targets_seed)

# One autometric log file per phase and per run, in a folder of its own so a
# dataset run never overwrites the mock-community costs.
autometric_dir_bio <- here(
  "data/data_final/autometric",
  paste0("dada2_", bio_dataset)
)

tar_plan(
  tar_target(
    name = file_sam_data_csv,
    command = here("data/data_raw/metadata", bio_cfg$sam_data),
    format = "file"
  ),
  tar_target(
    name = file_refseq_taxo,
    command = here(
      "data/data_raw/refseq/dada2_format",
      paste0(seed_taxonomy_db, ".fasta")
    ),
    format = "file"
  ),
  tar_target(
    name = fastq_files_folder,
    command = here(bio_cfg$raw_dir),
    format = "file"
  ),

  #> Match samples names from fastq files and metadata sam_data
  #> ———————————————————
  tar_target(
    s_d,
    sam_data_matching_names(
      path_sam_data = file_sam_data_csv,
      path_raw_seq = fastq_files_folder,
      sample_col_name = sample_col_name,
      pattern_remove_fastq_files = bio_cfg$name_strip,
      prefix = "samp_"
    )
  ),

  #> Paired end analysis
  #> ———————————————————

  ## > Remove primers
  tar_target(
    cutadapt,
    with_autometric(
      "cutadapt",
      cutadapt_remove_primers(
        path_to_fastq = fastq_files_folder,
        pattern = "fastq",
        # cutadapt_remove_primers() has its own pattern_R1 / pattern_R2
        # (defaults "_R1" / "_R2") which it passes on to list_fastq_files().
        # Without them it never sees the ENA's `_1` / `_2` files and stops with
        # "None file in the folder ... match the pattern_R1 _R1".
        pattern_R1 = bio_cfg$pattern_r1,
        pattern_R2 = bio_cfg$pattern_r2,
        primer_fw = bio_cfg$fw_primer,
        primer_rev = bio_cfg$rev_primer,
        folder_output = here(
          "data/data_intermediate",
          paste0("seq_wo_primers_", bio_dataset)
        ),
        nproc = n_threads,
        return_file_path = TRUE,
        args_before_cutadapt = cutadapt_conda_prelude
      ),
      dir = autometric_dir_bio
    ),
    format = "file"
  ),
  tar_target(data_raw, {
    cutadapt
    list_fastq_files(
      path = here(
        "data/data_intermediate",
        paste0("seq_wo_primers_", bio_dataset)
      ),
      pattern_R1 = bio_cfg$pattern_r1,
      pattern_R2 = bio_cfg$pattern_r2
    )
  }),

  ## > Classical dada2 pipeline
  tar_target(data_fnfs, data_raw$fnfs),
  tar_target(data_fnrs, data_raw$fnrs),
  ### Pre-filtered data with low stringency
  tar_target(
    filtered,
    with_autometric(
      "filtered",
      filter_trim(
        output_fw = here(
          "data/data_intermediate",
          paste0("filterAndTrim_fwd_", bio_dataset)
        ),
        output_rev = here(
          "data/data_intermediate",
          paste0("filterAndTrim_rev_", bio_dataset)
        ),
        fw = data_fnfs,
        rev = data_fnrs,
        multithread = n_threads,
        compress = TRUE,
        trimLeft = 1,
        trimRight = 1
      ),
      dir = autometric_dir_bio
    )
  ),

  ### Dereplicate fastq files
  tar_target(derep_fs, derepFastq(filtered[[1]]), format = "qs"),
  tar_target(derep_rs, derepFastq(filtered[[2]]), format = "qs"),
  ### Learns the error rates
  tar_target(
    err_fs,
    with_autometric(
      "err_fs",
      learnErrors(derep_fs, multithread = n_threads),
      dir = autometric_dir_bio
    ),
    format = "qs"
  ),
  tar_target(
    err_rs,
    with_autometric(
      "err_rs",
      learnErrors(derep_rs, multithread = n_threads),
      dir = autometric_dir_bio
    ),
    format = "qs"
  ),
  ### Make amplicon sequence variants
  tar_target(
    ddF,
    with_autometric(
      "ddF",
      dada(derep_fs, err_fs, multithread = n_threads),
      dir = autometric_dir_bio
    ),
    format = "qs"
  ),
  tar_target(
    ddR,
    with_autometric(
      "ddR",
      dada(derep_rs, err_rs, multithread = n_threads),
      dir = autometric_dir_bio
    ),
    format = "qs"
  ),
  ### Merge paired sequences
  tar_target(
    merged_seq,
    mergePairs(
      dadaF = ddF,
      dadaR = ddR,
      derepF = derep_fs,
      derepR = derep_rs
    ),
    format = "qs"
  ),
  ### Build a table of ASV x Samples
  tar_target(seq_tab_Pairs, makeSequenceTable(merged_seq)),

  #> end Paired-end analysis
  #> ———————————————————————

  ## > Filtering sequences
  ### Remove chimera
  tar_target(seqtab_wo_chimera, chimera_removal_vs(seq_tab_Pairs)),
  ### Remove sequences based on length
  tar_target(
    seqtab,
    seqtab_wo_chimera[, nchar(colnames(seqtab_wo_chimera)) >= seq_len_min]
  ),

  ## > Load sample data and rename samples
  tar_target(
    sam_tab,
    rename_samples(
      sample_data(s_d$sam_data),
      names_of_samples = s_d$sam_data$samples_names_common
    )
  ),
  tar_target(
    samp_n_otu_table,
    s_d$sam_names_matching$common_names[match(
      rownames(seqtab),
      s_d$sam_names_matching$raw_fastq
    )]
  ),
  tar_target(
    asv_tab,
    otu_table(
      rename_samples(
        otu_table(
          seqtab[
            !(duplicated(samp_n_otu_table) |
              duplicated(samp_n_otu_table, fromLast = TRUE)),
          ],
          taxa_are_rows = FALSE
        ),
        names_of_samples = samp_n_otu_table[
          !(duplicated(samp_n_otu_table) |
            duplicated(samp_n_otu_table, fromLast = TRUE))
        ]
      ),
      taxa_are_rows = FALSE
    )
  ),

  ## > Seed taxonomy. The benchmark re-assigns in store_assign_taxo_<dataset>;
  ## this only fills the tax_table, which add_shuffle_seq_pq() requires.
  tar_target(
    tax_tab,
    with_autometric(
      "tax_tab",
      assignTaxonomy(
        seqtab,
        refFasta = file_refseq_taxo,
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
      ),
      dir = autometric_dir_bio
    )
  ),

  ## > Create the phyloseq object 'd_asv'
  tar_target(
    d_asv,
    add_dna_to_phyloseq(
      phyloseq(
        asv_tab,
        sam_tab,
        tax_table(
          as.matrix(tax_tab, dimnames = rownames(tax_tab))
        )
      )
    )
  ),
  ## > Create post-clustering ASV into OTU using vsearch (Q5)
  tar_target(
    d_vs,
    asv2otu(
      d_asv,
      method = "vsearch",
      tax_adjust = 0
    )
  ),
  ## > Clean post-clustering OTU using mumu
  tar_target(d_asv_mumu, mumu_pq(d_asv)$new_physeq),
  tar_target(d_vs_mumu, mumu_pq(d_vs)$new_physeq),

  tar_target(
    track_df,
    track_wkflow(
      list(
        "Raw Forward sequences" = unlist(list_fastq_files(
          fastq_files_folder,
          paired_end = FALSE,
          pattern_R1 = bio_cfg$pattern_r1
        )),
        "Forward wo primers" = unlist(list_fastq_files(
          here(
            "data/data_intermediate",
            paste0("seq_wo_primers_", bio_dataset)
          ),
          paired_end = FALSE,
          pattern_R1 = bio_cfg$pattern_r1
        )),
        "Forward sequences" = ddF,
        "Paired sequences" = seq_tab_Pairs,
        "Paired sequences without chimera" = seqtab_wo_chimera,
        "Paired sequences without chimera and longer than seq_len_min" = seqtab,
        "ASV denoising" = d_asv,
        "OTU after vsearch reclustering at 97%" = d_vs,
        "ASV mumu cleaning algorithm" = d_asv_mumu,
        "OTU vs after mumu cleaning algorithm" = d_vs_mumu
      )
    )
  ),

  tar_target(
    benchmark_costs_dada2,
    {
      track_df # aggregate only once every logged phase has run
      read_autometric_dir(autometric_dir_bio) |>
        summarise_autometric_costs()
    }
  ),

  tar_target(
    session_info,
    {
      track_df
      list(
        pqverse = pqverse_versions(),
        session = sessioninfo::session_info()
      )
    }
  )
)
