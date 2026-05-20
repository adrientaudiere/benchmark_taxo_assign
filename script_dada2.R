library("conflicted")
library("MiscMetabar")
library("targets")
library("tarchetypes")
library("here")
library("tibble")
library("tidyr")
library("autometric")

if (tar_active()) {
  log_start(path = "data/data_final/autometric_log_dada2.txt", seconds = 1)
}

here::i_am("script_dada2.R")
source(here("config.R"))
source(here("R/functions.R"))
lapply(list.files("~/Nextcloud/IdEst/Projets/MiscMetabar/R/", full.names = TRUE),
       source)

tar_option_set(seed = targets_seed)


tar_plan(
  tar_target(
    name = file_sam_data_csv,
    command = here("data/data_raw/metadata", sam_data_file_name),
    format = "file"
  ),
  tar_target(
    name = file_refseq_taxo,
    command = here("data/data_raw/refseq/", refseq_file_name),
    format = "file"
  ),
  tar_target(
    name = fastq_files_folder,
    command = here("data/data_raw/rawseq"),
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
      pattern_remove_fastq_files = "_R.*",
      prefix = "samp_"
    )
  ),

  #> Paired end analysis
  #> ———————————————————

  ## > Remove primers
  tar_target(
    cutadapt,
    {
      autometric::log_phase_set("cutadapt")
      cutadapt_remove_primers(
        path_to_fastq = fastq_files_folder,
        pattern = "fastq",
        primer_fw = fw_primer_sequences,
        primer_rev = rev_primer_sequences,
        folder_output = here("data/data_intermediate/seq_wo_primers/"),
        nproc = n_threads,
        return_file_path = TRUE,
        args_before_cutadapt = cutadapt_conda_prelude
      )
    },
    format = "file"
  ),
  tar_target(data_raw, {
    cutadapt
    list_fastq_files(path = here::here("data/data_intermediate/seq_wo_primers/"),
                     pattern_R1="_R1",
                     pattern_R2 = "_R2")
  }),

  ## > Classical dada2 pipeline
  tar_target(data_fnfs, data_raw$fnfs),
  tar_target(data_fnrs, data_raw$fnrs),
  ### Pre-filtered data with low stringency
  tar_target(
    filtered,
    {
      autometric::log_phase_set("filtered")
      filter_trim(
        output_fw = paste(
          getwd(),
          here("/data/data_intermediate/filterAndTrim_fwd"),
          sep = ""
        ),
        output_rev = paste(
          getwd(),
          here("/data/data_intermediate/filterAndTrim_rev"),
          sep = ""
        ),
        fw = data_fnfs,
        rev = data_fnrs,
        multithread = n_threads,
        compress = TRUE,
        trimLeft = 1,
        trimRight = 1
      )
    }
  ),

  ### Dereplicate fastq files
  tar_target(derep_fs, derepFastq(filtered[[1]]), format = "qs"),
  tar_target(derep_rs, derepFastq(filtered[[2]]), format = "qs"),
  ### Learns the error rates
  tar_target(
    err_fs,
    {
      autometric::log_phase_set("err_fs")
      learnErrors(derep_fs, multithread = n_threads)
    },
    format = "qs"
  ),
  tar_target(
    err_rs,
    {
      autometric::log_phase_set("err_rs")
      learnErrors(derep_rs, multithread = n_threads)
    },
    format = "qs"
  ),
  ### Make amplicon sequence variants
  tar_target(
    ddF,
    {
      autometric::log_phase_set("ddF")
      dada(derep_fs, err_fs, multithread = n_threads)
    },
    format = "qs"
  ),
  tar_target(
    ddR,
    {
      autometric::log_phase_set("ddR")
      dada(derep_rs, err_rs, multithread = n_threads)
    },
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
  ### Build a a table of ASV x Samples
  tar_target(seq_tab_Pairs, makeSequenceTable(merged_seq)),

  #> end Paired-end analysis
  #> ———————————————————————

  ## > Filtering sequences
  ### Remove chimera
  tar_target(seqtab_wo_chimera, chimera_removal_vs(seq_tab_Pairs)),
  ### Remove sequences based on length
  tar_target(seqtab, seqtab_wo_chimera[, nchar(colnames(seqtab_wo_chimera)) >= seq_len_min]),

  ## > Load sample data and rename samples
  tar_target(
    sam_tab,
    rename_samples(
      sample_data(s_d$sam_data),
      names_of_samples = s_d$sam_data$samples_names_common
    )
  ),
  tar_target(samp_n_otu_table,
             s_d$sam_names_matching$common_names[match(rownames(seqtab), s_d$sam_names_matching$raw_fastq)]),
  tar_target(asv_tab, otu_table(
    rename_samples(otu_table(seqtab[!(duplicated(samp_n_otu_table) |
                                   duplicated(samp_n_otu_table, fromLast = TRUE)), ],
                             taxa_are_rows = FALSE),
                   names_of_samples = samp_n_otu_table[!(duplicated(samp_n_otu_table) |
                                                        duplicated(samp_n_otu_table, fromLast = TRUE))]),
    taxa_are_rows = FALSE
  )),

  tar_target(
    tax_tab,
    {
      autometric::log_phase_set("tax_tab")
      assignTaxonomy(
        seqtab,
        refFasta = file_refseq_taxo,
        taxLevels = c(
          "Kingdom",
          "Phyla",
          "Class",
          "Order",
          "Family",
          "Genus",
          "Species"
        ),
        multithread = n_threads
      )
    }
  ),

  ## > Create the phyloseq object 'data_phyloseq' with
  ###   (i) table of asv,
  ###   ii) taxonomic table,
  ###   (iii) sample data and
  ###   (iv) references sequences

  tar_target(d_asv, add_dna_to_phyloseq(
    phyloseq(asv_tab, sam_tab, tax_table(
      as.matrix(tax_tab, dimnames = rownames(tax_tab))
    ))
  )),
  ## > Create post-clustering ASV into OTU using vsearch
  tar_target(d_vs, asv2otu(
    d_asv,
    method = "vsearch", tax_adjust = 0
  )),
  ## > Create post-clustering ASV into OTU using vsearch
  tar_target(d_idtaxa, asv2otu(
    d_asv,
    method = "clusterize", tax_adjust = 0
  )),
  ## > Clean post-clustering OTU using mumu
  tar_target(d_asv_mumu, mumu_pq(d_asv)$new_physeq),
  tar_target(d_vs_mumu, mumu_pq(d_vs)$new_physeq),
  tar_target(d_idtaxa_mumu, mumu_pq(d_idtaxa)$new_physeq),

  tar_target(track_df, track_wkflow(
    list(
      "Raw Forward sequences" = unlist(list_fastq_files(fastq_files_folder,
                                                        paired_end = FALSE, pattern_R1 = "_R1")),
      "Forward wo primers" = unlist(list_fastq_files(here::here("data/data_intermediate/seq_wo_primers/"),
                                                     paired_end = FALSE, pattern_R1 = "_R1")),
      "Forward sequences" = ddF,
      "Paired sequences" = seq_tab_Pairs,
      "Paired sequences without chimera" = seqtab_wo_chimera,
      "Paired sequences without chimera and longer than 200bp" = seqtab,
      "ASV denoising" = d_asv,
      "OTU after vsearch reclustering at 97%" = d_vs,
      "OTU after idtaxa reclustering at 97%" = d_idtaxa,
      "ASV mumu cleaning algorithm" = d_asv_mumu,
      "OTU vs after mumu cleaning algorithm" = d_vs_mumu,
      "OTU idtaxa after mumu cleaning algorithm" = d_idtaxa_mumu
    )
  )),

  tar_target(
    benchmark_costs_dada2,
    {
      track_df  # force aggregation after all pipeline targets complete
      log_df <- autometric::log_read(
        here("data/data_final/autometric_log_dada2.txt")
      )
      log_df |>
        dplyr::filter(!is.na(phase), phase != "") |>
        dplyr::group_by(phase) |>
        dplyr::summarise(
          wall_time_s      = as.numeric(max(time) - min(time)),
          peak_resident_mb = max(resident, na.rm = TRUE),
          mean_cpu_pct     = mean(cpu, na.rm = TRUE),
          n_samples        = dplyr::n(),
          .groups = "drop"
        )
    }
  )
)

