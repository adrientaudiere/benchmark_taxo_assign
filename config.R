# Shared constants for the benchmark pipelines. Sourced by
# script_assign_taxo_parallel.R. script_dada2.R still has inline copies; switch
# it to `source("config.R")` when next touching that pipeline.

# ITS-1F (Gardes & Bruns 1993) / ITS2 (White 1990), per Pauvert et al. 2018.
fw_primer_sequences <- "CTTGGTCATTTAGAGGAAGTAA"
rev_primer_sequences <- "GCTGCGTTCTTCATCGATGC"

# Glomeromycota primers — kept for reference, not used by the main pipeline.
fw_primer_AM  <- "AAGCTCGTAGTTGAATTTCG"    # AMV4.5NF, Sato et al. 2005
rev_primer_AM <- "CCCAACTATCCCTATTAATCAT"  # AMDGR,    Sato et al. 2005

n_threads   <- 4
seq_len_min <- 200
prop_fake   <- 0.5

refseq_file_name   <- "Unite_Fungi.fasta"
sam_data_file_name <- "sam_data.csv"
sample_col_name    <- "Sample_names"
fake_ref_fasta     <- "data/data_raw/fake_ref/fake_ref_asv_100.fasta"
taxo_mock_csv      <- "data/data_raw/metadata/taxo_mock.csv"

# Prelude that activates the cutadapt conda env. Override per-machine if your
# conda lives elsewhere.
cutadapt_conda_prelude <-
  "source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv && "

# Number of crew workers for parallel assignments. Lower this if your machine
# can't host n_workers * n_threads cores.
n_workers <- 1

targets_seed <- 22

cv_fold_number <- 10L
cv_fold_tested <- 2L    # smoke-test default; set to cv_fold_number for the publication run
cv_max_seq     <- 100L  # NULL for the full DB; set low (e.g. 200) for fast smoke-tests
