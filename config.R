# Shared constants for the benchmark pipelines. Sourced by pipelines/dada2.R,
# pipelines/assign_taxo.R, pipelines/cross_val.R, make_databases.R and the
# analysis notebooks.

# ITS-1F (Gardes & Bruns 1993) / ITS2 (White 1990), per Pauvert et al. 2018.
fw_primer_sequences <- "CTTGGTCATTTAGAGGAAGTAA"
rev_primer_sequences <- "GCTGCGTTCTTCATCGATGC"

# Glomeromycota primers — kept for reference, not used by the main pipeline.
fw_primer_AM  <- "AAGCTCGTAGTTGAATTTCG"    # AMV4.5NF, Sato et al. 2005
rev_primer_AM <- "CCCAACTATCCCTATTAATCAT"  # AMDGR,    Sato et al. 2005

n_threads   <- 4
seq_len_min <- 200
prop_fake   <- 0.5

# Reference databases benchmarked (decision 7: one SSU variant only). Used by
# R/values_map.R for both the assignment and the cross-validation grids.
db_list <- c(
  "Unite",
  "Unite_Fungi",
  "EUK_ITS_v2",
  "EUK_ITS_v2_Fungi",
  "EUK_ITS_v2_Fungi_cut",
  # "EUK_SSU_v2",
  # "EUK_SSU_v2_Fungi",
  # "EUK_SSU_v2_cut",
  "EUK_SSU_v2_Fungi_cut"
)

# Which simplification step each database represents (Q2 figures).
db_meta <- tibble::tribble(
  ~db,                    ~db_base,  ~simplification,
  "Unite",                "Unite",   "full",
  "Unite_Fungi",          "Unite",   "Fungi",
  "EUK_ITS_v2",           "EUK_ITS", "full",
  "EUK_ITS_v2_Fungi",     "EUK_ITS", "Fungi",
  "EUK_ITS_v2_Fungi_cut", "EUK_ITS", "Fungi+cut",
  "EUK_SSU_v2_Fungi_cut", "EUK_SSU", "Fungi+cut"
)

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
n_workers <- 3

targets_seed <- 22

cv_fold_number <- 10L
cv_fold_tested <- cv_fold_number  # publication value; lower to 2L for a smoke test
cv_max_seq     <- NULL            # publication value (full DB); set low (e.g. 200L) for a smoke test

# When TRUE, the pipelines use the mini_* variants of the reference databases
# for fast smoke-testing. Derived from the {targets} project name
# (_targets.yaml): the *_mini projects write to their own stores, so a smoke
# test can never overwrite a production store. Target names do not change.
#   Sys.setenv(TAR_PROJECT = "assign_taxo_mini"); targets::tar_make()
mini_db <- grepl("_mini$", Sys.getenv("TAR_PROJECT", unset = ""))
message("config.R: mini_db = ", mini_db,
        if (mini_db) " (smoke-test databases mini_*)" else " (full databases)")
