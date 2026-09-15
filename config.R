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

# Minimum BLAST query cover (%) for blastn assignments and cross-validation
# (decision 18). The ASVs start with 30-45 bp of the 18S end that ITS
# references do not contain, so the assign_blastn() default of 95 rejects most
# true hits (critique_fable/07_results_validity.md, B18).
blastn_min_cover <- 80

# ---- Reference databases ----------------------------------------------------
# How they are built and how to move to a new release:
# docs/reference_databases.md. In short: edit the two tables below (new
# versioned names), run make_databases.R, then rerun the pipelines.

# General FASTA releases, downloaded by make_databases.R with
# dbpq::download_unite_db() / dbpq::download_eukaryome_db() and stored as
# data/data_raw/refseq/sources/<source>.fasta. `source` is the versioned name
# carried by every file and target derived from the release.
#   UNITE "eukaryotes 2" = all eukaryotes WITH singletons (s_all archive);
#   UNITE "eukaryotes"   = all eukaryotes without singletons (all archive).
#   UNITE DOIs resolve to an HTML page: the url comes from the PlutoF API
#   (https://api.plutof.ut.ee/v1/public/dois/?identifier=<doi>).
#   EUKARYOME urls are listed on https://eukaryome.org/generalfasta/.
reference_sources <- tibble::tribble(
  ~source,                ~provider,   ~release,     ~doi,                   ~url,
  "Unite_s_all_20250219", "unite",     "19.02.2025", "10.15156/BIO/3301232", "https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz",
  "Unite_all_20250219",   "unite",     "19.02.2025", "10.15156/BIO/3301231", "https://s3.hpc.ut.ee/plutof-public/original/e861a3d6-54f4-42dc-882a-5f129beac39a.tgz",
  "EUK_ITS_v2.1",         "eukaryome", "2.1",        NA_character_,          "https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v2.1.zip"
)

# Databases benchmarked: one source and one simplification each.
#   "full"      = the whole release;
#   "Fungi"     = records whose kingdom is exactly Fungi;
#   "Fungi+cut" = "Fungi", then trimmed to the amplicon region with cutadapt
#                 (config primers; records without primer sites kept whole).
# Decision 7 (revised 2026-09-15): no SSU database. ITS1 is not part of the
# 18S gene, so an SSU record only holds the 18S start of the ASVs.
benchmark_dbs <- tibble::tribble(
  ~source,                ~simplification,
  "Unite_s_all_20250219", "full",
  "Unite_s_all_20250219", "Fungi",
  "Unite_all_20250219",   "full",
  "Unite_all_20250219",   "Fungi",
  "EUK_ITS_v2.1",         "full",
  "EUK_ITS_v2.1",         "Fungi",
  "EUK_ITS_v2.1",         "Fungi+cut"
)
db_suffix <- c(full = "", Fungi = "_Fungi", "Fungi+cut" = "_Fungi_cut")
benchmark_dbs$db <- paste0(benchmark_dbs$source, db_suffix[benchmark_dbs$simplification])

# Used by R/values_map.R (assignment and cross-validation grids) and by the
# Q2 figures (which simplification step each database represents).
db_list <- benchmark_dbs$db
db_meta <- tibble::tibble(
  db             = benchmark_dbs$db,
  db_base        = benchmark_dbs$source,
  simplification = benchmark_dbs$simplification
)

# Seed taxonomy of the DADA2 pipeline (the benchmark itself re-assigns).
seed_taxonomy_db <- "Unite_s_all_20250219_Fungi"
# Release whose non-Fungi records feed the external negative controls.
fake_ref_source <- "Unite_s_all_20250219"
# (method, database) preferred by the `preference` consensus strategy (Q3).
preference_method <- "sintax"
preference_db     <- "EUK_ITS_v2.1"
sam_data_file_name <- "sam_data.csv"
sample_col_name    <- "Sample_names"
fake_ref_fasta     <- "data/data_raw/fake_ref/fake_ref_asv_100.fasta"
taxo_mock_csv      <- "data/data_raw/metadata/taxo_mock.csv"

# Prelude that activates the cutadapt conda env. Override per-machine if your
# conda lives elsewhere.
cutadapt_conda_prelude <-
  "source ~/miniforge3/etc/profile.d/conda.sh && conda activate cutadaptenv && "

# Primer trimming of the references (ROADMAP S6.9) and of the cross-validation
# queries (D1a), with dbpq::cutadapt_rm_primers_db(): ITS1F trimmed when found,
# reverse complement of ITS2 trimmed when found. primer_min_overlap = length
# of ITS2, so a few bases at the end of a record are not taken for a primer.
# The _Fungi_cut databases keep the records without primer sites
# (cut_discard_untrimmed = FALSE; the dbpq default discards them).
primer_min_overlap    <- nchar(rev_primer_sequences)
cut_discard_untrimmed <- FALSE
# Records shorter than cut_min_length after trimming are dropped from the
# _Fungi_cut databases: ITS2-only records starting at the 5.8S site are left
# empty or nearly (317 empty, 1530 < 50 bp in EUK_ITS_v2.1_Fungi_cut; ROADMAP
# S6.9). Same value as the query length filter of cross_val().
cut_min_length        <- 50L

# ITS extraction of the ASVs with ITSx (ROADMAP Q4). ITSx 1.1.3 lives in the
# `itsxenv` conda env (bioconda). Only ITS1 is kept (the ASVs start with 45 bp
# of 18S and end in the 5.8S). ITSx only trims the flanking regions here, it is
# not used as a Fungi filter: "all" uses every organism profile so non-fungal
# ASVs are trimmed too, and undetected ASVs keep their full sequence.
itsx_conda_prelude <-
  "source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxenv && "
itsx_region          <- "ITS1"
itsx_organism_groups <- "all"
# Databases assigned from the ITSx input too (ROADMAP S8.5, lever L7): the
# Fungi-filtered ones (_Fungi, _Fungi_cut); the full releases get the raw ASVs only.
itsx_db_list <- benchmark_dbs$db[benchmark_dbs$simplification != "full"]

# Number of crew workers for parallel assignments. Lower this if your machine
# can't host n_workers * n_threads cores.
n_workers <- 3

# Threads of each assignment computation (ROADMAP S8.3): dada2 runs on its own
# crew worker, sintax / lca / blastn on the n_workers workers. Keep
# assign_threads_dada2 + n_workers * assign_threads_fast <= the cores (6 here).
assign_threads_dada2 <- 3
assign_threads_fast  <- 1

# Threads of each cross-validation target (ROADMAP S8.4), passed as nproc to
# sintax / lca / blastn and as multithread to dada2. The CV pipeline runs
# n_workers targets at once: keep n_workers * cv_threads <= the cores.
cv_threads <- 2

targets_seed <- 22

# When TRUE, the pipelines use the mini_* variants of the reference databases
# for fast smoke-testing. Derived from the {targets} project name
# (_targets.yaml): the *_mini projects write to their own stores, so a smoke
# test can never overwrite a production store. Target names do not change.
#   Sys.setenv(TAR_PROJECT = "assign_taxo_mini"); targets::tar_make()
mini_db <- grepl("_mini$", Sys.getenv("TAR_PROJECT", unset = ""))
message("config.R: mini_db = ", mini_db,
        if (mini_db) " (smoke-test databases mini_*)" else " (full databases)")

cv_fold_number <- 10L
# Publication values (ROADMAP S8.4, 2026-09-14): every fold of a 10-fold CV
# (as Bokulich et al. 2018) on 5000 sequences drawn at random from each
# database. The *_mini projects use smoke-test values: 2 folds, 200 sequences.
cv_fold_tested <- if (mini_db) 2L else cv_fold_number
cv_max_seq     <- if (mini_db) 200L else 5000L
# CV queries are trimmed to the amplicon with the config primers (ROADMAP D1a,
# 2026-09-15); records without the ITS2 site are never queried but stay in the
# training part. cv_oversample records are drawn per wanted query so that every
# database gets cv_max_seq queries (22 % of EUK_ITS_v2.1 records lack the site).
cv_oversample <- 2
