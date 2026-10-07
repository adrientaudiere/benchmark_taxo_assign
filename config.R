# Shared constants for the benchmark pipelines. Sourced by pipelines/dada2.R,
# pipelines/assign_taxo.R, pipelines/cross_val.R, make_databases.R and the
# analysis notebooks.

# ITS-1F (Gardes & Bruns 1993) / ITS2 (White 1990), per Pauvert et al. 2018.
fw_primer_sequences <- "CTTGGTCATTTAGAGGAAGTAA"
rev_primer_sequences <- "GCTGCGTTCTTCATCGATGC"


n_threads <- 4
seq_len_min <- 200
prop_fake <- 0.5

# Minimum BLAST query cover (%) for blastn assignments and cross-validation
# (decision 18). The ASVs start with 30-45 bp of the 18S end that ITS
# references do not contain, so the assign_blastn() default of 95 rejects most
# true hits (reviews/critique_fable/07_results_validity.md, B18).
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
  ~source                , ~provider   , ~release     , ~doi                   , ~url                                                                                   ,
  "Unite_s_all_20250219" , "unite"     , "19.02.2025" , "10.15156/BIO/3301232" , "https://s3.hpc.ut.ee/plutof-public/original/b02db549-5f04-43fc-afb6-02888b594d10.tgz" ,
  "Unite_all_20250219"   , "unite"     , "19.02.2025" , "10.15156/BIO/3301231" , "https://s3.hpc.ut.ee/plutof-public/original/e861a3d6-54f4-42dc-882a-5f129beac39a.tgz" ,
  "EUK_ITS_v2.1"         , "eukaryome" , "2.1"        , NA_character_          , "https://sisu.ut.ee/wp-content/uploads/sites/643/General_EUK_ITS_v2.1.zip"
)

# Databases benchmarked: one source and one simplification each
# (docs/objectives_design.md §1.3, decisions 1, 3, 13, 14, 21, 24, 25).
#   "full"      = the whole release;
#   "Fungi"     = records whose kingdom is exactly Fungi;
#   "Fungi+rep" = "Fungi" plus non-fungal representative records drawn from
#                 the "full" release of the same source (rep_share of the
#                 fungal record count, spread over the retained kingdoms);
#   "Fungi+cut" = "Fungi", then trimmed to the amplicon region with cutadapt
#                 (config primers; records without primer sites kept whole).
#                 Not in the grid since 2026-09-17: see zoom_dbs.
# Decision 7 (revised 2026-09-15): no SSU database. ITS1 is not part of the
# 18S gene, so an SSU record only holds the 18S start of the ASVs.
benchmark_dbs <- tidyr::expand_grid(
  source = c("Unite_s_all_20250219", "Unite_all_20250219", "EUK_ITS_v2.1"),
  simplification = c("full", "Fungi", "Fungi+rep")
)
# Databases built for the final zoom only (trimmed references on the best
# database x method combinations, Pauvert mock only; decisions 24, 25).
# derive_all_variants() builds them, but they are in no assignment grid.
zoom_dbs <- tibble::tribble(
  ~source        , ~simplification ,
  "EUK_ITS_v2.1" , "Fungi+cut"
)
db_suffix <- c(
  full = "",
  Fungi = "_Fungi",
  "Fungi+rep" = "_Fungi_rep",
  "Fungi+cut" = "_Fungi_cut"
)
benchmark_dbs$db <- paste0(
  benchmark_dbs$source,
  db_suffix[benchmark_dbs$simplification]
)
zoom_dbs$db <- paste0(zoom_dbs$source, db_suffix[zoom_dbs$simplification])

# Used by R/values_map.R (assignment grid) and by the Q2 figures (which
# simplification step each database represents).
db_list <- benchmark_dbs$db
# Cross-validation grid: unchanged until its databases are chosen from the mock
# results (objectives_design decision 27), so the CV store keeps its targets.
cv_db_list <- c(
  benchmark_dbs$db[benchmark_dbs$simplification %in% c("full", "Fungi")],
  zoom_dbs$db
)
db_meta <- tibble::tibble(
  db = benchmark_dbs$db,
  db_base = benchmark_dbs$source,
  simplification = benchmark_dbs$simplification
)

# Seed taxonomy of the DADA2 pipeline (the benchmark itself re-assigns).
seed_taxonomy_db <- "Unite_s_all_20250219_Fungi"
# Release whose non-Fungi records feed the external negative controls.
fake_ref_source <- "Unite_s_all_20250219"
# Release the in silico reads of D1b are generated from (objectives_design
# decision 5). Same release as `fake_ref_source` today, but for another reason:
# read from `sources/` so its k__/p__ headers still carry the lineage that is
# the ground truth (analysis/in_silico_simulation.qmd).
insilico_source <- "Unite_s_all_20250219"
# Non-fungal kingdoms retained as groups of representative sequences in the
# `fungi + rep` databases (docs/objectives_design.md §1.3, decisions 13 and
# 21): at least this many records in the release, once the labels matching
# `rep_kingdom_exclude` are removed. The external negative controls are drawn
# only from the kingdoms retained in every release listed here, so that each
# control can be named at Kingdom on every `fungi + rep` database.
rep_kingdom_min_records <- c(
  Unite_all_20250219 = 10L,
  Unite_s_all_20250219 = 10L,
  EUK_ITS_v2.1 = 100L
)
rep_kingdom_exclude <- "^Fungi$|^cf\\.|\\.reg|^_|spike|Incertae_sedis|^Tartumycota$"
# Number of representative records of a `fungi + rep` database, as a share of
# its fungal records (decision 3), rounded up.
rep_share <- 0.05
# A record whose identity with one of the external negative controls reaches
# this value is never drawn as a representative (developer, 2026-09-22): a
# control must not have a close relative in a database whose composition we
# choose. Measured on the first build: 6, 7 and 10 of the 100 controls had a
# representative at >= 97 % identity in the UNITE, UNITE_sh and EUKARYOME
# `fungi + rep` databases. The whole releases keep the controls (they come from
# `fake_ref_source`, so `ext_correct` stays computable there) and the `_Fungi`
# databases hold none.
rep_max_identity_to_external <- 0.97
# Longest record kept in every database of a source, the whole release and
# its `_Fungi` / `_Fungi+rep` variants, representatives included (OD 32,
# 2026-09-24, for the variants; OD 35, 2026-10-01, for the whole release). A
# 32 568 bp Laccaria genomic clone made sintax name every shuffled control on
# `EUK_ITS_v2.1_Fungi` (SA5), and the 32 759 bp plant mitochondrial genomes of
# the whole release drew 55 of 97 to a genus with the patched dada2. Only
# EUKARYOME carries such records: UNITE's longest is 3 173 bp, so a cap there
# would change nothing.
max_record_length <- c(EUK_ITS_v2.1 = 5000L)
# EUKARYOME kingdom names mapped to their UNITE spelling before comparing.
kingdom_aliases <- c(Straminipila = "Stramenopila")
# (method, database) preferred by the `preference` consensus strategy (Q3).
preference_method <- "sintax"
preference_db <- "EUK_ITS_v2.1"
sam_data_file_name <- "sam_data.csv"
sample_col_name <- "Sample_names"

# Sanger sequences of the strains of each mock community, from which every ASV
# and OTU gets its own truth (ROADMAP 0.6, R/mock_truth.R): the table of the
# authors, its separator, the columns holding the strain identifier, the taxon
# name and the sequence, and the cache of GBIF lineages written by
# refresh_gna_lineage_cache(). `dataset` is "pauvert" for the production
# projects, otherwise the key of `bio_datasets`.
mock_sanger_tables <- tibble::tribble(
  ~dataset            , ~path                                                         , ~sep , ~unit_col     , ~name_col        , ~seq_col          , ~lineage_cache                                             ,
  "pauvert"           , "data/data_raw/metadata/pauvert2019_table_s1_mock_sanger.csv" , ";"  , "Mock_Strain" , "Species"        , "Sanger_Sequence" , "data/data_raw/metadata/pauvert_gna_lineage.csv"           ,
  "tedersoo_illumina" , "data/data_raw/metadata/tedersoo_mock_table_s1_sanger.csv"    , ","  , "Specimen"    , "Identification" , "Sanger_read"     , "data/data_raw/metadata/tedersoo_illumina_gna_lineage.csv" ,
  "rolling_mock"      , "data/data_raw/metadata/rolling2022_mock_expected_its1.csv"   , ","  , "Strain"      , "Species"        , "Expected_ITS1"   , "data/data_raw/metadata/rolling_mock_gna_lineage.csv"
)
# rolling_mock has no Sanger reads: its "sequence" is the expected ITS1
# amplicon of each strain published by the authors
# (github.com/thierroll/dada2_custom_fungal, refseq.txt, kept as
# data/data_raw/metadata/rolling2022_refseq_authors.txt). The authors write it
# from the 5.8S end; the csv holds its reverse complement, in the orientation
# of the reads after cutadapt --revcomp (ITS-1-F side first), since
# sanger_hits() searches the plus strand only. Twelve sequences: the two
# C. parapsilosis strains and the two S. cerevisiae ITS1 variants each count
# as a strain.

# ---- Biological-community datasets (D1c) ------------------------------------
# One {targets} project pair per dataset (developer decision 2026-09-16):
# `dada2_<dataset>` builds the phyloseq from fastq, `assign_taxo_<dataset>`
# assigns it. The dataset key is derived from TAR_PROJECT exactly as mini_db
# is, so the production projects are untouched and each dataset owns its store.
#
# Primers VERIFIED 2026-09-16 against `primers.fas` of the study's own
# repository (Mycology-Microbiology-Center/fullITS-multiplatform-eval), whose
# 62 bytes read exactly:
#     >ITS9mun_ITS4ngsUni
#     GTACACACCGCCCGTCG...GCATATHANTAAGSGSAGGcg
# That is cutadapt's linked-adapter form, FORWARD...revcomp(REVERSE), so the
# reverse primer is stored here in its own reading sense (the reverse
# complement of `GCATATHANTAAGSGSAGGcg`), as `rev_primer_sequences` above holds
# the ITS2 of White 1990 in its reading sense. The NextITS defaults quoted in
# its documentation drop the leading G of ITS9mun and the leading CG of
# ITS4ngsUni: do not copy them from there.
#
# File patterns are per-dataset because naming is not uniform. The ENA serves
# these runs as `<run>_1.fastq.gz` / `_2.fastq.gz`, but the local copies were
# renamed to `_R1` / `_R2` to match the mock-community convention, and the
# columns below describe WHAT IS ON DISK, not what the ENA serves. Check the
# folder before changing them: both `cutadapt_remove_primers()` and
# `list_fastq_files()` take these patterns, and a mismatch stops the pipeline
# on its first target with "None file in the folder ... match the pattern_R1".
#
# `fw_primer = NA` means the reads carry no primer site and primer removal is
# skipped (OD 29). Only the in silico set is in that case: InSilicoSeq built
# its reads from reference records that begin inside ITS1, so ITS1F is in 0 of
# the 100 input records and cutadapt with --discard-untrimmed would drop
# everything. `pipelines/dada2_bio.R` branches on it.
bio_datasets <- tibble::tribble(
  ~dataset            , ~raw_dir                                 , ~sam_data                        , ~fw_primer               , ~rev_primer              , ~pattern_r1   , ~pattern_r2   , ~name_strip        , ~itsx_region     , ~otu_input , ~min_occur  , ~cutadapt_args ,
  "tedersoo_illumina" , "data/data_raw/rawseq_tedersoo_illumina" , "sam_data_tedersoo_illumina.csv" , "GTACACACCGCCCGTCG"      , "CGCCTSCSCTTANTDATATGC"  , "_R1"         , "_R2"         , "_R[12]\\.fastq.*" , "ITS1,ITS2,full" , TRUE       , NA_integer_ , ""             ,
  "insilico"          , "data/data_raw/insilico"                 , "sam_data_insilico.csv"          , NA_character_            , NA_character_            , "_R1"         , "_R2"         , "_R[12]\\.fastq.*" , "ITS1,ITS2,full" , TRUE       , NA_integer_ , ""             ,
  "gssp_air"          , "data/data_raw/bio_gssp_air"             , NA_character_                    , NA_character_            , NA_character_            , NA_character_ , NA_character_ , NA_character_      , "ITS2"           , FALSE      , 10L         , ""             ,
  "korhonen_wood"     , "data/data_raw/bio_korhonen_wood"        , NA_character_                    , NA_character_            , NA_character_            , NA_character_ , NA_character_ , NA_character_      , "ITS2"           , FALSE      ,  1L         , ""             ,
  "rolling_mock"      , "data/data_raw/rawseq_rolling_mock"      , "sam_data_rolling_mock.csv"      , "CTTGGTCATTTAGAGGAAGTAA" , "GTTCAAAGAYTCGATGATTCAC" , "_R1"         , "_R2"         , "_R[12]\\.fastq.*" , "ITS1"           , TRUE       , NA_integer_ , "--revcomp"
)
# rolling_mock (OD 42): the fungal ITS1 mock of Rolling et al. 2022 (JCI
# Insight 7:e151663; SRA PRJNA779405), 11 strains of 9 medically relevant
# species, MiSeq 2 x 300. Three runs with the standard ITS-1-F, used as
# replicates of one community (presence only): balanced SRR16904671, extreme1
# SRR16904670, extreme2 SRR16904669; the wobble-primer and S. cerevisiae
# enriched runs are left out. Primers ITS-1-F and 5.8S-1R as printed in the
# paper's Methods, each in its own reading sense; the library mixes both
# orientations, hence "--revcomp".
# cutadapt_args: extra arguments of the primer removal
# (MiscMetabar::cutadapt_remove_primers(cutadapt_args = ), pipelines/dada2_bio.R),
# "" for none. "--revcomp" (cutadapt >= 4.6) for a library whose R1 holds
# either primer: cutadapt also tries each pair with R1 and R2 exchanged and
# writes the matching pairs in the forward orientation. Measured on the
# Rolling et al. 2022 mock (SRR16904671, analysis/sandbox/check_rolling_orientation.R):
# half of the pairs are reversed, 49.9 % are kept without it and 98.7 % with
# it, the same pairs and reads as two passes with the primers swapped.
# The columns `taxo_mock` (a set-membership truth, replaced by the per-unit
# truth of R/mock_truth.R) and `cut_suffix` (amplicon-specific trimmed
# references, NA in every row: out of the grid, decisions 24, 25) were
# removed on 2026-10-07 (code review §5); they are in git history.
# itsx_region: ITSx regions of the dataset, comma-separated (ROADMAP 0.5, OD
# 34): "ITS1", "ITS2" or "full" (ITS1 + 5.8S + ITS2); each one is an input of
# the grid. A dataset whose amplicon covers the whole ITS gets all three.

# Datasets that publish their OTUs with representative sequences skip dada2
# (OD 30): their `dada2_<dataset>` project runs pipelines/load_bio.R, which
# builds `d_asv` from the published tables (R/bio_loaders.R) instead of from
# fastq, so the fastq columns above are NA. min_occur: an OTU is kept when
# it occurs in at least that many samples (OD 30b, revised 2026-09-26: the
# common OTUs, not a random draw); NA for the fastq datasets.

# otu_input: FALSE for the datasets whose units are already published OTUs
# (OD 30c): re-clustering them at 97 % answers nothing, so assign_taxo leaves
# out its OTU input (d_vs, the otu_* rows of values_map) and the ASV clusters
# of Q5. TRUE keeps the three inputs, as in the production project.

# Dataset of the running project: the TAR_PROJECT name without its
# `dada2_` / `assign_taxo_` prefix ("assign_taxo", "assign_taxo_mini", a dataset
# name...); `bio_cfg` is NULL unless it names a row of bio_datasets.
bio_dataset <- sub(
  "^(dada2|assign_taxo)_",
  "",
  Sys.getenv("TAR_PROJECT", unset = "")
)
bio_cfg <- if (bio_dataset %in% bio_datasets$dataset) {
  as.list(bio_datasets[match(bio_dataset, bio_datasets$dataset), ])
} else {
  NULL
}
if (!is.null(bio_cfg)) {
  message("config.R: biological dataset = ", bio_dataset)
}

# DADA2 store that assign_taxo reads its `d_asv` from. The production project
# keeps store_dada2; a dataset project reads its own store.
dada2_store <- if (is.null(bio_cfg)) {
  "store_dada2"
} else {
  paste0("store_dada2_", bio_dataset)
}
fake_ref_fasta <- "data/data_raw/fake_ref/fake_ref_asv_100.fasta"
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
primer_min_overlap <- nchar(rev_primer_sequences)
cut_discard_untrimmed <- FALSE
# Records shorter than cut_min_length after trimming are dropped from the
# _Fungi_cut databases: ITS2-only records starting at the 5.8S site are left
# empty or nearly (317 empty, 1530 < 50 bp in EUK_ITS_v2.1_Fungi_cut; ROADMAP
# S6.9). Same value as the query length filter of cross_val().
cut_min_length <- 50L

# ITS extraction of the ASVs with ITSx (ROADMAP Q4). ITSx 1.1.3 lives in the
# `itsxenv` conda env (bioconda). Only ITS1 is kept (the ASVs start with 45 bp
# of 18S and end in the 5.8S). ITSx only trims the flanking regions here, it is
# not used as a Fungi filter: "all" uses every organism profile so non-fungal
# ASVs are trimmed too, and undetected ASVs keep their full sequence.
itsx_conda_prelude <-
  "source ~/miniforge3/etc/profile.d/conda.sh && conda activate itsxenv && "
# Region kept by ITSx, per dataset (ROADMAP 0.5): ITS1 for the ITS1F / ITS2
# amplicon of the Pauvert mock, the region of `bio_datasets$itsx_region` for a
# dataset project (the full ITS for Tedersoo, whose ASVs carry 154 bp of 18S
# and 39 bp of 28S around it; docs/experimental_design.md §5.1).
itsx_region <- if (is.null(bio_cfg)) "ITS1" else bio_cfg$itsx_region
# ITSx regions assigned for the running dataset (OD 34): one input per region,
# e.g. "ITS1,ITS2,full" for the datasets whose amplicon covers the whole ITS.
itsx_regions <- strsplit(itsx_region, ",", fixed = TRUE)[[1]]
# Whether assign_taxo runs its OTU input (`bio_datasets$otu_input`, OD 30c).
run_otu_input <- if (is.null(bio_cfg)) TRUE else bio_cfg$otu_input
# Units of each ITSx region (OD 34): the ASVs, and their 97 % clusters when
# the dataset has an OTU input.
itsx_units <- if (run_otu_input) c("asv", "otu") else "asv"
itsx_organism_groups <- "all"

# Number of crew workers for parallel assignments. Lower this if your machine
# can't host n_workers * n_threads cores.
# Lowered from 3 to 2 on 2026-09-22 after the in silico run exhausted the
# machine (0 GB available, 21 GB of swap, no target completed in 45 s) while a
# dada2 worker held 46.8 GB beside three fast workers. The pipeline is
# unchanged; what changed is the headroom, since the desktop now holds more
# than it did during the mock run of the same afternoon, which completed at a
# 58.1 GB peak. The fast workers are the only lever that costs nothing: they
# are vsearch and BLAST+ at 1.5-6 GB each, while the single dada2 worker is
# irreducible. See ROADMAP S8.6.
n_workers <- 2
# One computation at a time (developer, 2026-09-29): when TRUE, assign_taxo
# runs every computation, dada2 included, on a single crew worker that is
# relaunched after each task, and n_workers is ignored. The run of 2026-09-28
# 22:01 was killed by the kernel OOM killer after 24 min: dada2 x full
# EUK_ITS_v2.1 (56.4 GB peak) ran beside two sintax computations on the same
# database, beyond the 62 GB of RAM and 19 GB of swap. Also the concurrency of
# the final timing run (OD 31, ROADMAP M2.4).
single_worker <- TRUE

# Threads of each assignment computation, per method (M2.4, OD 31 revised by
# OD 34: one computation at a time, single_worker = TRUE). dada2, lca and
# blastn get every core of the machine, which is an input of the problem;
# sintax stays at 1 thread, the only setting at which its seeded draw is
# reproducible (sintax_seed below). Until 2026-10-07: 3 threads for dada2 and
# 1 for the others, sized for dada2 running beside n_workers fast workers.
assign_threads <- c(
  dada2 = parallel::detectCores(),
  sintax = 1L,
  lca = parallel::detectCores(),
  blastn = parallel::detectCores()
)
# lca search settings (developer, 2026-09-26, OD 33): the common LCA of QIIME 2's
# classify-consensus-vsearch, the "balanced" configuration Bokulich et al. 2018
# optimised on 16S and fungal ITS, instead of MiscMetabar's stampa defaults
# (top_hits_only = TRUE, id = 0.5, maxaccepts = 0). stampa was near-exhaustive
# on real ITS queries (4.5-4.75 s each against 102 137 records, 25-40 times
# slower) and named 62 % of the non-fungal controls; its results are kept in
# HISTORY.md. The consensus threshold is swept by values_map (lca_cutoff).
lca_top_hits_only <- FALSE
lca_id <- 0.8
lca_maxaccepts <- 10
# Pinned rather than left to the default of MiscMetabar::assign_vsearch_lca()
# (32, the value OD 33 states), so a changed default in the checkout cannot
# silently change the design (code review §3.5).
lca_maxrejects <- 32L
# The settings of the blastn and sintax searches, pinned for the same reason
# (code review §3.5, developer 2026-10-07): the values are the defaults of
# BLAST+ 2.17.0 megablast and of vsearch --sintax (MiscMetabar passes
# --sintax_random), written out so that another BLAST+ or MiscMetabar cannot
# change the search. blast_pq() adds -num_threads and -outfmt itself.
blastn_search_args <- paste(
  "-task megablast -word_size 28 -evalue 10",
  "-reward 1 -penalty -2 -gapopen 0 -gapextend 0 -dust '20 64 1'"
)
blastn_args <- paste(blastn_search_args, "-max_target_seqs 500")
sintax_cmd_args <- "--sintax_random --strand plus --wordlength 8"

# Threads of each cross-validation target (ROADMAP S8.4), passed as nproc to
# sintax / lca / blastn and as multithread to dada2. The CV pipeline runs
# cv_n_workers_fast targets at once: keep cv_n_workers_fast * cv_threads
# within the cores.
cv_threads <- 2

# Fast workers of the cross-validation pipeline only (decision 2026-09-16).
# Kept separate from n_workers, which assign_taxo shares, because the CV is
# memory-bound and assign_taxo is not: since cv_reduce_reference is FALSE, a CV
# target holds a whole EUKARYOME reference. Peaks measured on EUK_ITS_v2.1:
# dada2 >= 45 GB (run killed before the end), lca 6.2, sintax 5.9, blastn 3.5.
# The single dada2 worker is resident almost continuously, so 3 fast workers
# would reach 60.6 GB on a 62 GB machine, against 57.1 GB for 2.
# See docs/diagnostics_cv_and_databases.md section 7.
cv_n_workers_fast <- 2

targets_seed <- 22
# Seed of the vsearch random number generator for sintax (ROADMAP SA6):
# --sintax_random draws the bootstrap k-mers and breaks ties at random, so two
# unseeded runs against EUK_ITS_v2.1_Fungi gave a different genus to 11 of 71
# in silico ASVs. Reproducible with a single thread, which
# assign_threads[["sintax"]] = 1 guarantees.
sintax_seed <- targets_seed

# When TRUE, the pipelines use the mini_* variants of the reference databases
# for fast smoke-testing. Derived from the {targets} project name
# (_targets.yaml): the *_mini projects write to their own stores, so a smoke
# test can never overwrite a production store. Target names do not change.
#   Sys.setenv(TAR_PROJECT = "assign_taxo_mini"); targets::tar_make()
mini_db <- grepl("_mini$", Sys.getenv("TAR_PROJECT", unset = ""))
message(
  "config.R: mini_db = ",
  mini_db,
  if (mini_db) " (smoke-test databases mini_*)" else " (full databases)"
)

cv_fold_number <- 3L
# Publication values (ROADMAP S8.4, revised 2026-09-16): every fold of a 3-fold
# CV on 5000 queries drawn at random from each database. Three folds, not ten:
# since cv_reduce_reference is FALSE the reference is the whole database and a
# fold only removes its own queries (5000 of 156 820 on Unite_all_20250219), so
# the extra folds would only rebuild the index for nothing (docs/
# diagnostics_cv_and_databases.md section 7). The *_mini projects test 2 folds
# of 200 sequences.
cv_fold_tested <- if (mini_db) 2L else cv_fold_number
cv_max_seq <- if (mini_db) 200L else 5000L
# Reference of a CV target (ROADMAP B24, decision 2026-09-16): FALSE keeps the
# whole database and removes only the tested queries, as Bokulich et al. 2018,
# instead of reducing the reference to the drawn pool. The queries are drawn
# exactly as before; only the reference changes. TRUE restores the behaviour of
# the runs made before 2026-09-16, whose results measured a reference of about
# 5 400 records whatever the database.
cv_reduce_reference <- FALSE
# CV queries are trimmed to the amplicon with the config primers (ROADMAP D1a,
# 2026-09-15); records without the ITS2 site are never queried but stay in the
# training part. cv_oversample records are drawn per wanted query so that every
# database gets cv_max_seq queries (22 % of EUK_ITS_v2.1 records lack the site).
cv_oversample <- 2
