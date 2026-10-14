# Settings of the cross-validation (OD 28: parked), moved out of config.R on
# 2026-10-14. Source it after config.R: it reads benchmark_dbs, zoom_dbs and
# mini_db. Like every file of Archives/cross_val/, it keeps its
# project-relative place: restored, it sits at the project root.

# Cross-validation grid: unchanged until its databases are chosen from the mock
# results (objectives_design decision 27), so the CV store keeps its targets.
cv_db_list <- c(
  benchmark_dbs$db[benchmark_dbs$simplification %in% c("full", "Fungi")],
  zoom_dbs$db
)

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
