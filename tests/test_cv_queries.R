# Trimmed cross-validation queries of R/cv_queries.R. The subsample rule is
# tested without cutadapt; the cutadapt run is skipped when the cutadaptenv
# conda env is missing. Run from the project root:
#   Rscript tests/test_cv_queries.R

library("testthat")
library("here")

here::i_am("tests/test_cv_queries.R")
source(here("config.R"))
source(here("R/load_pqverse.R"))
source(here("R/cv_queries.R"))

test_that("record_ids returns the sintax identifiers and checks the paired lengths", {
  path <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">A_sp|KF1|SH1.10FU|refs;tax=k:Fungi,g:X",
      "ACGTACGT",
      ">EUK2;tax=k:Fungi,g:Y",
      "ACGTA"
    ),
    path
  )
  expect_equal(record_ids(path), c("A_sp|KF1|SH1.10FU|refs", "EUK2"))
  expect_equal(
    record_ids(path, widths = c(8L, 5L)),
    c("A_sp|KF1|SH1.10FU|refs", "EUK2")
  )
  expect_error(record_ids(path, widths = c(5L, 8L)), "same records")
})

test_that("cv_select_queries keeps the pool up to the max_seq-th query", {
  pool <- paste0("r", 1:8)
  sel <- cv_select_queries(pool, c("r2", "r3", "r5", "r6", "r8"), max_seq = 3)
  expect_equal(sel$pool, paste0("r", 1:5))
  expect_equal(sel$queries, c("r2", "r3", "r5"))
  expect_equal(cv_select_queries(pool, pool)$queries, pool)
  expect_message(
    few <- cv_select_queries(pool, c("r4", "r1"), max_seq = 3),
    "2 queries only"
  )
  expect_equal(few$pool, pool)
  expect_equal(few$queries, c("r1", "r4"))
})

test_that("trim_cv_queries cuts to the amplicon and drops records without the ITS2 site", {
  cutadapt_ok <- system2(
    "bash",
    c("-c", shQuote(paste0(cutadapt_conda_prelude, "cutadapt --version"))),
    stdout = FALSE,
    stderr = FALSE
  ) ==
    0
  skip_if_not(cutadapt_ok, "cutadapt (cutadaptenv conda env) not available")
  load_pqverse("dbpq")
  its2_rc <- as.character(
    Biostrings::reverseComplement(Biostrings::DNAString(rev_primer_sequences))
  )
  insert <- paste(rep("ACGTTGACCATGGTACCAGT", 5), collapse = "")
  dna <- Biostrings::DNAStringSet(c(
    "id1;tax=k:Fungi,g:A" = paste0(
      "TTTTT",
      fw_primer_sequences,
      insert,
      its2_rc,
      "GGGGG"
    ),
    "id2;tax=k:Fungi,g:B" = paste0(insert, its2_rc, "GGGGG"),
    "id3;tax=k:Fungi,g:C" = paste0(insert, substr(its2_rc, 1, 8)),
    "id4;tax=k:Fungi,g:D" = insert
  ))
  res <- trim_cv_queries(
    dna,
    primer_fw = fw_primer_sequences,
    primer_rev = rev_primer_sequences,
    min_overlap = primer_min_overlap,
    prelude = cutadapt_conda_prelude
  )
  expect_equal(names(res), c("id1;tax=k:Fungi,g:A", "id2;tax=k:Fungi,g:B"))
  expect_true(all(as.character(res) == insert))
})
