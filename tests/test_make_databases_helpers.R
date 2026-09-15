# The shell-based helpers of make_databases.R (record-wise filtering, head,
# balanced fake-reference subset) on a tiny fasta. Run from the project root:
#   Rscript tests/test_make_databases_helpers.R

library("testthat")
library("here")

here::i_am("tests/test_make_databases_helpers.R")
source(here("R/load_pqverse.R"))
load_pqverse("dbpq")
source(here("make_databases.R"))

# Four records, the third one wrapped on two lines (multi-line fasta).
tiny_fasta <- function() {
  path <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">r1|k__Fungi;p__Ascomycota;c__A;o__O;f__F;g__G;s__S1",
    "ACGTACGTAC",
    ">r2|k__Metazoa;p__Chordata;c__B;o__O;f__F;g__G;s__S2",
    "TTGACCGGTT",
    ">r3|k__Fungi;p__Basidiomycota;c__C;o__O;f__F;g__G;s__S3",
    "GGGCC",
    "CAAAT",
    ">r4|k__Metazoa;p__Chordata;c__D;o__O;f__F;g__G;s__S4",
    "CCCGGGAAAT"
  ), path)
  path
}

read_records <- function(path) {
  seqs <- Biostrings::readDNAStringSet(path)
  setNames(as.character(seqs), names(seqs))
}

test_that("derive_no_pattern drops matching records and joins multi-line sequences", {
  out <- tempfile(fileext = ".fasta")
  derive_no_pattern(tiny_fasta(), out, pattern = "Fungi")
  recs <- read_records(out)
  expect_equal(length(recs), 2)
  expect_false(any(grepl("Fungi", names(recs))))
  expect_equal(unname(recs[grepl("^r2", names(recs))]), "TTGACCGGTT")
})

test_that("derive_mini keeps the first n records, not the first n lines", {
  out <- tempfile(fileext = ".fasta")
  derive_mini(tiny_fasta(), out, n = 3)
  recs <- read_records(out)
  expect_equal(length(recs), 3)
  expect_equal(unname(recs[grepl("^r3", names(recs))]), "GGGCCCAAAT")
})

test_that("skip_if_exists is idempotent unless force = TRUE", {
  out <- tempfile(fileext = ".fasta")
  derive_mini(tiny_fasta(), out, n = 1)
  expect_message(derive_mini(tiny_fasta(), out, n = 4), "already exists")
  expect_equal(length(read_records(out)), 1)
  derive_mini(tiny_fasta(), out, n = 4, force = TRUE)
  expect_equal(length(read_records(out)), 4)
})

test_that("derive_sintax keeps the kingdom of UNITE general-release headers", {
  input <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">Fusarium_sp|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_sp",
    "ACGTACGTAC",
    ">Agaricus_sp|KF2|SH2.10FU|refs|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_sp",
    "TTGACCGGTT"
  ), input)
  out <- tempfile(fileext = ".fasta")
  derive_sintax(input, out)
  headers <- grep("^>", readLines(out), value = TRUE)
  expect_length(headers, 2)
  expect_true(all(grepl(";tax=k:Fungi,p:", headers, fixed = TRUE)))
  expect_false(any(grepl("|k__", headers, fixed = TRUE)))
})

test_that("derive_dada2 keeps the kingdom of UNITE general-release headers", {
  input <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">Fusarium_sp|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_sp",
    "ACGTACGTAC"
  ), input)
  out <- tempfile(fileext = ".fasta")
  derive_dada2(input, out)
  expect_equal(
    readLines(out, n = 1),
    ">Fungi;Ascomycota;Sordariomycetes;Hypocreales;Nectriaceae;Fusarium;Fusarium_sp;"
  )
})

test_that("derive_sintax and derive_dada2 keep only the name of EUKARYOME qualified ranks", {
  input <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">EUK1;k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Russulales;f__Russulaceae;g__Lactarius(Fungi);s__rufus",
    "ACGT",
    ">EUK2;k__Fungi;p__Ascomycota;c__Pichiomycetes;o__Pichiales;f__Pichiaceae;g__(Candida];s__californica",
    "ACGT",
    ">EUK3;k__Fungi;p__Mortierellomycota;c__Mortierellomycetes;o__Mortierellales;f__Mortierellaceae;g__Mortierella.s.str;s__unclassified",
    "ACGT",
    ">EUK4;k__Fungi;p__Ascomycota;c__Pezizomycetes;o__Pezizales;f__Pezizaceae.nom.prov;g__Peziza.s.str.;s__badia",
    "ACGT",
    ">EUK5;k__Alveolata;p__Ciliophora;c__Spirotrichea;o__Sporadotrichida;f__Gonostomatidae(Sporadotrichida);g__(Candida);s__argentea",
    "ACGT",
    ">EUK6;k__Fungi;p__Ascomycota;c__Dothideomycetes;o__Dothideales;f__Dothideales.fam.incertae.sedis;g__Dothideales.gen05;s__unclassified",
    "ACGT"
  ), input)
  expected_genus <- c("Lactarius", "Candida", "Mortierella", "Peziza", "Candida", "Dothideales.gen05")
  expected_family <- c(
    "Russulaceae", "Pichiaceae", "Mortierellaceae", "Pezizaceae",
    "Gonostomatidae", "Dothideales.fam.incertae.sedis"
  )

  out_dir <- tempfile("derive_")
  sintax <- file.path(out_dir, "sintax.fasta")
  derive_sintax(input, sintax)
  sintax_headers <- grep("^>", readLines(sintax), value = TRUE)
  expect_equal(sub(".*,g:([^,]*),.*", "\\1", sintax_headers), expected_genus)
  expect_equal(sub(".*,f:([^,]*),.*", "\\1", sintax_headers), expected_family)

  dada2 <- file.path(out_dir, "dada2.fasta")
  derive_dada2(input, dada2)
  dada2_ranks <- strsplit(sub("^>", "", grep("^>", readLines(dada2), value = TRUE)), ";")
  expect_equal(vapply(dada2_ranks, \(x) x[6], character(1)), expected_genus)
  expect_equal(vapply(dada2_ranks, \(x) x[5], character(1)), expected_family)
  expect_setequal(list.files(out_dir), c("sintax.fasta", "dada2.fasta"))
})

test_that("derive_kingdom_only keeps kingdom Fungi only, in sintax and dada2 formats", {
  sintax_in <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">a;tax=k:Fungi,p:Ascomycota,f:Nectriaceae", "ACGT",
    ">b;tax=k:Metazoa,p:Cnidaria,f:Fungiidae", "ACGT",
    ">c;tax=k:cf.Fungi,p:unclassified", "ACGT"
  ), sintax_in)
  sintax_out <- tempfile(fileext = ".fasta")
  derive_kingdom_only(sintax_in, sintax_out, format = "sintax")
  expect_equal(grep("^>", readLines(sintax_out), value = TRUE),
               ">a;tax=k:Fungi,p:Ascomycota,f:Nectriaceae")

  dada2_in <- tempfile(fileext = ".fasta")
  writeLines(c(
    ">Fungi;Ascomycota;Nectriaceae;", "ACGT",
    ">Metazoa;Cnidaria;Fungiidae;", "ACGT"
  ), dada2_in)
  dada2_out <- tempfile(fileext = ".fasta")
  derive_kingdom_only(dada2_in, dada2_out, format = "dada2")
  expect_equal(grep("^>", readLines(dada2_out), value = TRUE),
               ">Fungi;Ascomycota;Nectriaceae;")
})

test_that("download_reference_source stores the general FASTA and its provenance", {
  work <- tempfile("sources_")
  dir.create(work)
  testthat::local_mocked_bindings(
    download_unite_db = function(dest_dir, url, extract, ...) {
      path <- file.path(dest_dir, "sh_general_release_dynamic_s_all_19.02.2025.fasta")
      writeLines(c(">a|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota", "ACGT"), path)
      path
    },
    .package = "dbpq"
  )
  src <- tibble::tibble(
    source = "Unite_test", provider = "unite", release = "19.02.2025",
    doi = "10.15156/BIO/0", url = "https://example.org/x.tgz"
  )
  out <- download_reference_source(src, dir = work)
  expect_equal(basename(out), "Unite_test.fasta")
  expect_false(dir.exists(file.path(work, "Unite_test_download")))
  prov <- utils::read.csv(file.path(work, "Unite_test.provenance.csv"))
  expect_equal(prov$original_file, "sh_general_release_dynamic_s_all_19.02.2025.fasta")
  expect_equal(prov$md5, unname(tools::md5sum(out)))
  expect_message(download_reference_source(src, dir = work), "already exists")
})

test_that("extract_single_fasta handles a zip holding a 7z archive (EUKARYOME v2.1)", {
  skip_if(Sys.which("7z") == "" || Sys.which("zip") == "", "7z or zip not available")
  src <- tempfile("euk_src_")
  dir.create(src)
  fasta <- file.path(src, "General_EUK_TEST_v2.1.fasta")
  writeLines(c(">EUK1;k__Fungi;p__Ascomycota", "ACGT"), fasta)
  seven <- file.path(src, "General_EUK_TEST_v2.1.7z")
  system2("7z", c("a", "-y", shQuote(seven), shQuote(fasta)), stdout = FALSE)
  zip_file <- file.path(src, "General_EUK_TEST_v2.1.zip")
  utils::zip(zip_file, files = seven, flags = "-j -q")

  exdir <- tempfile("euk_extract_")
  dir.create(exdir)
  out <- extract_single_fasta(zip_file, exdir)
  expect_equal(basename(out), "General_EUK_TEST_v2.1.fasta")
  expect_equal(readLines(out, n = 1), ">EUK1;k__Fungi;p__Ascomycota")
  expect_false(file.exists(file.path(exdir, "General_EUK_TEST_v2.1.7z")))
})

test_that("derive_fake_ref takes one record per phylum first, then fills up", {
  out <- tempfile(fileext = ".fasta")
  derive_fake_ref(input = tiny_fasta(), output = out, n = 3, seed = 1)
  recs <- read_records(out)
  expect_equal(length(recs), 3)
  phyla <- stringr::str_match(names(recs), "p__\\s*(.*?)\\s*;c__")[, 2]
  expect_equal(length(unique(phyla)), 3)

  out_all <- tempfile(fileext = ".fasta")
  derive_fake_ref(input = tiny_fasta(), output = out_all, n = 10, seed = 1)
  expect_equal(length(read_records(out_all)), 4)
})
