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
  writeLines(
    c(
      ">r1|k__Fungi;p__Ascomycota;c__A;o__O;f__F;g__G;s__S1",
      "ACGTACGTAC",
      ">r2|k__Metazoa;p__Chordata;c__B;o__O;f__F;g__G;s__S2",
      "TTGACCGGTT",
      ">r3|k__Fungi;p__Basidiomycota;c__C;o__O;f__F;g__G;s__S3",
      "GGGCC",
      "CAAAT",
      ">r4|k__Metazoa;p__Chordata;c__D;o__O;f__F;g__G;s__S4",
      "CCCGGGAAAT"
    ),
    path
  )
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
  input <- tiny_fasta()
  out <- tempfile(fileext = ".fasta")
  derive_mini(input, out, n = 1)
  expect_message(derive_mini(input, out, n = 4), "already exists")
  expect_equal(length(read_records(out)), 1)
  derive_mini(input, out, n = 4, force = TRUE)
  expect_equal(length(read_records(out)), 4)
})

test_that("skip_if_exists reports an output older than its input, keeps it", {
  input <- tiny_fasta()
  out <- tempfile(fileext = ".fasta")
  derive_mini(input, out, n = 1)
  Sys.setFileTime(out, Sys.time() - 3600)
  expect_warning(derive_mini(input, out, n = 4), "STALE")
  expect_equal(length(read_records(out)), 1)
})

test_that("atomic_output leaves the old file when the writer fails", {
  out <- tempfile(fileext = ".fasta")
  writeLines(c(">a", "ACGT"), out)
  expect_error(atomic_output(out, function(part) {
    writeLines(">b", part)
    stop("interrupted")
  }))
  expect_equal(readLines(out), c(">a", "ACGT"))
  expect_false(file.exists(sub("\\.fasta$", ".part.fasta", out)))
  atomic_output(out, \(part) writeLines(c(">b", "TTTT"), part))
  expect_equal(readLines(out), c(">b", "TTTT"))
})

test_that("derive_sintax keeps the kingdom of UNITE general-release headers", {
  input <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">Fusarium_sp|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_sp",
      "ACGTACGTAC",
      ">Agaricus_sp|KF2|SH2.10FU|refs|k__Fungi;p__Basidiomycota;c__Agaricomycetes;o__Agaricales;f__Agaricaceae;g__Agaricus;s__Agaricus_sp",
      "TTGACCGGTT"
    ),
    input
  )
  out <- tempfile(fileext = ".fasta")
  derive_sintax(input, out)
  headers <- grep("^>", readLines(out), value = TRUE)
  expect_length(headers, 2)
  expect_true(all(grepl(";tax=k:Fungi,p:", headers, fixed = TRUE)))
  expect_false(any(grepl("|k__", headers, fixed = TRUE)))
})

test_that("derive_dada2 keeps the kingdom of UNITE general-release headers", {
  input <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">Fusarium_sp|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_sp",
      "ACGTACGTAC"
    ),
    input
  )
  out <- tempfile(fileext = ".fasta")
  derive_dada2(input, out)
  expect_equal(
    readLines(out, n = 1),
    ">Fungi;Ascomycota;Sordariomycetes;Hypocreales;Nectriaceae;Fusarium;Fusarium_sp;"
  )
})

test_that("derive_sintax and derive_dada2 keep only the name of EUKARYOME qualified ranks", {
  input <- tempfile(fileext = ".fasta")
  writeLines(
    c(
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
    ),
    input
  )
  expected_genus <- c(
    "Lactarius",
    "Candida",
    "Mortierella",
    "Peziza",
    "Candida",
    "Dothideales.gen05"
  )
  expected_family <- c(
    "Russulaceae",
    "Pichiaceae",
    "Mortierellaceae",
    "Pezizaceae",
    "Gonostomatidae",
    "Dothideales.fam.incertae.sedis"
  )

  out_dir <- tempfile("derive_")
  sintax <- file.path(out_dir, "sintax.fasta")
  derive_sintax(input, sintax)
  sintax_headers <- grep("^>", readLines(sintax), value = TRUE)
  expect_equal(sub(".*,g:([^,]*),.*", "\\1", sintax_headers), expected_genus)
  expect_equal(sub(".*,f:([^,]*),.*", "\\1", sintax_headers), expected_family)

  dada2 <- file.path(out_dir, "dada2.fasta")
  derive_dada2(input, dada2)
  dada2_ranks <- strsplit(
    sub("^>", "", grep("^>", readLines(dada2), value = TRUE)),
    ";"
  )
  expect_equal(vapply(dada2_ranks, \(x) x[6], character(1)), expected_genus)
  expect_equal(vapply(dada2_ranks, \(x) x[5], character(1)), expected_family)
  expect_setequal(list.files(out_dir), c("sintax.fasta", "dada2.fasta"))
})

test_that("derive_kingdom_only keeps kingdom Fungi only, in sintax and dada2 formats", {
  sintax_in <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">a;tax=k:Fungi,p:Ascomycota,f:Nectriaceae",
      "ACGT",
      ">b;tax=k:Metazoa,p:Cnidaria,f:Fungiidae",
      "ACGT",
      ">c;tax=k:cf.Fungi,p:unclassified",
      "ACGT"
    ),
    sintax_in
  )
  sintax_out <- tempfile(fileext = ".fasta")
  derive_kingdom_only(sintax_in, sintax_out, format = "sintax")
  expect_equal(
    grep("^>", readLines(sintax_out), value = TRUE),
    ">a;tax=k:Fungi,p:Ascomycota,f:Nectriaceae"
  )

  dada2_in <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">Fungi;Ascomycota;Nectriaceae;",
      "ACGT",
      ">Metazoa;Cnidaria;Fungiidae;",
      "ACGT"
    ),
    dada2_in
  )
  dada2_out <- tempfile(fileext = ".fasta")
  derive_kingdom_only(dada2_in, dada2_out, format = "dada2")
  expect_equal(
    grep("^>", readLines(dada2_out), value = TRUE),
    ">Fungi;Ascomycota;Nectriaceae;"
  )
})

test_that("download_reference_source stores the general FASTA and its provenance", {
  work <- tempfile("sources_")
  dir.create(work)
  testthat::local_mocked_bindings(
    download_unite_db = function(dest_dir, url, extract, ...) {
      path <- file.path(
        dest_dir,
        "sh_general_release_dynamic_s_all_19.02.2025.fasta"
      )
      writeLines(c(">a|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota", "ACGT"), path)
      path
    },
    .package = "dbpq"
  )
  src <- tibble::tibble(
    source = "Unite_test",
    provider = "unite",
    release = "19.02.2025",
    doi = "10.15156/BIO/0",
    url = "https://example.org/x.tgz"
  )
  out <- download_reference_source(src, dir = work)
  expect_equal(basename(out), "Unite_test.fasta")
  expect_false(dir.exists(file.path(work, "Unite_test_download")))
  prov <- utils::read.csv(file.path(work, "Unite_test.provenance.csv"))
  expect_equal(
    prov$original_file,
    "sh_general_release_dynamic_s_all_19.02.2025.fasta"
  )
  expect_equal(prov$md5, unname(tools::md5sum(out)))
  expect_message(download_reference_source(src, dir = work), "already exists")
})

test_that("extract_single_fasta handles a zip holding a 7z archive (EUKARYOME v2.1)", {
  skip_if(
    Sys.which("7z") == "" || Sys.which("zip") == "",
    "7z or zip not available"
  )
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
  # Any seed: the single-record phylum (Basidiomycota, record 3) must be
  # drawn, and the file must hold min(n, records) records (code review §1.2).
  input <- tiny_fasta()
  for (seed in 1:20) {
    out <- tempfile(fileext = ".fasta")
    derive_fake_ref(input = input, output = out, n = 3, seed = seed)
    recs <- read_records(out)
    phyla <- stringr::str_match(names(recs), "p__\\s*(.*?)\\s*;c__")[, 2]
    expect_equal(length(recs), 3, info = paste("seed", seed))
    expect_setequal(
      phyla,
      c("Ascomycota", "Chordata", "Basidiomycota")
    )
  }

  out_all <- tempfile(fileext = ".fasta")
  derive_fake_ref(input = input, output = out_all, n = 10, seed = 1)
  expect_equal(length(read_records(out_all)), 4)
})

test_that("derive_fake_ref draws only from `kingdoms` when given", {
  out <- tempfile(fileext = ".fasta")
  derive_fake_ref(
    input = tiny_fasta(),
    output = out,
    n = 10,
    seed = 1,
    kingdoms = "Metazoa"
  )
  expect_setequal(sub("\\|.*", "", names(read_records(out))), c("r2", "r4"))
  expect_error(
    derive_fake_ref(
      input = tiny_fasta(),
      output = tempfile(fileext = ".fasta"),
      kingdoms = "Rhizaria"
    ),
    "no record"
  )
})

test_that("kingdom_counts reads sintax headers and applies aliases", {
  path <- tempfile(fileext = ".fasta")
  writeLines(
    c(
      ">a;tax=k:Fungi,p:Ascomycota",
      "ACGT",
      ">b;tax=k:Straminipila,p:Oomycota",
      "ACGT",
      ">c;tax=k:Stramenopila,p:Oomycota",
      "ACGT"
    ),
    path
  )
  counts <- kingdom_counts(path, aliases = c(Straminipila = "Stramenopila"))
  expect_equal(counts[["Stramenopila"]], 2L)
  expect_equal(counts[["Fungi"]], 1L)
})

test_that("retained_kingdoms applies the threshold and the exclusions", {
  counts <- c(
    Fungi = 500L,
    Metazoa = 120L,
    Picozoa = 7L,
    cf.Fungi = 300L,
    Eukaryota.reg527 = 140L,
    `_mitochondrion` = 200L,
    Eukaryota_kgd_Incertae_sedis = 1348L,
    Unispike1 = 50L
  )
  expect_equal(retained_kingdoms(counts, min_records = 10L), "Metazoa")
})

test_that("water_fill spreads the total evenly and gives small kingdoms all their records", {
  available <- c(A = 100L, B = 50L, C = 3L, D = 40L)
  alloc <- water_fill(available, 40L)
  expect_equal(sum(alloc), 40L)
  expect_equal(alloc[["C"]], 3L)
  expect_equal(unname(alloc[c("A", "B", "D")]), c(13L, 12L, 12L))
  expect_equal(water_fill(available, 1000L), available)
})

test_that("derive_fungi_rep appends drawn non-fungal records in both formats", {
  skip_if_not(dbpq::is_vsearch_installed(), "vsearch not available")
  work <- withr::local_tempdir()
  old_dirs <- list(
    sintax_dir = sintax_dir,
    dada2_dir = dada2_dir,
    sources_dir = sources_dir
  )
  withr::defer(list2env(old_dirs, envir = globalenv()))
  list2env(
    list(
      sintax_dir = file.path(work, "sintax"),
      dada2_dir = file.path(work, "dada2"),
      sources_dir = file.path(work, "sources")
    ),
    envir = globalenv()
  )
  lapply(c(sintax_dir, dada2_dir, sources_dir), dir.create)

  set.seed(3)
  kingdoms <- c(
    rep("Fungi", 40),
    rep("Metazoa", 30),
    rep("Viridiplantae", 5),
    rep("_mitochondrion", 20),
    rep("Straminipila", 12)
  )
  seqs <- vapply(
    seq_along(kingdoms),
    \(i) {
      paste(
        sample(c("A", "C", "G", "T"), 20 + i, replace = TRUE),
        collapse = ""
      )
    },
    character(1)
  )
  # rec44 is rec43 with one substitution: a control at ~99 % identity, which
  # must not be drawn either (config.R::rep_max_identity_to_external).
  seqs[44] <- paste0("T", substring(seqs[43], 2))
  ids <- paste0("rec", seq_along(kingdoms))
  writeLines(
    c(rbind(
      paste0(">", ids, ";tax=k:", kingdoms, ",p:P,c:C,o:O,f:F,g:G,s:S"),
      seqs
    )),
    refseq_path("SRC", "sintax")
  )
  writeLines(
    c(rbind(paste0(">", kingdoms, ";P;C;O;F;G;S;"), seqs)),
    refseq_path("SRC", "dada2")
  )
  fungi <- kingdoms == "Fungi"
  writeLines(
    c(rbind(
      paste0(">", ids, ";tax=k:", kingdoms, ",p:P,c:C,o:O,f:F,g:G,s:S")[fungi],
      seqs[fungi]
    )),
    refseq_path("SRC_Fungi", "sintax")
  )
  writeLines(
    c(rbind(paste0(">", kingdoms, ";P;C;O;F;G;S;")[fungi], seqs[fungi])),
    refseq_path("SRC_Fungi", "dada2")
  )
  external <- file.path(work, "external.fasta")
  writeLines(c(">ext1", seqs[43]), external)

  derive_fungi_rep(
    "SRC",
    output_db = "SRC_Fungi_rep",
    share = 0.1,
    min_records = 10L,
    external_fasta = external,
    mini_n = 10L
  )

  out_sintax <- Biostrings::readDNAStringSet(refseq_path(
    "SRC_Fungi_rep",
    "sintax"
  ))
  out_dada2 <- Biostrings::readDNAStringSet(refseq_path(
    "SRC_Fungi_rep",
    "dada2"
  ))
  expect_length(out_sintax, 44)
  expect_equal(
    unname(as.character(out_sintax)),
    unname(as.character(out_dada2))
  )
  expect_equal(
    names(out_sintax)[1:40],
    paste0(ids[fungi], ";tax=k:Fungi,p:P,c:C,o:O,f:F,g:G,s:S")
  )
  rep_kingdoms <- sub("^.*tax=k:([^,]*),.*$", "\\1", names(out_sintax)[41:44])
  expect_equal(as.vector(table(rep_kingdoms)), c(2L, 2L))
  expect_setequal(unique(rep_kingdoms), c("Metazoa", "Straminipila"))
  expect_equal(sub(";.*$", "", names(out_dada2)[41:44]), rep_kingdoms)
  expect_false(seqs[43] %in% as.character(out_sintax))
  expect_false(seqs[44] %in% as.character(out_sintax))

  reps <- utils::read.csv(file.path(sources_dir, "SRC_Fungi_rep_reps.csv"))
  expect_equal(nrow(reps), 4)
  expect_setequal(reps$kingdom, c("Metazoa", "Stramenopila"))

  expect_length(
    Biostrings::readDNAStringSet(refseq_path(
      "SRC_Fungi_rep",
      "sintax",
      mini = TRUE
    )),
    11
  )
  expect_message(
    derive_fungi_rep(
      "SRC",
      output_db = "SRC_Fungi_rep",
      external_fasta = external
    ),
    "already exists"
  )

  # A failed rebuild keeps the previous files.
  rebuilt <- c(
    refseq_path("SRC_Fungi_rep", "sintax"),
    refseq_path("SRC_Fungi_rep", "dada2"),
    file.path(sources_dir, "SRC_Fungi_rep_reps.csv")
  )
  expect_error(derive_fungi_rep(
    "SRC",
    output_db = "SRC_Fungi_rep",
    share = 0.1,
    min_records = 10L,
    external_fasta = file.path(work, "absent.fasta"),
    mini_n = 10L,
    force = TRUE
  ))
  expect_true(all(file.exists(rebuilt)))

  # A `_Fungi` file newer than the output is reported, not rebuilt.
  Sys.setFileTime(refseq_path("SRC_Fungi", "sintax"), Sys.time() + 60)
  expect_warning(
    derive_fungi_rep(
      "SRC",
      output_db = "SRC_Fungi_rep",
      external_fasta = external
    ),
    "STALE"
  )
})
