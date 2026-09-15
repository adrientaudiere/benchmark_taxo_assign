# create_fake_pq_from_refseq() wraps a set of reference sequences in a
# degenerate phyloseq (one zero-count sample) so assign_*() can consume a CV
# fold. Run from the project root:
#   Rscript tests/test_create_fake_pq_from_refseq.R

library("testthat")
library("here")
library("phyloseq")

here::i_am("tests/test_create_fake_pq_from_refseq.R")
source(here("R/load_pqverse.R"))
load_pqverse("MiscMetabar")
source(here("R/create_fake_pq_from_refseq.R"))

seqs <- Biostrings::DNAStringSet(c(
  "ACGTACGTACGTACGTACGT",
  "TTGACCGGTTAACCGGTTAA",
  "GGGCCCAAATTTGGGCCCAA"
))
names(seqs) <- c(
  "ref1;tax=k:Fungi,p:Ascomycota,c:Sordariomycetes,o:Hypocreales,f:Nectriaceae,g:Fusarium,s:Fusarium_oxysporum",
  "ref2;tax=k:Fungi,p:Basidiomycota,c:Agaricomycetes,o:Agaricales,f:Agaricaceae,g:Agaricus,s:Agaricus_bisporus",
  "ref3;tax=k:Fungi,p:Ascomycota,c:Eurotiomycetes,o:Eurotiales,f:Aspergillaceae,g:Penicillium,s:Penicillium_sp"
)

test_that("taxonomy is parsed from sintax headers", {
  pq <- create_fake_pq_from_refseq(seqs)
  expect_s4_class(pq, "phyloseq")
  expect_equal(ntaxa(pq), 3)
  expect_equal(taxa_names(pq), names(seqs))
  expect_equal(rank_names(pq),
               c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
  expect_true(grepl("Fusarium", tax_table(pq)[1, "Genus"]))
  expect_true(grepl("Agaricales", tax_table(pq)[2, "Order"]))
  expect_equal(nsamples(pq), 1)
  expect_equal(sum(otu_table(pq)), 0)
  expect_equal(as.character(refseq(pq)[[2]]), as.character(seqs[[2]]))
})

test_that("dada2-format headers (UNITE and EUKARYOME styles) are parsed too", {
  dada2_seqs <- seqs
  names(dada2_seqs) <- c(
    "Fusarium_oxysporum|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum",
    "Fungi;Basidiomycota;Agaricomycetes;Agaricales;Agaricaceae;Agaricus;Agaricus_bisporus;",
    "Fungi;Ascomycota;Eurotiomycetes;Eurotiales;Aspergillaceae;Penicillium;"
  )
  pq <- create_fake_pq_from_refseq(dada2_seqs)
  tt <- as.matrix(unclass(tax_table(pq)))
  expect_equal(unname(tt[1, ]), c("Fungi", "Ascomycota", "Sordariomycetes", "Hypocreales",
                                  "Nectriaceae", "Fusarium", "Fusarium_oxysporum"))
  expect_equal(unname(tt[2, "Kingdom"]), "Fungi")
  expect_equal(unname(tt[2, "Species"]), "Agaricus_bisporus")
  expect_equal(unname(tt[3, "Genus"]), "Penicillium")
  expect_true(tt[3, "Species"] %in% c("", NA))
})

test_that("headers carry the taxonomy when the names are identifiers", {
  ided <- seqs
  names(ided) <- c("id1", "id2", "id3")
  pq <- create_fake_pq_from_refseq(ided, headers = c(
    "Fungi;Ascomycota;Sordariomycetes;Hypocreales;Nectriaceae;Fusarium;Fusarium_oxysporum;",
    "Fungi;Basidiomycota;Agaricomycetes;Agaricales;Agaricaceae;Agaricus;Agaricus_bisporus;",
    "Fungi;Basidiomycota;Agaricomycetes;Agaricales;Agaricaceae;Agaricus;Agaricus_bisporus;"
  ))
  tt <- as.matrix(unclass(tax_table(pq)))
  expect_equal(taxa_names(pq), c("id1", "id2", "id3"))
  expect_equal(unname(tt[1, "Genus"]), "Fusarium")
  expect_equal(unname(tt[3, "Species"]), "Agaricus_bisporus")
})

test_that("parse_header_taxonomy keeps sintax behaviour and pads missing ranks", {
  m <- parse_header_taxonomy(c("id1;tax=k:Fungi,p:Ascomycota", "A;B;C;"), 4)
  expect_equal(m[1, ], c("k:Fungi", "p:Ascomycota", "", ""))
  expect_equal(m[2, ], c("A", "B", "C", ""))
})

test_that("a fasta path and taxonomy_in_names = FALSE are accepted", {
  path <- tempfile(fileext = ".fasta")
  Biostrings::writeXStringSet(seqs, path)
  pq <- create_fake_pq_from_refseq(path, taxonomy_in_names = FALSE)
  expect_equal(ntaxa(pq), 3)
  expect_true(all(as.matrix(unclass(tax_table(pq))) == "FAKE"))
})
