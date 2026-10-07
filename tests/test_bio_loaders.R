# Loaders of the published D1c datasets (R/bio_loaders.R, OD 30). Run from the
# project root:
#   Rscript tests/test_bio_loaders.R

library("testthat")
library("here")

here::i_am("tests/test_bio_loaders.R")
source(here("R/bio_loaders.R"))

test_that("common_units() keeps the units occurring in min_occur samples", {
  otu <- matrix(
    c(5, 1, 0, 0, 3, 0, 0, 0, 2),
    nrow = 3,
    dimnames = list(c("S1", "S2", "S3"), c("OTU_b", "OTU_c", "OTU_a"))
  )
  expect_equal(common_units(otu), c("OTU_a", "OTU_b", "OTU_c"))
  expect_equal(common_units(otu, min_occur = 2), "OTU_b")
  expect_equal(common_units(otu, min_occur = 4), character())
})

test_that("load_gssp_air() keeps the common OTUs and drops empty samples", {
  dir <- tempfile()
  dir.create(dir)
  writeLines(
    c(
      "OTU;nsample;nread;kingdom;phylum;class;order;family;genus;species;sequence",
      "OTU1;1;5;Fungi;Asco;C;O;F;Gen_1;pseudospecies_1;ACGTACGT",
      "OTU2;1;3;Fungi;Basidio;C;O;F;Gen_2;Gen_two_9;GGGTTTAA",
      "OTU3;1;2;Fungi;Asco;C;O;F;Gen_3;pseudospecies_3;TTTTCCCC"
    ),
    file.path(dir, "taxonomy.csv")
  )
  writeLines(
    c(
      '"";"OTU1";"OTU2";"OTU3"',
      '"S1";5;0;0',
      '"S2";1;3;0',
      '"S3";0;0;2'
    ),
    file.path(dir, "otu.table.csv")
  )
  writeLines(
    c(
      '"sample.id";"site";"lat"',
      '"S1";"A";52,5',
      '"S2";"B";10,25',
      '"S3";"C";-3,5'
    ),
    file.path(dir, "metadata.csv")
  )
  all_otu <- load_gssp_air(dir)
  expect_equal(phyloseq::taxa_names(all_otu), c("OTU1", "OTU2", "OTU3"))
  expect_equal(
    as.character(all_otu@refseq),
    c(OTU1 = "ACGTACGT", OTU2 = "GGGTTTAA", OTU3 = "TTTTCCCC")
  )
  pq <- load_gssp_air(dir, min_occur = 2)
  expect_equal(phyloseq::taxa_names(pq), "OTU1")
  expect_equal(phyloseq::sample_names(pq), c("S1", "S2"))
  expect_equal(colnames(pq@tax_table)[c(1, 7)], c("Kingdom", "Species"))
  expect_type(pq@sam_data$lat, "double")
})

test_that("load_korhonen_wood() pools both deposits on the hash and drops absent OTUs", {
  dir <- tempfile()
  dec <- file.path(dir, "doi_10_5061_dryad_3j9kd51s7__v1")
  spr <- file.path(dir, "doi_10_5061_dryad_44j0zpcng__v1")
  dir.create(dec, recursive = TRUE)
  dir.create(spr, recursive = TRUE)
  put <- function(folder, file, lines) {
    writeLines(
      lines,
      file.path(folder, paste0("Korhonen_et_al_2024_", file, ".csv"))
    )
  }
  put(
    dec,
    "Fungi_ITS2_OTU_taxonomy",
    c(
      "OTU_name;OTU_ID;RepSeq;k;k%;p;p%;c;c%;o;o%;f;f%;g;g%;s;s%",
      "OTU_1;>h1;AAAACCCC;k__Fungi;100;p__Asco;90;;;;;;;;;;",
      "OTU_2;>h2;CCCCGGGG;k__Fungi;100;p__Basidio;90;;;;;;;;;;",
      "OTU_3;>h3;GGGGTTTT;k__Fungi;100;;;;;;;;;;;;"
    )
  )
  put(
    dec,
    "Fungi_presabs_DFWD_replicates",
    c("SampleNo;OTU_1;OTU_2;OTU_3", "1;1;0;0", "2;0;1;0")
  )
  put(
    dec,
    "Fungi_DFWD_SampleNo_SiteNo_table",
    c(
      "SampleNo;SiteNo;SampleCode;Valid/Outlier",
      "1;7;7-A1;valid",
      "2;7;7-A2;outlier"
    )
  )
  put(dec, "Fungi_presabs_BirchLog_samples", c("BirchLogNo;OTU_2", "1;1"))
  put(
    dec,
    "BirchLog_sample_metadata",
    c("BirchLogNo;SampleCode;SiteNo", "1;02-KO;16")
  )
  put(
    spr,
    "Fungal_OTU_taxonomy_and_metadata",
    c(
      "OTU_ID;RepSeq;k;k%;p;p%;c;c%;o;o%;f;f%;g;g%;s;s%;Comment",
      "h2;CCCCGGGG;Fungi;100;Basidio;90;;;;;;;;;;;",
      "h4;TTTTAAAA;Fungi;100;Mucoro;90;;;;;;;;;;;"
    )
  )
  put(
    spr,
    "Fungal_OTU_presabs",
    c("SampleID;h2;h4", "Fungi_Trunk_03;1;1", "Fungi_Top_03;0;0")
  )

  pq <- load_korhonen_wood(dir)
  expect_equal(sort(phyloseq::taxa_names(pq)), c("h1", "h2", "h4"))
  expect_equal(
    sort(phyloseq::sample_names(pq)),
    c("BirchLog_01", "DFWD_001", "DFWD_002", "Fungi_Trunk_03")
  )
  expect_equal(as.character(pq@tax_table["h1", "Phylum"]), "Asco")
  expect_equal(as.character(pq@tax_table["h4", "Phylum"]), "Mucoro")
  expect_equal(as.character(pq@refseq["h2"]), c(h2 = "CCCCGGGG"))
  sam <- as(phyloseq::sample_data(pq), "data.frame")
  expect_equal(sam["DFWD_002", "outlier"], TRUE)
  expect_equal(sam["DFWD_001", "site"], "site_7")
  expect_equal(sam["Fungi_Trunk_03", "substrate"], "spruce_trunk")
  expect_equal(sam["Fungi_Trunk_03", "site"], "tree_3")
})
