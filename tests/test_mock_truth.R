# Per-unit truth of a mock community (R/mock_truth.R, ROADMAP 0.6). The name
# helpers and the lineage reader are tested without vsearch; the matching rule
# needs vsearch and is skipped without it. Run from the project root:
#   Rscript tests/test_mock_truth.R

library("testthat")
library("here")

here::i_am("tests/test_mock_truth.R")
source(here("R/load_pqverse.R"))
load_pqverse("dbpq")
source(here("R/mock_truth.R"))

set.seed(7)
rand_seq <- function(n) {
  paste(sample(c("A", "C", "G", "T"), n, replace = TRUE), collapse = "")
}
mutate_one <- function(seq, pos = 150) {
  substr(seq, pos, pos) <- c(A = "C", C = "G", G = "T", T = "A")[
    substr(seq, pos, pos)
  ]
  seq
}

strain_seqs <- c(
  s_alpha = rand_seq(400),
  s_beta = rand_seq(400),
  s_gamma = rand_seq(400)
)
strain_seqs[["s_delta"]] <- strain_seqs[["s_gamma"]]

sanger <- tibble::tibble(
  strain = names(strain_seqs),
  name = c(
    "Alpha_unus",
    "Beta_duo",
    "Gamma_tres",
    "Gamma_quattuor"
  ),
  sequence = unname(strain_seqs)
)

lineage <- tibble::tibble(
  strain = sanger$strain,
  name = sanger$name,
  Kingdom = "Fungi",
  Phylum = c("Ascomycota", "Ascomycota", "Basidiomycota", "Basidiomycota"),
  Class = c(
    "Sordariomycetes",
    "Sordariomycetes",
    "Agaricomycetes",
    "Agaricomycetes"
  ),
  Order = c("Hypocreales", "Hypocreales", "Agaricales", "Agaricales"),
  Family = c("Nectriaceae", "Nectriaceae", "Agaricaceae", "Agaricaceae"),
  Genus = genus_name(sanger$name),
  Species = binomial_name(sanger$name),
  Genus_accepted = c("Alpha", "Betula", "Gamma", "Gamma"),
  Species_accepted = c(
    "Alpha_unus",
    "Betula_duo",
    "Gamma_tres",
    "Gamma_quattuor"
  )
)

test_that("binomial_name and genus_name read both spellings", {
  expect_equal(
    binomial_name(c("Sidera americana", "Agaricus_essettei")),
    c("Sidera_americana", "Agaricus_essettei")
  )
  expect_equal(
    binomial_name("Fusarium fujikuroi species complex"),
    "Fusarium_fujikuroi"
  )
  expect_equal(
    genus_name(c("Sidera americana", "Agaricus_essettei")),
    c("Sidera", "Agaricus")
  )
  expect_true(is.na(binomial_name(NA_character_)))
})

test_that("species_name drops the names that identify no species", {
  # "Phoma sp" is an identification to the genus: the unit must stop at Genus,
  # not carry a species truth no assignment can match (chapter 01 turns every
  # "<Genus>_sp" assignment into NA too).
  expect_equal(
    species_name(c("Phoma sp", "Alternaria sp.", "Fusarium spp", "Phoma")),
    rep(NA_character_, 4)
  )
  expect_equal(
    species_name(c("Sidera americana", "Agaricus_essettei")),
    c("Sidera_americana", "Agaricus_essettei")
  )
  expect_true(is.na(species_name(NA_character_)))
})

test_that("mock_lineage stops on a name missing from the cache", {
  cache <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      submitted_name = "Alpha unus",
      accepted_name = "Alpha unus",
      Kingdom = "Fungi",
      Phylum = "Ascomycota",
      Class = "Sordariomycetes",
      Order = "Hypocreales",
      Family = "Nectriaceae",
      Genus = "Alpha"
    ),
    cache,
    row.names = FALSE
  )
  expect_error(mock_lineage(sanger, cache), "missing from")
  one <- mock_lineage(sanger[1, ], cache)
  expect_equal(one$Species, "Alpha_unus")
  expect_equal(one$Family, "Nectriaceae")
  expect_error(
    mock_lineage(sanger, tempfile()),
    "run refresh_gna_lineage_cache"
  )
})

test_that("lineage_disagreements lists the ranks where the authors' table differs", {
  published <- data.frame(
    Class = c(
      "Sordariomycetes",
      "Dothideomycetes",
      "Agaricomycetes",
      "Agaricomycetes"
    ),
    Family = c("Nectriaceae", "Nectriaceae", "Agaricaceae", "")
  )
  dis <- lineage_disagreements(lineage, published)
  expect_equal(nrow(dis), 1)
  expect_equal(dis$rank, "Class")
  expect_equal(dis$strain, "s_beta")
  expect_equal(dis$published, "Dothideomycetes")
  expect_equal(dis$gbif, "Sordariomycetes")
})

test_that("a unit takes the truth of the strain it matches, in both inclusion directions", {
  skip_if_not(dbpq::is_vsearch_installed(), "vsearch not available")
  units <- Biostrings::DNAStringSet(c(
    exact = strain_seqs[["s_alpha"]],
    one_error = mutate_one(strain_seqs[["s_beta"]]),
    inside = substr(strain_seqs[["s_alpha"]], 30, 380),
    tie = strain_seqs[["s_gamma"]],
    unrelated = rand_seq(400)
  ))
  truth <- mock_truth_table(units, sanger, lineage)

  expect_equal(truth$taxon, names(units))
  expect_equal(truth$truth_depth[truth$taxon == "exact"], "Species")
  expect_equal(truth$Species[truth$taxon == "exact"], "Alpha_unus")
  expect_equal(truth$Species_accepted[truth$taxon == "exact"], "Alpha_unus")

  # 399 of 400 bases identical: above the 99.5 % identity of the rule.
  expect_equal(truth$truth_depth[truth$taxon == "one_error"], "Species")
  expect_equal(truth$Species_accepted[truth$taxon == "one_error"], "Betula_duo")

  # The unit sits inside the Sanger read: covered at 100 %, the read is not.
  expect_equal(truth$truth_depth[truth$taxon == "inside"], "Species")
  expect_equal(truth$Species[truth$taxon == "inside"], "Alpha_unus")

  # Two strains share this sequence and differ at Species: truth stops at Genus
  # and the unit leaves the matrix below it.
  tie <- truth[truth$taxon == "tie", ]
  expect_equal(tie$n_strains, 2L)
  expect_equal(tie$truth_depth, "Genus")
  expect_equal(tie$Genus, "Gamma")
  expect_true(is.na(tie$Species))
  expect_equal(tie$Genus_accepted, "Gamma")
  expect_true(is.na(tie$Species_accepted))

  # No match: the unit stays in the table with no truth (Q2b).
  none <- truth[truth$taxon == "unrelated", ]
  expect_equal(none$n_strains, 0L)
  expect_true(is.na(none$truth_depth))
  expect_true(all(is.na(none[truth_ranks])))
})

test_that("a unit matching two strains of the same species is not a tie", {
  skip_if_not(dbpq::is_vsearch_installed(), "vsearch not available")
  same_species <- sanger
  same_species$name[4] <- "Gamma_tres"
  same_lineage <- lineage
  same_lineage$name[4] <- "Gamma_tres"
  same_lineage$Species[4] <- "Gamma_tres"
  same_lineage$Species_accepted[4] <- "Gamma_tres"

  truth <- mock_truth_table(
    Biostrings::DNAStringSet(c(both = strain_seqs[["s_gamma"]])),
    same_species,
    same_lineage
  )
  expect_equal(truth$n_strains, 2L)
  expect_equal(truth$truth_depth, "Species")
  expect_equal(truth$Species, "Gamma_tres")
})

test_that("a stricter cover rule drops the unit included in a longer Sanger read", {
  skip_if_not(dbpq::is_vsearch_installed(), "vsearch not available")
  units <- Biostrings::DNAStringSet(c(
    inside = substr(strain_seqs[["s_alpha"]], 30, 380)
  ))
  strict <- mock_truth_table(units, sanger, lineage, min_cover = 1.001)
  expect_true(is.na(strict$truth_depth))
})

test_that("external_control_truth reads the lineage from the control names", {
  taxa <- c(
    "Taxa_1",
    "fake_3",
    paste0(
      "external_Oxytrichidae_sp|UDB0661760|SH1048417.10FU|reps|",
      "k__Alveolata;p__Ciliophora;c__Spirotrichea;o__Sporadotrichida;",
      "f__Oxytrichidae;g__Oxytrichidae_gen_Incertae_sedis;s__Oxytrichidae_sp"
    ),
    "external_no_lineage"
  )
  ext <- external_control_truth(taxa)
  expect_equal(nrow(ext), 2)
  expect_equal(ext$Kingdom, c("Alveolata", NA))
  expect_equal(ext$Order[1], "Sporadotrichida")
  # "s__Oxytrichidae_sp" names no species: the control has no species truth.
  expect_true(is.na(ext$Species[1]))
  expect_equal(ext$Genus[1], "Oxytrichidae_gen_Incertae_sedis")
  expect_true(all(is.na(ext[2, truth_ranks])))
})

test_that("a strain named to the genus only stops the truth at Genus", {
  skip_if_not(dbpq::is_vsearch_installed(), "vsearch not available")
  genus_only <- sanger[1, ]
  genus_only$name <- "Alpha sp"
  genus_lineage <- lineage[1, ]
  genus_lineage$name <- "Alpha sp"
  genus_lineage$Species <- species_name("Alpha sp")
  genus_lineage$Species_accepted <- species_name("Alpha sp")
  truth <- mock_truth_table(
    Biostrings::DNAStringSet(c(exact = strain_seqs[["s_alpha"]])),
    genus_only,
    genus_lineage
  )
  expect_equal(truth$truth_depth, "Genus")
  expect_equal(truth$Genus, "Alpha")
  expect_true(is.na(truth$Species))
})
