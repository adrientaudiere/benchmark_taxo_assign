# 07 — Validity of the stored results

Found on 2026-09-10 while validating the restructured pipelines on a
`assign_taxo_mini` run and rendering the analysis chapters on it. Unlike
critiques 01–06, these defects change **scientific results**, not only
maintainability. All of them predate the 2026-09-10 restructuring; the
restructuring copied B12 faithfully and made the others visible.

ROADMAP ids in brackets.

## Summary

| Id | Defect | Since | Affects | Code fixed | Data / rerun needed |
|---|---|---|---|---|---|
| B12 | `values_map`: `ifelse(mini_db, ...)` recycles its first value, so **every target used the first database** | commit `17ecc9e`, 2026-05-26 | whole DB axis of `store_assign_taxo` and `store_cross_val` | yes | rerun both pipelines |
| B13 | UNITE sintax files lose the kingdom (`…\|k__X;tax=p:…`), so sintax / lca / blastn ranks are **shifted by one** | DB derivation of 2026-05-26 | every sintax / lca / blastn result on `Unite`, `Unite_Fungi` and their `mini_*` | yes (`derive_sintax()`) | regenerate 4 fasta files, then rerun |
| B14 | cross-validation truth table empty for dada2-format references, so **dada2 scores 0** | since `cross_val()` exists | dada2 rows of `store_cross_val` | yes | rerun CV |
| B15 | cross-validation lca ranks named `Class_sintax` etc. | since `cross_val()` exists | lca rows of `store_cross_val` (figures facet them apart) | yes, plus normalisation in chapter 06 | rerun CV (or keep the normalisation) |

Consequence: **no figure or table computed from `store_assign_taxo` or
`store_cross_val` since 2026-05-26 can be used**, including the Q1, Q2, Q3,
M2 and D1a items ticked in the ROADMAP done log. The 2025-02 sequential run
predates B12; whether its UNITE sintax file already had the B13 header
problem is unknown (the reference files were regenerated on 2026-05-26).

## B12 — every target used the first database [S1.5]

`db_name = ifelse(mini_db, paste0("mini_", db), db)` inside `mutate()`:
`mini_db` is a scalar, so `ifelse()` returns a single value, recycled to all
rows.

```
# A tibble: 3 × 3
  db                   mini_TRUE  mini_FALSE
1 Unite                mini_Unite Unite
2 EUK_ITS_v2           mini_Unite Unite
3 EUK_SSU_v2_Fungi_cut mini_Unite Unite
```

Evidence in the stores:

- lca and blastn `Genus_*` columns are **identical** across the six
  databases, in both `store_assign_taxo` (2026-05) and
  `store_assign_taxo_mini` (2026-09-10); sintax and dada2 columns differ only
  by their stochastic part (`--sintax_random`, bootstraps).
- Re-running sintax on five mock ASVs against `mini_Unite` reproduces the
  stored `sintax__EUK_ITS_v2___0.5` column exactly; against `mini_EUK_ITS_v2`
  it does not.
- dada2 "EUK" columns carry `k__`-prefixed values although the EUKARYOME
  dada2 files have no prefixes.

Fix: `if (mini_db) ... else ...` in both builders of `R/values_map.R`;
`tests/test_values_map.R` now requires one distinct reference file per
database in both modes (the previous test only checked a single database).

## B13 — UNITE sintax conversion drops the kingdom [S1.6]

`dbpq::format2sintax()` expects `>ID;k__Kingdom;p__...`. UNITE
general-release headers are `>Name|Acc|SH|type|k__Kingdom;p__...`, so the
parser takes `Name|Acc|SH|type|k__Kingdom` as the identifier:

```
input : Fusarium_sp|KF1|SH1.10FU|refs|k__Fungi;p__Ascomycota;c__Sordariomycetes
output: Fusarium_sp|KF1|SH1.10FU|refs|k__Fungi;tax=p:Ascomycota,c:Sordariomycetes
```

All 266 589 records of `sintax_format/Unite.fasta` and 168 049 of
`Unite_Fungi.fasta` have this form; the EUKARYOME sintax files are correct
(`;tax=k:` in every record). MiscMetabar's `assign_sintax()`,
`assign_vsearch_lca()` and `assign_blastn()` split the prediction string
**by position**, so a missing first rank shifts every column: `Kingdom_*`
holds phyla, `Genus_*` holds species, `Species_*` is empty. On the mini run
the mock F1 of sintax, lca and blastn is 0 at Kingdom, Phylum and Genus,
while dada2 (which reads the dada2-format file) scores 1, 1 and 0.55.

Fix applied: `derive_sintax()` in `make_databases.R` rewrites `|k__` into
`;k__` before calling `format2sintax()` (tested on a UNITE-style fasta in
`tests/test_make_databases_helpers.R`).

**Update 2026-09-11.** `dbpq::format2dada2()` has the same defect on UNITE
general headers (its dada2 output starts at the phylum), so `derive_dada2()`
applies the same rewrite (`general_headers_fixed()`, tested). The files are
no longer patched in place: `make_databases.R` now downloads the general
FASTA releases and derives every file under new versioned names
(`Unite_s_all_20250219`, `EUK_ITS_v2.1`, …; see
`docs/reference_databases.md` and ROADMAP S7.x). The old unversioned files
are left on disk as legacy and nothing reads them.

Upstream, two reports are worth filing:

- **dbpq** — `format2sintax()` should accept UNITE general-release headers
  (identifier fields separated by `|`).
- **MiscMetabar** — the three `assign_*()` functions should parse the
  prediction by rank letter (`k:`, `p:`, …) rather than by position, so a
  missing rank yields NA instead of silently shifting every column.

## B14 — dada2 truth table empty in cross-validation [S1.4]

`create_fake_pq_from_refseq()` only parsed `;tax=` headers. The dada2
branch of `cross_val()` reads `dada2_format/` files, whose headers are
`…|k__Fungi;p__…` (UNITE) or `Fungi;Ascomycota;…;` (EUKARYOME): the truth
table was entirely empty and every comparison failed. In `store_cross_val`,
dada2 `good_classifications` is 0.000 for every database and rank.

Fix: `parse_header_taxonomy()` in `R/create_fake_pq_from_refseq.R` reads
both styles; the dada2 branch of `cross_val()` also strips `k__` prefixes
from `assignTaxonomy()` output. Tested with both header styles.

## B15 — lca rank names in cross-validation [S1.4]

`assign_vsearch_lca()` names its columns `<rank>_sintax` (default
`suffix = "_sintax"`), so `cv_results` holds `Class_sintax` … for lca. In
chapter 06, those names became NA factor levels and produced duplicated keys
(the render error "argument non numérique pour un opérateur binaire").

Fix: `cross_val()` renames the lca columns; chapter 06 normalises the
suffix for stores built before the fix and asserts that no rank is NA.

## B16 — "Fungi" filter lets non-Fungi records in [S7.2]

`derive_pattern_only(pattern = "Fungi")` kept every record whose header
contains `Fungi` anywhere. Records whose kingdom is not Fungi, in the files
built on 2026-05-26:

| File | Records | Kingdom not Fungi | Examples |
|---|---|---|---|
| `sintax_format/EUK_ITS_v2_Fungi.fasta` | 1 045 718 | 965 | `k:Metazoa,…,f:Fungiidae` (corals) |
| `sintax_format/EUK_SSU_v2_Fungi.fasta` | 97 790 | 22 | `k:cf.Fungi`, `k:_pseudogene.Fungi` |
| `sintax_format/Unite_Fungi.fasta` | 168 049 | 19 | `k__Metazoa;…;f__Fungiidae`, `Baikalospongia_fungiformis` |

Small in proportion, but they are wrong answers the classifiers can return
with high confidence. Fix: `derive_kingdom_only()` anchors the pattern on the
kingdom field of each header format (`tax=k:Fungi,`, `^>Fungi;`,
`k__Fungi;`), tested with Fungiidae and `cf.Fungi` records.

## B17 — EUKARYOME name qualifiers [S7.8]

Found on the 2026-09-11 `assign_taxo_mini` run: the sintax targets on
`EUK_ITS_v2.1_Fungi_cut` (0.4, 0.5, 0.6) and `EUK_SSU_v2.1_Fungi_cut` (0.6)
warned "NAs introduits lors de la conversion automatique" in
`assign_sintax()`. EUKARYOME qualifies homonyms (`g__Lactarius(Fungi)`,
`f__Gonostomatidae(Sporadotrichida)`), brackets some genera (`g__(Candida)`,
`g__(Candida]`) and adds `.s.str` or `.nom.prov` suffixes. v2.1 has 11 413
ITS and 5 586 SSU headers with parentheses; v2.0 had 9 641 in ITS, so the
2026-05 stores were affected too. Only the `Fungi_cut` mini files warned
because their first 10 000 records happen to include such genera; every full
EUKARYOME database does.

- **sintax**: vsearch prints `g:Lactarius(Fungi)(0.97)`. `assign_sintax()`
  splits on the first `(` (`too_many = "drop"`): the genus is `Lactarius` but
  its bootstrap is NA, so `min_bootstrap` never removes it (optimistic bias
  at genus rank). `g:(Candida)` yields an empty genus.
- **dada2 and lca** return `Lactarius(Fungi)` (8 genus columns of the mini
  `d_all_taxo`), which never matches the mock.
- **Mock genera concerned**: `Lactarius` and `Cryptococcus` (`(Fungi)`),
  `Mortierella` and `Peziza` (`.s.str`): 4 of 98.

Fix: `general_headers_fixed()` in `make_databases.R` keeps the name only
before the dbpq conversion (tested on real header forms). The EUKARYOME dada2
and sintax files must be rebuilt with `force = TRUE`, and the `EUK_` targets
invalidated, since database paths were not file targets yet (they are since S7.9). Upstream: a
dbpq report (proposals_for_dbpq.md §0).

## B18 — blastn query cover and the 18S flank of the ASVs [S1.3]

`assign_blastn()` keeps a hit only if its alignments cover `min_cover = 95` %
of the query (BLAST `qcovs`): the share of the ASV included in the reference,
whatever the reference length. Neither `pipelines/assign_taxo.R` nor
`cross_val()` changes this default. The mock ASVs (ITS1F–ITS2, primers
removed) start with 30–45 bp of the 18S end (`AGTCGTAACAAGG…GGATCATTA`, motif
at median position 37 in 196 of 200 ASVs), while ITS references start at
ITS1, so that flank never aligns.

Best hit of each ASV against the full `Unite_s_all_20250219_Fungi` (sintax
format, 2026-09-11): median identity 98.8 %, median query cover 88 %; 156 of
200 below 95 %, all with the unaligned part at the ASV start (median 32 bp),
147 of them aligned from reference position 1.

| Full database | ASVs with a hit at ≥ 95 % id, bit ≥ 50, e ≤ 1e-30 | … and cover ≥ 95 % | Median best cover |
|---|---|---|---|
| `Unite_s_all_20250219_Fungi` | 200 | 44 | 88 % |
| `EUK_ITS_v2.1_Fungi` | 200 | 88 | 93 % |
| `EUK_ITS_v2.1_Fungi_cut` | 123 | 84 | 100 % |

On the `mini_*` databases the same filters leave 1 ASV (UNITE), which is why
the smoke tests show almost no blastn assignment. With `min_cover = 95`, the
blastn results on uncut databases depend on how much 18S flank their records
keep (EUKARYOME more than UNITE) as much as on the classifier.

Fix (decision 18, revised 2026-09-11): `min_cover = 80`
(`config.R::blastn_min_cover`) for the blastn assignment and CV targets; 194
of 200 ASVs keep a hit at that cover on the full UNITE Fungi database. ITS
extraction of the ASVs with ITSx is proposed as a new benchmark question (Q4).

## B19 — sintax and lca never queried the negative controls [S6.7]

Found on 2026-09-11 while splitting the assignments into computations (S8.2).
`assign_sintax()` and `assign_vsearch_lca()` write their query fasta with
`write_temp_fasta(clean_pq = TRUE)` by default, and `clean_pq()` removes the
taxa without reads. The negative controls (`fake_*` shuffled sequences,
`external_*` non-Fungi sequences) have no reads by construction, so these two
methods never saw them:

- `assign_sintax(d_asv_for_assignation, behavior = "return_matrix")` on the
  mini store queried 200 of the 400 taxa (0 controls); with
  `clean_pq = FALSE`, 400 (200 controls).
- Kingdom assigned to the controls in the mini `d_all_taxo`: dada2 94 % of
  `external_*` and 79 % of `fake_*`; blastn 3.5 % and 0 %; sintax and lca 0 %
  and 0 %.

The controls are the TN denominator of `tc_metrics_mock()`: for sintax and
lca every control counted as a correct rejection, so their TN, TNR,
specificity, MCC and ACC were overstated in every store so far, and the
comparison with dada2 and blastn was biased in their favour.

Fix: `compute_assignment()` (`R/assign_compute.R`) calls both methods with
`clean_pq = FALSE`; sintax and lca now query 400 sequences instead of 200.
Upstream: a MiscMetabar note (S6.6), since `clean_pq = TRUE` silently changes
the query set of `assign_*()` when a phyloseq object holds empty taxa.

## B20 — MCC lost to integer overflow [S6.8]

Found on 2026-09-14 by rendering the chapters on the mini stores with R
warnings shown (`_quarto.yml` hides them): 14 × "NA produit par débordement
d'entier par le haut" in `(TP + FP) * (TP + FN) * (FP + TN) * (TN + FN)`.
`comparpq::tc_metrics_mock_vec()` counts with `sum()`, which returns integers,
and the product of the four margins exceeds `.Machine$integer.max` as soon as
each margin is above about 216 taxa, which the 200 ASVs + 199 controls reach.

- Mini `res_comp_tax`: 223 of 1344 MCC values are `NA`; 14 are overflows
  (finite in double precision, from −0.40 to 0.45, e.g.
  `rel_majority_consensus` at `Gen_sp`, `itsx_dada2__Unite_all_20250219_Fungi___0.4`
  at `Genus`), 209 have a zero margin (MCC undefined, a legitimate `NA`).
- The overflow hits the rows with the most balanced confusion matrices, so the
  missing MCC values are not missing at random: Q1–Q3 MCC summaries were biased.

Fix: comparpq `tc_metrics_mock_vec()` computes the MCC in double precision
(dev checkout: code, regression test, NEWS). The 209 zero-margin cells (120
ranks where the method assigns nothing, TP + FP = 0; 80 Kingdom cells where
every control gets a kingdom, TN + FN = 0; 9 Phylum cells, FP + TN = 0) now
get MCC = 0 (Chicco & Jurman 2020, developer decision of 2026-09-14): dropping
them would have removed the failures of a method from its mean.

## B21 — The `_Fungi_cut` databases are not amplicon regions [S6.9]

Found on 2026-09-15 after the developer noticed that cutadapt drops far more
records than expected. `dbpq::cutadapt_rm_primers_db()` runs
`cutadapt -e 0.1 -a '<primer_fw>...<primer_rev>' --discard-untrimmed` with
the reverse primer as written. On a forward-oriented reference the ITS2 site
is its reverse complement (`GCATCGATGAAGAACGCAGC`), so the 3′ adapter never
matches. Synthetic records with cutadapt 2.6: a record is kept as soon as one
adapter matches, and only a record carrying the reverse primer *as written*
is trimmed at both ends.

| | `EUK_ITS_v2.1_Fungi` → `_cut` | `EUK_SSU_v2.1_Fungi` → `_cut` |
|---|---|---|
| records | 1 076 870 → 51 602 (4.8 %) | 141 884 → 116 746 (82.3 %) |
| median length | 523 → 558 bp | 1218 → 46 bp |
| `_cut` records still holding the ITS2 site | 86.9 % | 0.2 % |

Primer sites on 3000 random records (≤ 2 mismatches): ITS1F in 0.1 %
(`Unite_s_all_20250219_Fungi`), 2.0 % (`EUK_ITS_v2.1_Fungi`), 2.7 %
(`EUK_ITS_v2.1`); reverse complement of ITS2 in 94.5 %, 93.7 % and 78.4 %;
reverse primer as written in 1 of 3000. The reference ITS records start inside
the ITS region, after ITS1F, so `--discard-untrimmed` keeps only the few that
start in the 18S; the SSU records end at the 18S, so trimming at ITS1F leaves
the 45–47 bp of 18S downstream of the primer.

Consequences: every Q2 result on `EUK_ITS_v2.1_Fungi_cut` compares a
primer-biased 5 % subset, not a trimmed database; `EUK_SSU_v2.1_Fungi_cut` is
a database of 46 bp 18S fragments, which can only match the 18S start of the
ASVs (B18); the 100 % blastn cover on `EUK_ITS_v2.1_Fungi_cut` (B18) and the
`_cut` genus-retention figures of S8.4 come from these files. Decision 7
assumed the ITS region lies inside the SSU gene; it does not.

Fix (2026-09-15, developer decisions): dbpq `cutadapt_rm_primers_db()`
reverse-complements `primer_rev` by default and marks each primer required or
optional explicitly; this project trims at both sites when found and keeps the
records without them (`config.R::cut_discard_untrimmed = FALSE`,
`primer_min_overlap = 20`). On 3000 random records every record is kept and
94.1 % (`EUK_ITS_v2.1_Fungi`) or 95.4 % (`Unite_s_all_20250219_Fungi`) are
trimmed, to a median of 211 or 216 bp. The SSU database is removed (decision 7
revised). `EUK_ITS_v2.1_Fungi_cut` must be rebuilt with `force = TRUE`.
Scripts: `check_cut_loss2.R`, `check_primer_sites.R`, `verify_cut_fix.R`
(session scratchpad).

## B22 — dada2 cross-validation draws one record per taxonomy string [D1a]

Found on 2026-09-15 while smoke-testing the trimmed CV queries. `cross_val()`
runs `dna[!duplicated(names(dna))]` before drawing the subsample. dada2 reads
`dada2_format/`, whose headers are the taxonomy only, so the deduplication
keeps the first record of each taxonomy string; sintax, lca and blastn read
`sintax_format/`, whose headers start with a record identifier.

| database | records | dada2 records left |
|---|---|---|
| `Unite_s_all_20250219_Fungi` | 168 030 | 38 310 (22.8 %) |
| `EUK_ITS_v2.1_Fungi` | 1 076 870 | 61 710 (5.7 %) |
| `EUK_ITS_v2.1` | 1 877 003 | 206 399 (11 %) |
| `mini_EUK_ITS_v2.1_Fungi_cut` ∩ `_Fungi` | 10 000 | 44 |

Consequences: dada2 is cross-validated on another population of records than
the three other methods; in the standard variant every training record has a
distinct taxonomy string, so the exact taxonomy of a query is never in its
training part. Every dada2 CV result stored before 2026-09-15 carries this
bias.

Fix (2026-09-15): the dada2 and sintax files of a database hold the same
records in the same order (same counts, lengths, sequences and taxonomy
strings on the 7 databases and their `mini_*` files). `cross_val()` now names
the dada2 records with the sintax identifiers (`id_fasta`, `query_id_fasta`,
`record_ids()` checks the lengths record by record) and keeps the dada2
headers for the truth table and the training fasta. Smoke test on the mini
files: 366 usable dada2 queries of 400 drawn on `Unite_s_all_20250219_Fungi`.

## Validation plan before trusting a new production run

1. Build the versioned reference files: `download_reference_sources()` then
   `derive_all_variants()` (docs/reference_databases.md), and check the first
   header of each new file.
2. `run_project("assign_taxo_mini")`: check that `Genus_lca__*` columns now
   differ between databases, that sintax / lca reach a non-zero F1 at
   Kingdom and Phylum in `01_load_and_clean.qmd`, that `tar_meta()` shows no
   "NAs introduits" warning and that no `Genus_*` value contains `(` (B17).
3. `run_project("cross_val_mini")`: check that dada2 is no longer 0 and that
   lca ranks are plain rank names.
4. Only then run the production projects and regenerate every figure.
