# Hleap et al. 2021 metrics: what the paper says, what the code does, what comparpq does

Written 2026-09-17 at the developer's request ("je veux être certain de faire
comme lui"), to settle the FP / TNR question of `docs/objectives_design.md`
(Q-a) on primary sources.

Sources read:

- **Paper**: Hleap J.S., Littlefair J.E., Steinke D., Hebert P.D.N., Cristescu
  M.E. (2021) *Assessment of current taxonomic assignment strategies for
  metabarcoding eukaryotes*, Molecular Ecology Resources. PDF from the
  developer's Zotero library (`5KEAUMA8`), sections "Inclusion of true
  negatives", "Metrics for benchmarking", "Curating mock community
  composition".
- **Code**: <https://github.com/jshleap/TA_pipes>, commit `adc8876`
  (2021-02-03). The metrics are computed by `score()` in
  `optimize_n_score.py` (lines 211–274), called by `benchmark_run.py` through
  `opt_n_score(query, truthfile)` for every method, mock and taxonomic level.
- **comparpq**: `tc_metrics_mock_vec()`, `comparpq/R/compare_taxo.R:161-241`.

---

## 1. The paper's definitions (verbatim)

> "The TP are sequences known to be derived from a species present in the mock
> community and assigned to this taxon. FP are sequences that were assigned to
> a species not present in the community, while FN are species known to be
> present but that had no sequence assigned to them. Finally, the TN are
> sequences present in the dataset that should not be assigned to any species.
> To compute TN, we determined how many of the randomized sequences were
> assigned to any of the 7 taxonomic levels and subtracted that value from the
> total randomized sequences in each mock community."

Composite metrics (paper): FDR = FP / (FP + TP), TPR = TP / (TP + FN),
PPV = TP / (TP + FP), F1 = 2 · TPR · PPV / (TPR + PPV), MCC (standard formula).
**There is no TNR / specificity in Hleap et al.**

Negative controls (paper): 10 % of the denoised, chimera-free sequences,
shuffled with `esl-shuffle` (HMMER 3.2.1) **keeping mono- and di-nucleotide
composition**. No external (non-target) sequences.

Per-sequence truth (paper, "Curating mock community composition"): every ASV
got its own label by a BLAST search against the reference sequences of the
mock members (≥ 99.5 % identity and 100 % query cover), then phylogenetic
placement for the rest, with manual curation; unidentified sequences were
removed (`Mocks/unidentified.tar.gz`). The truth files (`Mocks/*_realized.txt`)
hold **one lineage per ASV** (e.g. `sq10;size=231;` → *Notropis atherinoides*).
ASVs with fewer than 8 reads were removed beforehand.

## 2. What the code actually computes

`score()`, per taxonomic level, on the outer merge of truth and prediction by
sequence name (`mer`), with `shuffled` = the shuffled sequences:

```python
tp = mer[~shuffled & (mer[tr] == mer[pr])].shape[0]
fp = mer[(~shuffled & (mer[tr] != mer[pr]) & (~mer[pr].isnull())) |
         (shuffled & ~mer[pr].isnull())].shape[0]
fn = mer[~shuffled & (mer[tr] != mer[pr]) & (mer[pr].isnull())].shape[0]
tn = nshuf - mer[shuffled & ~mer[pr].isnull()].shape[0]
```

and asserts `tp + fp + tn + fn == all_calls`. Each sequence falls in exactly
**one** cell, at each level:

| Sequence | Prediction at the level | Cell |
|---|---|---|
| real | equal to **its own** truth label | TP |
| real | another taxon | FP |
| real | NA (or absent from the output) | FN |
| shuffled | any taxon | **FP** |
| shuffled | NA | TN |

MCC: standard formula, `-inf` when a margin is zero. F1 as in the paper.

### Two places where the code differs from the paper text

1. **FN.** Paper: *species* present with no sequence assigned. Code: *real
   sequences left unassigned* at that level. The code is what produced the
   published numbers (the `tp + fp + tn + fn == all_calls` assertion only holds
   with sequence counts).
2. **TN.** Paper: shuffled sequences not assigned at *any* of the 7 levels.
   Code: shuffled sequences not assigned *at the scored level*. Identical at
   Kingdom (a method that assigns a lower rank also fills Kingdom), different
   at lower ranks: in the code a shuffled sequence named only to Phylum is a TN
   at Genus.

Minor: line 222 (`truth.copy().loc[...] = 'shuffled'`) writes into a copy and
has no effect; harmless, the shuffled labels are already set at load time
(line 93).

### Scores are computed on a held-out part of the mock

Paper ("Optimization of parameters and evaluation of accuracy") and code
agree. For each method and mock, a grid search picks the parameters that
**maximise MCC** (F1 for HMMUFOTU; ties broken toward the least stringent set,
e.g. lowest `p_id` for BLAST) on **one third** of the mock sequences, drawn at
random (`split_train()`, lines 601–619, redrawn until both parts hold shuffled
sequences). The reported metrics are computed **on the other two thirds**,
with those parameters, at every level (`optimize_it()`, lines 736–746). The
search is run at the Species level (`main(..., taxlevel = 6)`).

The grids are wide: BLAST `p_id` 70–100 by 5, six e-values, six
`max_target_seqs`; LCA adds `m_hit`, `p_hit`, `n_hit`; IDTAXA `thresh`,
`bootst`, `mind`.

So Hleap reports **tuned** methods evaluated out of sample, not a fixed
parameter grid shown in full. Our design (every parameter row reported on the
whole mock) answers a different question — how sensitive each method is to
its parameters — and is closer to Hleap's own R² analysis of the training
folds than to his headline figures.

## 3. comparpq against Hleap

| | Hleap (code) | `comparpq::tc_metrics_mock_vec()` | Same? |
|---|---|---|---|
| Unit | sequences, everywhere | TP / FP: real sequences; FN: truth-table rows; TN: controls | **no** |
| Truth | one lineage **per sequence** | a **set** of values per rank (all strains of the mock) | **no** |
| TP | real sequence = its own label | real sequence whose value is **in the set** | **no** |
| FP, real sequences | real sequence given another taxon | real sequence given a value **outside the set** | **no** (misassignment inside the mock is invisible) |
| FP, controls | **shuffled sequence given a taxon counts as FP** | controls excluded from FP | **no** |
| FN | real sequence left NA | truth row whose value never appears in the column | **no** |
| TN | shuffled left NA at the level | controls (`fake_` + `external_`) left NA at the rank | yes for `fake_`; `external_` is our extension |
| TNR | **not defined** | TN / (TN + FN) | not in Hleap |
| FDR, TPR, PPV, F1 | as paper | same formulas | formulas yes, inputs differ |
| MCC zero margin | `-inf` | 0 | no (convention only) |
| Shuffle | di-nucleotide preserving (`esl-shuffle`), 10 % of the sequences | `sample()` on each sequence (mono-nucleotide only), `prop_fake = 0.5` | **no** |
| Low-abundance filter | ASVs < 8 reads removed before assignment | none (decision 11 proposes > 20 reads for complementary metrics) | no |

**Conclusion: comparpq follows neither the code nor the text of Hleap.** The
closest reading of comparpq is the paper's FN sentence (taxa never found), but
counted on strain rows, not species.

## 4. What "doing it like Hleap" requires here

1. **A per-ASV truth for each mock.** This is the prerequisite: without it, TP /
   FP cannot be per sequence. Hleap built it by BLAST (≥ 99.5 %, 100 % cover)
   against the reference sequences of the mock members. For the Pauvert mock,
   `taxo_mock.csv` names Sanger sequences per strain (`MockStrain` and `ITS1`
   columns, e.g. `Agaricus_essettei_M187_sangerITS1`), but **no Sanger fasta was
   found in the project** (searched 2026-09-17). For Tedersoo, Table S1 is not
   fetched yet. This also changes decision 11 (set-membership TP was accepted as
   a limited bias) and makes `misassign_seq` computable on the mocks
   (= FP on real sequences / real sequences).
>> We must take the same approach. Tedersoo table S1 do contain the sequences (accesible at
  https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2F1755-0998.70189&file=men70189-sup-0002-TableS1.xlsx) and
  Pauvert table S1 (List of the fungal strains of the mock community and their ITS Sanger
sequences.) but I can't found it on the web. And can't we also use the mock_hleap2021 
fastq files ?

2. **Scoring as in `score()`**: one cell per sequence per rank; shuffled
   sequences assigned → FP; real sequences left NA → FN.
3. **No TNR in the headline metrics** (Hleap has none). If kept, TN / (TN +
   shuffled assigned) is the specificity on controls; it is our addition.
4. **`external_` controls** do not exist in Hleap. Decision 6 (TN only on
   "fungi" databases, `ext_fungi` / `ext_correct` elsewhere) is an extension to
   state as such in the Methods.
5. **Shuffle method**: `add_shuffle_seq_pq()` shuffles single nucleotides; Hleap
   preserves di-nucleotide composition. A mono-nucleotide shuffle destroys more
   k-mer structure, so shuffled sequences are easier to reject for k-mer methods
   (sintax, dada2): TN is likely higher than with Hleap's controls.
>> Add a note for the discussion
6. **Filter and proportion**: < 8 reads removed and 10 % shuffled in Hleap,
   against no filter and 50 % here. The proportion changes FP and TN counts,
   hence PPV, F1 and MCC.
>> We will easily compute metrics with filter at c(2:10, 12, 15, 20, 30, 40, 50, 100) reads and
see how it change the metrics (almost free cost). Good idea, add it to roadmap

7. **Scored level for TN**: follow the code (per level) or the text (any
   level); they only differ below Kingdom.   
8. **Parameters**: Hleap tunes on one third of the mock (max MCC) and reports
   the other two thirds. Doing the same means a split of each mock and a
   "best parameters" selection before the headline figures; our sweeps
   (bootstrap, `lca_cutoff`, blastn `min_id` × vote) already provide the grid,
   at no extra compute.
>> I prefer my method than Hleap, we have different purpose. 

Points 1 and 5–8 are choices for the developer; nothing is changed in comparpq
or in the project.

## 5. Decisions (2026-09-17)

The developer's answers are kept above after **>>**; what they settle:

| Point of §4 | Decision | Recorded in |
|---|---|---|
| 1. Per-ASV truth | **Same approach as Hleap.** Pauvert: Table S1 (189 strains, Sanger ITS sequences) stored as `data/data_raw/metadata/pauvert2019_table_s1_mock_sanger.csv`. Tedersoo: Table S1 (103 specimens, Sanger reads) downloaded by the developer, stored as `data/data_raw/metadata/tedersoo_mock_table_s1.xlsx` and `tedersoo_mock_table_s1_sanger.csv` | ROADMAP M1.4 |
| 2. Scoring | Detailed in §6; Q2a–Q2d decided | ROADMAP M1.4 |
| 3. TNR | **Kept** in the metric set | — |
| 5. Shuffle method | Note for the discussion | ROADMAP M1.6 |
| 6. Low-abundance filter | Sweep of *n* ∈ {2, …, 10, 12, 15, 20, 30, 40, 50, 100} reads, post-hoc, **on the mocks and the biological datasets only** (not on cross-validation nor in silico) | ROADMAP M1.5 |
| 6. Share of shuffled sequences | **Kept** at 50 % (`config.R::prop_fake`) | — |
| 7. TN level | **Per rank, as in the code**: a shuffled sequence named only to Phylum is FP at Phylum and TN at Genus | — |
| 8. Parameters | **Keep our design** (every parameter row reported on the whole mock): the purpose differs from Hleap's tuned, held-out scores | — |
| `mock_hleap2021` | Folder **deleted** 2026-09-17 (developer) | ROADMAP |

### Where the Pauvert Table S1 came from

The ScienceDirect supplementary file `1-s2.0-S1754504818302800-mmc2.csv`,
named in the `ADDENDUM` of the authors' code deposit
(<https://doi.org/10.15454/VKTWKR>, which explains how they built
`mock_sanger.fasta` from it), downloaded from
`https://ars.els-cdn.com/content/image/1-s2.0-S1754504818302800-mmc2.csv`
(md5 `78797995ca2058d395886403da090d3f`). Columns: `Mock_Strain`,
`Unique_ITS1_Sanger_Sequence` (175 Yes), `Present_In_Raw_MiSeq_Data` (160
Yes), `Putative_Cause_Absence`, `Phylum` … `Species`, `Sanger_Sequence`
(189 non-empty, median 542 bp). The raw reads of the same study are the
Dataverse deposit <https://doi.org/10.15454/8CVWRR> (`ITSEnz.Dik{A,B,C}`).

### Why `mock_hleap2021` cannot replace them

`data/data_raw/mock_hleap2021/Mocks/` holds Hleap's own mocks, already
denoised (`*_realized.fa` + per-ASV truth `.txt`), not fastq. They are **COI**
mocks of animals (fish, insects, zooplankton: Leray fragment and
MLepF1/LepR1), so they fall outside the scope of decision 12 (fungal
taxonomic assignment only) and cannot be assigned against UNITE or EUKARYOME
ITS. Their use would only be to check that our scoring code reproduces
Hleap's numbers on his own data, which would also need his predictions.

The folder was deleted on 2026-09-17 (developer); nothing in the pipelines
read it.

### Tedersoo Table S1

`men70189-sup-0002-tables1.xlsx`, supplementary Table S1 of
doi:10.1111/1755-0998.70189, downloaded by hand (Wiley answers HTTP 403 to
scripts), md5 `d9a4fcb8b295534c43ad586902f7bd8a`. One sheet: 103 specimens
(`TUF…` vouchers), `Identification` (**a binomial only, no lineage**),
`Sanger read` (lower case, 466–953 bp, median 599; 101 distinct sequences),
then read counts per taxon for 24 platform × pipeline × PCR-cycle columns.
11 specimens have 0 reads in both Illumina ASV columns. The csv copy
upper-cases the sequences and names the count columns
`<platform>__<pipeline>__PCR_<n>_cyc` (the last four PacBio columns repeat
the names of the four before them in the sheet and get a `__rep1` suffix).

---

## 6. Point 2 in detail: scoring like `score()` here

### The rule, per rank

Every unit that is scored (ASV or OTU of the mock, and control) falls into
exactly one cell at each rank:

| Unit | Prediction at the rank | Cell |
|---|---|---|
| real, with a truth at this rank | equal to **its own** truth | TP |
| real, with a truth at this rank | another name | FP |
| real, with a truth at this rank | NA | FN |
| real, **no truth at all** (no Sanger match, Q2b) | any name / NA | FP / FN |
| real, truth stopped above this rank (tie, Q2c) | — | left out at this rank |
| `fake_` | any name | FP |
| `fake_` | NA | TN |
| `external_`, "fungi" database (decision 6) | any name / NA | FP / TN |
| `external_`, "all" or "fungi + rep" database | not in the confusion matrix | scored by `ext_fungi`, `ext_correct` |

So TP + FP + FN + TN = number of scored units at that rank, and every
metric (PPV, TPR, F1, MCC, TNR) is computed on that one unit.

Worked example, one Pauvert ASV of *Botrytis cinerea* assigned
*Botrytis pseudocinerea*: today it is a TP at Species (the name belongs to the
mock); with this rule it is an **FP at Species** and a TP at Genus. A
`fake_` sequence assigned down to *Ascomycota* only: today a missed TN at
every rank; with this rule an FP at Kingdom and Phylum and a TN from Class
down.

### What changes with respect to today

- **FN** no longer counts truth-table rows never found, but real units left
  unassigned. A strain absent from the reads (29 of the 189 Pauvert strains
  are marked `Present_In_Raw_MiSeq_Data = No`) no longer lowers TPR: recovering
  the mock composition is a denoising question, outside decision 12.
- **Misassignment inside the mock becomes visible** (decision 11 becomes
  moot) and `misassign_seq` = FP on real units / real units.
- **TNR is kept** (point 3) as TN / (TN + FP) with every FP, controls and
  real misassignments together (decision 7, confirmed as Q2d).
- The OTU input needs its own truth (the centroid sequence of each OTU,
  matched the same way).
- `comparpq::tc_metrics_mock_vec()` takes a set of truth values per rank; the
  new rule needs **a truth per taxon** (a table indexed by taxon name).

### How many units get a truth (measured 2026-09-17)

`vsearch --usearch_global --id 0.9 --top_hits_only` of the ASVs of `d_asv`
against the Sanger sequences (identity as vsearch's default, terminal gaps
excluded; qcov = query cover). **Corrected the same day**: a first count read
the identities with `awk` under a French locale, which parsed `99.5` as `99`
and undercounted every threshold (it gave 138 / 98 and 75 / 0 at 99.5 %).

| Best identity ≥ | Pauvert (196 ASVs): truth | of which qcov = 100 % | tied on several Sanger | Tedersoo (238 ASVs): truth | of which qcov = 100 % | tied on several Sanger |
|---|---:|---:|---:|---:|---:|---:|
| 98 % | 180 (92 %) | 128 | 17 | 124 (52 %) | 0 | 4 |
| 98.5 % | 178 (91 %) | 127 | 16 | 124 (52 %) | 0 | 4 |
| 99 % | 172 (88 %) | 123 | 14 | 122 (51 %) | 0 | 4 |
| 99.5 % (Hleap) | 150 (77 %) | **105** | 11 | 115 (48 %) | **0** | 4 |

- 193 Pauvert ASVs and only **133 Tedersoo ASVs** have any hit at ≥ 90 %: 105
  Tedersoo ASVs match no Sanger read at all. They are discussed in the
  Tedersoo paper and are not investigated here (developer, 2026-09-17).
- **Tied** = the best identity is reached by several Sanger sequences at once,
  because some strains have identical ITS (six *Botrytis* strains in
  Pauvert). In every tie measured, the tied strains share the genus and differ
  at species: the truth stops at Genus (decision Q2c).
- Median query cover of the hits: 100 % for Pauvert, **76 %** for Tedersoo
  (none above 92 %): the Tedersoo ASVs (full ITS plus 18S and 28S flanks, 749
  bp median) are longer than the Sanger reads, so Hleap's 100 % query cover can
  never be reached there. Hleap's own queries (313 bp COI) sat inside full
  barcodes, where 100 % cover is natural.

### Decisions and open points

- **Q2a — Match rule.** Developer: a threshold per mock was meant, but it does
  not work for Tedersoo; proposed instead: **the Sanger sequence must be
  included in the ASV, with 0.5 % of errors allowed**. Generalised to both
  mocks as "the shorter of the two sequences is covered, identity ≥ 99.5 %",
  because the inclusion goes opposite ways: Pauvert ASVs (222 bp median) sit
  inside the Sanger reads (542 bp), Tedersoo Sanger reads (584 bp) sit inside
  the ASVs (750 bp). Measured below. **Confirmed by the developer
  (2026-09-17)**: a unit gets the truth of a Sanger sequence when identity ≥
  99.5 % **and** the shorter of the two sequences is covered at ≥ 99.5 %
  (vsearch `id`, `qcov`, `tcov`); ties take the lowest common rank (Q2c).
- **Q2b — Units without a truth: kept** in `NA_real` **and** in the confusion
  matrix (developer): the aim is a **relative** comparison of database × method
  couples, not Hleap's absolute scores. As in Hleap's code (`NaN != NaN`), such
  a unit is FN when left NA and FP when given any name, at every rank.
- **Q2c — Units matching several strains: truth = lowest common rank** of the
  tied strains (developer), and the unit is left out of the matrix at the
  ranks below.
- **Q2d — TNR = TN / (TN + FP)**, all FP, as written in decision 7
  (developer).
- **Names.** The Tedersoo truth is a binomial only: **the lineage (Kingdom →
  Genus) must be rebuilt** (developer, 2026-09-17; tool to choose: `taxize`
  0.10.1, `rgbif` 3.8.5 and the pqverse `taxinfo` are installed). Both truths
  must also use the names of the reference databases, since a synonym turns a
  correct answer into an FP. This already applies to today's set-based truth.
  Species names of the mocks found in the reference databases (exact match,
  2026-09-17):

  | Mock | Species | in `Unite_s_all_20250219` | in `EUK_ITS_v2.1` |
  |---|---:|---:|---:|
  | Pauvert | 181 | 134 | 148 |
  | Tedersoo | 103 | 71 | 94 |

  EUKARYOME stores the **epithet only** at species rank (`s:rufus` under
  `g:Lactarius`): the comparison above joins genus and epithet, as chapter 01
  does with its `Gen_sp_*` columns.

#### Names measured with the GNA verifier (2026-09-17)

Developer's proposal: link the mock species names to a classification with
`taxinfo::gna_verifier_pq()`, for Tedersoo and Pauvert alike. Two facts about
the tools:

- `gna_verifier_pq()` returns the **accepted name** (`currentCanonicalSimple`,
  synonyms resolved) but **no lineage**. It calls `taxize::gna_verifier()`,
  which drops the classification fields; the GNA verifier API itself does
  return them (`classificationPath`, `classificationRanks`, e.g.
  `Fungi|Basidiomycota|Agaricomycetes|Hymenochaetales|Rickenellaceae|Sidera|Sidera americana`
  for *Sidera americana* in the GBIF backbone).
- The lineage is that of the **accepted** name.

Results of the API (`POST /api/v1/verifications`), 181 Pauvert and 103
Tedersoo names:

| Source | Mock | Exact | Partial / fuzzy | Synonyms | With lineage |
|---|---|---:|---:|---:|---:|
| GBIF backbone (11) | Pauvert | 177 | 4 | 37 | 181 |
| GBIF backbone (11) | Tedersoo | 101 | 2 | 12 | 103 |
| Index Fungarum (5) | Pauvert | 176 | 5 | 35 | 165 |
| Index Fungarum (5) | Tedersoo | 100 | 3 | 13 | 101 |

The GBIF accepted name sits in **another genus** for 28 Pauvert and 10
Tedersoo names (e.g. *Botrytis calthae* → *Botryotinia*, *Agrocybe erebia* →
*Cyclocybe*). The Pauvert Table S1 lineage agrees with the GBIF lineage for
181 / 181 phyla, 173 / 181 classes, 163 / 180 orders, 149 / 180 families and
152 / 180 genera; some disagreements are wrong partial matches to check by
hand (*Colletotrichum cingulata* placed in Dothideomycetes, the genus-only name
*Fusarium* in Leotiomycetes).

**Which name the reference databases use** (exact species match, genus +
epithet for EUKARYOME):

| Mock | Database | Submitted name | GBIF accepted name | Either |
|---|---|---:|---:|---:|
| Pauvert (181) | `Unite_s_all_20250219` | 134 | 138 | **153** |
| Pauvert (181) | `EUK_ITS_v2.1` | 148 | 159 | **163** |
| Tedersoo (103) | `Unite_s_all_20250219` | 71 | 76 | **77** |
| Tedersoo (103) | `EUK_ITS_v2.1` | 94 | 95 | **96** |

Neither name alone is the databases' name: replacing the mock names by the
accepted names would fix some mismatches and create others. **Decided
(developer, 2026-09-17)**: keep the submitted name as the truth, **accept
either the submitted or the accepted name as correct** at Species and Genus,
and take the ranks above Genus from the GBIF lineage, for Tedersoo and for
Pauvert (checking the Table S1 disagreements above by hand).

The lineage comes from the new `classification_col = TRUE` option of
`taxinfo::gna_verifier_pq()` (added 2026-09-17, with `classification_ranks`;
columns `classificationPath`, `classificationRanks`,
`classificationKingdom` … `classificationGenus`). Checked on the 103 Tedersoo
names with `data_sources = 11`: 103 lineages, all in Fungi.

#### Q2a measured: containment rule (2026-09-17)

`vsearch --usearch_global --id 0.9 --maxaccepts 0 --maxrejects 0`, every hit
kept, fields `qcov` (share of the ASV aligned) and `tcov` (share of the Sanger
read aligned). Identity is vsearch's default (terminal gaps excluded).

| Rule | Pauvert (196 ASVs) | ties | Tedersoo (238 ASVs) | ties |
|---|---:|---:|---:|---:|
| id ≥ 99.5 %, ASV inside Sanger (qcov = 100 %) | 105 | 11 | 0 | 0 |
| id ≥ 99.5 %, Sanger inside ASV (tcov = 100 %) | 0 | 0 | 89 | 3 |
| id ≥ 99.5 %, shorter one covered = 100 % | 105 | 11 | 89 | 3 |
| **id ≥ 99.5 %, shorter one covered ≥ 99.5 %** | **149** | 11–12 | **99** | 3 |
| id ≥ 99.5 %, shorter one covered ≥ 98 % | 150 | 12 | 99 | 3 |
| id ≥ 99 %, shorter one covered = 100 % | 123 | 14 | 95 | 3 |
| id ≥ 99 %, shorter one covered ≥ 98 % | 172 | 16 | 105 | 3 |

- A strict 100 % cover loses 44 Pauvert ASVs whose shorter sequence is
  covered at 99.4–99.9 %: one or two terminal bases outside the Sanger read.
  Allowing the 0.5 % on the cover as well as on the identity recovers them.
- Ties stop at Genus in every case (no tie crosses a genus).

Reads behind each class, with the rule in bold above:

| | ASVs | % of reads | median reads per ASV |
|---|---:|---:|---:|
| Pauvert: truth | 149 | 95.6 | 276 |
| Pauvert: hit ≥ 90 % but no truth | 44 | 3.9 | 55.5 |
| Pauvert: no hit ≥ 90 % | 3 | 0.6 | 8 |
| Tedersoo: truth | 99 | 87.6 | 95 |
| Tedersoo: hit ≥ 90 % but no truth | 34 | 5.1 | 4 |
| Tedersoo: no hit ≥ 90 % | 105 | 7.3 | 3 |

The units without a truth are mostly rare in Tedersoo (median 3–4 reads), so
the abundance-filter sweep (M1.5) will show how much of the Tedersoo result
they drive (Q2b keeps them in the matrix).
