# Objectives design: inputs, additions, databases, methods, parameters, metrics

Written 2026-09-17 at the developer's request, then **updated three times the same
day with the developer's answers** (section 4, decisions 1–12, 13–20, 21–25, their
text kept verbatim). For each
objective of the benchmark, this document lays out **which input data, which
additions to those inputs, which databases, which assignment methods, which
parameters, and above all which metrics actually answer the question**.

**Scope (developer, 2026-09-17).** The benchmark tests **one step: taxonomic
assignment of Fungi, and nothing else**. Anything about community ecology is
out of scope: no simulated community structure (miaSim / Hubbell dropped), and
the biological datasets are read only through the distribution of `NA_real`,
`NA_fake` and `ext_fungi` across studies.

Every statement about existing code was read from the source on 2026-09-17
(file and line given). Status labels used throughout:

- **built** — exists and runs today;
- **exists, not wired** — the object or file exists but no pipeline uses it for
  this purpose;
- **decided** — settled by the developer (section 4), not implemented;
- **proposed** — a design choice still open (section 5).

Costs are not repeated here: see `docs/experimental_design.md` (measured cost,
datasets and their ITS region). Where a decision changes the number of
computations, the change is stated. **No code has been changed for any of the
decisions below.**

---

## 1. The six axes

### 1.1 Input data

| Code | Object (store) | Dataset | Region | Size | Truth | Status |
|---|---|---|---|---|---|---|
| `P.asv` | `d_asv` (`store_dada2`) | Pauvert mock | ITS1 | 196 ASVs (195 assigned, decision 20) | `taxo_mock.csv`: **189 strain rows**, 181 species, 97 genera, 67 families, 30 orders, 11 classes, 3 phyla | built, assigned |
| `P.otu` | `d_vs` (`store_dada2`, vsearch 97 %) | Pauvert mock | ITS1 | 143 OTUs | same | **decided: assigned directly** (today Q5 derives OTU taxonomy from member ASVs) |
| `P.itsx` | `d_asv_itsx` (`store_assign_taxo`) | Pauvert mock | ITS1 extracted | 195 | same | built, assigned |
| `T.asv` | `d_asv` (`store_dada2_tedersoo_illumina`) | Tedersoo Illumina mock | full ITS | 238 ASVs | **not fetched** (paper Table S1 / UNITE) | built, not assigned |
| `T.otu` | `d_vs` (`pipelines/dada2_bio.R:300`) | Tedersoo Illumina mock | full ITS | 180 OTUs (142 / 122 after MUMU) | same | **decided: assigned directly** |
| `T.itsx` | `d_asv_itsx` | Tedersoo Illumina mock | full ITS extracted | — | same | decided (ITSx everywhere); needs a per-dataset `itsx_region` (`experimental_design.md` §5.1) |
| `S.iss` | D1b: InSilicoSeq reads → DADA2 | in silico | ITS | 100 input records, 3 seeds | the input records (per-sequence) | fastq generated, no targets project. **Source decided: `Unite_s_all_20250219`** (with singletons), which is what the code already reads (`config.R:79`) |
| `B.*` | loaders not written | biological (GSSP air, Korhonen wood ×2 = ITS2; GloSED soil = full ITS; leaf undecided) | ITS2 / full | unmeasured | **none** | proposed |
| `CV` | records of the benchmarked database | cross-validation | ITS1 amplicon (queries trimmed ITS1F–ITS2) | 5 000 queries × folds | the database labels (per-sequence) | built |

**Dropped (scope decision):** `S.mia`, the miaSim / Hubbell branch of
`analysis/in_silico_simulation.qmd` (D1b-B). With it goes the `S.iss` − `S.mia`
contrast; D1b now rests on `S.iss` alone.

### 1.2 Additions to the inputs (negative controls)

| Code | Function | What is added | Correct answer |
|---|---|---|---|
| `+fake` | `comparpq::add_shuffle_seq_pq()`, prefix `fake_` | shuffled copies of the input's own sequences, `prop_fake = 0.5` (`config.R:15`): 97 on Pauvert, 119 on Tedersoo | **NA at every rank**, on every database |
| `+ext` | `comparpq::add_external_seq_pq()`, prefix `external_` | 100 **non-Fungi** records of `Unite_s_all_20250219` (`fake_ref_asv_100.fasta`, `derive_fake_ref()`) | **depends on the database** (decision 6) |

Status: built for `P.*` and `T.asv`; **decided for every input**, including the
OTU inputs (decision 4) and the biological datasets (implied by the scope
decision: `NA_fake` and `ext_fungi` need them).

**Scoring of `+ext`, decided (decision 6):**

| Database simplification | `external_` enter TN? | `external_` scored by |
|---|---|---|
| **all** | no | `ext_fungi`, `ext_correct` |
| **fungi** | **yes** (NA is the correct answer) | TN, and `ext_fungi` |
| **fungi + rep** | no (decisions 15, 22) | `ext_fungi`; `ext_correct` **at Kingdom only** |

`fake_` controls enter TN on every database.

What `ext_correct` is for (decision 15): **measuring how well a database ×
method combination lets non-fungal sequences be filtered out**, not how well
it names them.

**Consequence to keep in mind.** The TN denominator now depends on the
database: `fake_` only on "all" and "fungi + rep", `fake_` + `external_` on
"fungi". TNR stays
comparable across databases (it is a rate), but **ACC and MCC are not**, since
their margins change size. Q2 comparisons between simplification levels should
therefore rest on PPV, TPR, F1, `NA_real` and the control rates, not on ACC /
MCC.

On the ITSx input the `fake_` sequences are shuffles of the extracted region
(`d_asv_itsx` → `d_asv_itsx_shuffled`, CLAUDE.md), so they are shorter than on
the raw input; the `external_` sequences are identical in both inputs.

### 1.3 Databases

Decided grid: `[EUK | UNITE | UNITE_sh] × [all | fungi | fungi + rep]`, with
**UNITE_sh = UNITE with singletons** (decision 1) and **no cutadapt-trimmed
database** (decisions 24, 25).

| Source ↓ / Simplification → | **all** | **fungi** (kingdom exactly Fungi) | **fungi + rep** |
|---|---|---|---|
| **EUK** = `EUK_ITS_v2.1` | built, 1 877 003 records | built, 1 076 870 | **decided**, not built: + 53 844 reps (5 %) from `EUK_ITS_v2.1` |
| **UNITE** = `Unite_all_20250219` | built, 156 820 | built, 102 137 | **decided**, not built: + 5 107 reps from `Unite_all_20250219` |
| **UNITE_sh** = `Unite_s_all_20250219` | built, 266 589 | built, 168 030 | **decided**, not built: + 8 402 reps from `Unite_s_all_20250219` |

#### Trimmed databases: a final zoom, not a grid level

**No `_cut` database in the benchmark for now** (decisions 24, 25). Trimming
(`fungi + rep cut` = `_Fungi_cut` records + untrimmed reps, decision 16) will
be tried **only on the few best database × method combinations**, once the
main grid has identified them, and on the Pauvert mock only. This follows the
project's strategy: answer secondary questions with **targeted zooms** rather
than by multiplying the total computation.

`EUK_ITS_v2.1_Fungi_cut` stays on disk until then, and `pipelines/cross_val.R`
keeps its `_cut` query rule (`query_fasta`, `query_db` columns), which the zoom
will need.

> **Note for the zoom (decision 25).** A `fungi + rep cut` database mixes
> trimmed fungal records (median ≈ 210 bp) with untrimmed reps (median
> above 500 bp). Global-alignment identity in vsearch ignores terminal gaps
> and blastn's cover is computed on the query, so lca and blastn should not be
> biased by the length gap, but the k-mer methods (sintax, dada2) are
> untested. Comparing `ext_fungi` between `fungi + rep` and `fungi + rep cut`
> measures it directly.

#### `fungi + rep` construction

Reps are non-fungal records (decisions 3, 13, 14, 21):

- **number**: 5 % of the fungal records of the database;
- **source**: the "all" release of the same source (`EUK_ITS_v2.1`,
  `Unite_all_20250219`, `Unite_s_all_20250219`);
- **groups**: kingdoms with **≥ 10 records** for UNITE and UNITE_sh, **≥ 100
  records** for EUKARYOME, after removing the non-group labels (`cf.*`,
  `*.reg*`, `_mitochondrion`, `_Bacteria`, `_Archaea`, `_plastid`,
  `_pseudogene`, spikes, `Eukaryota_kgd_Incertae_sedis`; also `Tartumycota`,
  1 EUKARYOME record carrying a fungal phylum name in the kingdom field);
- **share**: equal across kingdoms when possible, computed as a water-filling
  rule — every kingdom gets the same cap, a smaller kingdom gives all its
  records and the shortfall is spread over the others;
- **exclusion**: the 100 `external_` control records are never drawn as reps.

Kingdom lists and allocations (sintax files, counted 2026-09-17; the
`external_` records are not subtracted, at most 100 records overall):

**UNITE (`Unite_all_20250219`)** — 15 kingdoms ≥ 10 records, target 5 107
(cap 832; 10 kingdoms give all their records). Dropped: Apusozoa (5),
Choanoflagellozoa (2), Planomonada (2).

| Kingdom | Records | Reps |
|---|---:|---:|
| Viridiplantae | 35 585 | 832 |
| Metazoa | 8 466 | 832 |
| Alveolata | 5 489 | 832 |
| Stramenopila | 1 781 | 831 |
| Rhizaria | 1 056 | 831 |
| Protista | 287 | 287 |
| Rhodoplantae | 211 | 211 |
| Amoebozoa | 126 | 126 |
| Cryptista | 95 | 95 |
| Euglenozoa | 62 | 62 |
| Parabasalia | 48 | 48 |
| Haptista | 36 | 36 |
| Heterolobosa | 34 | 34 |
| Ichthyosporia | 34 | 34 |
| Glaucocystoplantae | 16 | 16 |

**UNITE_sh (`Unite_s_all_20250219`)** — 16 kingdoms ≥ 10 records, target
8 402 (cap 1 306; 11 kingdoms give all their records). Dropped: Picozoa (7),
Choanoflagellozoa (6), Planomonada (4), Nucleariae (2).

| Kingdom | Records | Reps |
|---|---:|---:|
| Viridiplantae | 58 965 | 1 306 |
| Metazoa | 16 123 | 1 306 |
| Alveolata | 8 548 | 1 306 |
| Stramenopila | 3 194 | 1 305 |
| Rhizaria | 2 211 | 1 305 |
| Amoebozoa | 422 | 422 |
| Rhodoplantae | 363 | 363 |
| Protista | 336 | 336 |
| Cryptista | 251 | 251 |
| Euglenozoa | 172 | 172 |
| Heterolobosa | 85 | 85 |
| Parabasalia | 76 | 76 |
| Ichthyosporia | 61 | 61 |
| Haptista | 54 | 54 |
| Glaucocystoplantae | 43 | 43 |
| Apusozoa | 11 | 11 |

**EUK (`EUK_ITS_v2.1`)** — 39 kingdoms ≥ 100 records, target 53 844 (cap
2 631; 25 kingdoms give all their records). Dropped: 19 named kingdoms of
5–72 records (538 records in total).

| Kingdom | Records | Reps |
|---|---:|---:|
| Alveolata | 157 388 | 2 631 |
| Viridiplantae | 156 269 | 2 631 |
| Metazoa | 124 008 | 2 631 |
| Straminipila | 103 401 | 2 631 |
| Rhizaria | 95 960 | 2 631 |
| Amoebozoa | 78 224 | 2 631 |
| Euglenozoa | 23 548 | 2 631 |
| Choanoflagellozoa | 7 185 | 2 631 |
| Heterolobosa | 4 720 | 2 631 |
| Haptista | 4 232 | 2 631 |
| Nucleariae | 3 829 | 2 631 |
| Centroheliozoa | 3 270 | 2 631 |
| Ichthyosporia | 3 242 | 2 631 |
| Cryptista | 3 009 | 2 631 |
| Apusomonada | 1 695 | 1 695 |
| Breviatae | 1 470 | 1 470 |
| Nereolimi | 1 376 | 1 376 |
| Divergimonada | 1 301 | 1 301 |
| Rhodoplantae | 1 241 | 1 241 |
| Picomonada | 1 236 | 1 236 |
| Telonemae | 968 | 968 |
| Tartarolimi | 872 | 872 |
| Preaxostyla | 804 | 804 |
| Hemimastigophora | 728 | 728 |
| Collodictyozoa | 648 | 648 |
| Corallochytriozoa | 607 | 607 |
| Malawimonada | 443 | 443 |
| Aphanistiae | 435 | 435 |
| Jakobae | 418 | 418 |
| Brechocaryoniae | 382 | 382 |
| Planomonada | 360 | 360 |
| Anhoristiae | 344 | 344 |
| Effugaciae | 339 | 339 |
| Fornicata | 336 | 336 |
| Aquilomonada | 302 | 302 |
| Parabasalia | 203 | 203 |
| Edaphomonada | 181 | 181 |
| Filasteriae | 177 | 177 |
| Obamonada | 144 | 144 |

The two releases do not always use the same kingdom name for the same group
(`Stramenopila` / `Straminipila`; EUKARYOME holds both `Apusomonada` and
`Apusozoa`): harmless for the reps, since each database draws from its own
release, but it matters for the `external_` filter below.

Exclusion of the `external_` records from the reps: they come from
`Unite_s_all_20250219`, so for UNITE the exclusion can go by accession, but
**EUKARYOME records carry other identifiers: excluding them there needs a
sequence match**, not a name match.

**Built 2026-09-17** (`make_databases.R::derive_fungi_rep()`, ROADMAP 0.2).
Totals as planned: 5 107, 8 402 and 53 844 representatives from 15, 16 and 39
kingdoms. The exclusion is done by identical sequence for every release, which
removes 58 records of `Unite_all_20250219`, 100 of `Unite_s_all_20250219` and 13
of `EUK_ITS_v2.1` from the candidates **before** water-filling, so a few
allocations differ from the tables above: UNITE cap 833 (Rhizaria 832,
Rhodoplantae 210, Amoebozoa 124, Parabasalia 47, Heterolobosa 33,
Ichthyosporia 33); UNITE_sh cap 1 308 (Viridiplantae 1 309, Amoebozoa 417,
Rhodoplantae 362, Cryptista 249, Euglenozoa 171, Heterolobosa 84, Parabasalia
75, Ichthyosporia 60, Haptista 53); EUKARYOME unchanged. The drawn records are
listed in `data/data_raw/refseq/sources/<db>_reps.csv`.

Identical sequence is not the only way a control can sit among the reps, and
the rule was tightened on 2026-09-22 (developer): **a representative must stay
below 97 % identity with every `external_` control** (`config.R::rep_max_identity_to_external`,
vsearch). The rule applies to the representatives only. The whole releases keep
the controls — they are records of `fake_ref_source`, which is what makes
`ext_correct` computable there — and the `_Fungi` databases hold none.

Controls with a hit at ≥ 97 % identity in each database, before the fix
(vsearch, 100 controls):

| Database | ≥ 97 % | ≥ 99 % | = 100 % |
|---|---:|---:|---:|
| `Unite_s_all_20250219` | 100 | 100 | 100 |
| `Unite_all_20250219` | 58 | 58 | 58 |
| `EUK_ITS_v2.1` | 78 | 68 | 44 |
| the three `_Fungi` and `EUK_ITS_v2.1_Fungi_cut` | 0 | 0 | 0 |
| `Unite_s_all_20250219_Fungi_rep` | 6 | 0 | 0 |
| `Unite_all_20250219_Fungi_rep` | 7 | 0 | 0 |
| `EUK_ITS_v2.1_Fungi_rep` | 10 | 5 | 3 |

After the fix (databases rebuilt 2026-09-22, 47, 47 and 403 candidates
excluded): **0** control with a hit at >= 97 % in the three `_Fungi_rep`
databases; the other lines are unchanged.

#### Consequence for the `external_` controls (decision 21)

The `external_` controls must not come from a kingdom removed from the reps.
The current file (`data/data_raw/fake_ref/fake_ref_asv_100.fasta`, 100
records of `Unite_s_all_20250219`) breaks that rule for **9 records**:
`Eukaryota_kgd_Incertae_sedis` (4), and one each of Planomonada, Picozoa,
Nucleariae, Choanoflagellozoa (below 10 records in UNITE_sh) and Apusozoa
(11 in UNITE_sh, but 5 in UNITE). The other 91: Viridiplantae 30, Metazoa 22,
Stramenopila 12, Alveolata 7, Rhizaria 5, Amoebozoa 5, Cryptista 2, and one
each of Rhodoplantae, Protista, Parabasalia, Ichthyosporia, Heterolobosa,
Haptista, Glaucocystoplantae, Euglenozoa.

**Done 2026-09-17** (decision 26: kingdoms retained in **all three** lists).
`derive_fake_ref(kingdoms = external_control_kingdoms())` draws from 13
kingdoms: Alveolata, Amoebozoa, Cryptista, Euglenozoa, Haptista, Heterolobosa,
Ichthyosporia, Metazoa, Parabasalia, Rhizaria, Rhodoplantae, Stramenopila
(= EUKARYOME's Straminipila, `config.R::kingdom_aliases`), Viridiplantae. The
new file holds 100 records from 65 phyla (Viridiplantae 38, Metazoa 24,
Stramenopila 11, Alveolata 9, Amoebozoa 5, Rhizaria 5, Cryptista 2, one each
of the others). Regenerating it **invalidates every assignment of the store**
(file target `fake_ref_file`, CLAUDE.md), at no extra cost since the new grid
reruns everything.

#### Computations per dataset

4 methods; ITSx input on the non-full databases, i.e. fungi and fungi + rep
(ROADMAP S8.5); OTU input (`d_vs`, no MUMU, no ITSx, decision 19) on every
database. Same count for every dataset, now that no dataset has a cut level.

| Input | Databases | Computations |
|---|---:|---:|
| ASV | 9 | 36 |
| OTU (`d_vs`) | 9 | 36 |
| ITSx ASV | 6 | 24 |
| **Total** | | **96**, against 44 today |

The lca change of §1.5 (lca rerun once) is included, since the new grid reruns
every computation anyway.

**Cross-validation on fewer databases** (decision 24): the full 3 × 3 grid
would be 4 methods × 9 databases × 2 variants = 72 targets (56 today); the
subset is open (Q-l). Only fungal records are queried, on every simplification
(decision 20).

### 1.4 Methods

One computation per method × database × input serves all parameter rows
(`R/assign_compute.R:22-57`):

| Method | Computation (called once) | Parameter rows derived afterwards |
|---|---|---|
| `dada2` | `dada2::assignTaxonomy(minBoot = 0, outputBootstraps = TRUE)` | threshold on bootstraps |
| `sintax` | `MiscMetabar::assign_sintax(min_bootstrap = 0, clean_pq = FALSE)`, default `cmd_args = "--sintax_random"` | threshold on bootstraps |
| `lca` | `MiscMetabar::add_new_taxonomy_pq(method = "lca", clean_pq = FALSE)` | **copies** today; **decided: `lca_cutoff` sweep** (§1.5) |
| `blastn` | `MiscMetabar::blast_pq(score_filter = FALSE)` (raw hits) | score filters + vote via `assign_blastn(blast_table = )` |

### 1.5 Parameters

| Method | Sweep | Fixed (value, source) | Compute cost |
|---|---|---|---|
| `dada2` | `min_bootstrap` ∈ {0.4, 0.5, 0.6} (built) | other `assignTaxonomy()` arguments at their defaults | free |
| `sintax` | `min_bootstrap` ∈ {0.4, 0.5, 0.6} (built) | `--sintax_random` | free |
| `lca` | **`lca_cutoff` ∈ {0.8, 0.9, 1}** (decision 9) | `id = 0.5`, `top_hits_only = TRUE`, `maxaccepts = 0`, `maxrejects = 32` (`MiscMetabar/R/vsearch.R:1655-1659`) | **free after one change** — see below |
| `blastn` | `vote_algorithm` ∈ {rel_majority, abs_majority, unanimity} (built); **`min_id` ∈ {90, 92, 95, 97}** (decisions 10, 17) | `min_cover = 80` (`config.R`, decision 18); `min_bit_score = 50`, `min_e_value = 1e-30` (`MiscMetabar/R/blast.R:969-972`) | free: filters apply to the stored raw hits |

**lca: how the cutoff sweep can be free.** `assign_vsearch_lca()` already runs
one vsearch command that writes two files (`MiscMetabar/R/vsearch.R:1700-1726`):
`--lcaout` (the LCA vsearch computed with `--lca_cutoff`) and `--userout` with
`query+id+target` (every retained hit with its taxonomy header). With
`top_hits_only = TRUE` and `maxaccepts = 0` (unlimited), `--userout` holds
exactly the top hits the LCA is computed on, so the three cutoffs can be
derived in R from one stored hit table, as S8.2 does for blastn.

Checked with vsearch on 20 queries against
`mini_Unite_s_all_20250219_Fungi`: `--top_hits_only` also filters `--userout`
(20 rows with it, 137 254 without; no row below its query's best identity).
The test queries came from the database itself, so each had a single 100 %
hit: ties between several top hits were not exercised.

**`--lcaout` alone cannot give the other cutoffs** (decision 18 suggests
using it): it holds one lineage per query, already truncated at the cutoff of
the run. A lineage cut at 1 cannot be extended to what 0.8 would have kept.
The derivation needs the hit table; `--lcaout` serves as the equivalence check.
Proposed MiscMetabar shape, mirroring `blast_pq()` / `assign_blastn(blast_table = )`:
`assign_vsearch_lca()` returns or saves the `--userout` table, and accepts it
back through a `hits_table =` argument to apply any `lca_cutoff` without
running vsearch.

Three facts bound this:

- today the hit table goes to `tempdir()` and is discarded: the lca computation
  must change to keep it, which **reruns lca once per database × input**;
- `resolve_vector_ranks(nb_agree_threshold = )`, which the `vote_algorithm`
  branch uses, takes an **absolute count**, not the fraction `--lca_cutoff`
  expects (vsearch help: "fraction of matching hits required for LCA"), so the
  derivation needs its own fraction rule;
- equivalence is checkable: the R derivation at cutoff 1 must reproduce the
  stored `--lcaout` columns.

The derivation lives **in MiscMetabar** (decision 18), in the `hits_table` shape above (decision 23).

**Parameter rows per computation** after these decisions: dada2 3, sintax 3,
lca 3 (no longer copies), blastn 3 votes × 4 `min_id` = **12**, i.e. 21 rows
per method × database × input (12 today).

### 1.6 Metrics

#### What `comparpq::tc_metrics_mock_vec()` computes today

Per rank and per assignment column (`comparpq/R/compare_taxo.R:161-241`):

| Metric | Unit | Definition as implemented |
|---|---|---|
| TP | real sequences | real sequences whose value at the rank **belongs to the set** of truth values |
| FP | real sequences | real sequences with a non-NA value **not in** the truth set; **controls excluded** (line 190: `tax_table[!fake_taxa_cond, …]`) |
| FN | **truth-table rows** | truth rows whose value **never appears** in the column |
| TN | control sequences | controls (`fake_` + `external_`) left NA |
| PPV (= TAR) | ratio | TP / (TP + FP); TAR is an alias built in chapter 01 |
| FDR | ratio | FP / (FP + TP) |
| TPR (= TDR) | ratio | TP / (TP + FN); TDR is an alias built in chapter 01 |
| TNR | ratio | **TN / (TN + FN)** (line 209) — non-standard |
| F1_score | ratio | 2 TP / (2 TP + FP + FN) |
| ACC | ratio | (TP + TN) / (TP + TN + FP + FN) |
| MCC | ratio | standard formula; 0 by convention when a margin is zero |

Properties that decide what these metrics can answer:

1. **TP is set membership, not per-sequence identity.** An ASV of species A
   assigned to species B, when B is also in the mock, counts as a TP; the mock
   carries no per-ASV truth. **To discuss in the manuscript** (developer,
   decision 11): a limited bias, because mock species are chosen to be distinct
   enough to be assigned; mitigated by the abundance-filtered metrics below.
2. **An assignment on a control is counted nowhere** — not an FP, only a missing
   TN.
3. **F1, ACC and MCC mix three units**: sequences (TP, FP), truth rows (FN) and
   controls (TN). Comparisons hold within one input at fixed size; absolute
   values are not comparable across datasets nor across units (ASV vs OTU).
4. **TNR mixes controls and truth rows.**

#### FP and TNR: to discuss (decision 7, Q-a)

Decision 7 is to redefine TNR "as TN / (TN + FP)". **Verified in comparpq:
with the current FP this would not be a specificity**, because FP counts real
sequences only (property 2): the ratio would divide controls by real
sequences.

The developer's reference point (Q-a): Hleap et al. 2021, *"FP are sequences
that were assigned to a species not present in the community"* (quoted by the
developer). **Superseded by `docs/hleap_2021_metrics.md`**, which reads the paper and the `TA_pipes` code: the two options below predate it.
Read literally, that definition **includes the controls**: a `fake_` sequence
given any species is "assigned to a species not present in the community".
comparpq departs from it by excluding the controls from FP. Two coherent
options follow:

| | **A — Hleap FP** | **B — FP on real sequences only (today)** |
|---|---|---|
| FP | real sequences outside the truth set **+ scored controls given a value** (`ctrl_assigned`) | real sequences outside the truth set |
| TNR | TN / (TN + FP) — the formula of decision 7 | TN / (TN + `ctrl_assigned`) = TN / n_scored_controls |
| PPV, F1, MCC | **fall when controls are assigned**: over-classification of noise is penalised in the headline metrics | blind to controls; `ctrl_assigned` reported beside them |
| Property 2 | fixed | still true |
| Caveat | F1 / PPV now depend on the number of controls injected (`prop_fake`, 100 `external_`), a design constant: fine within a dataset, one more reason not to compare absolute values across datasets. TNR's denominator mixes real misassignments and controls, so it is a rate over "sequences that should not have received that name", not over controls | the specificity is clean, but a method that names every shuffled sequence keeps a high F1 |

In both options, "scored controls" follows decision 6: `fake_` everywhere,
`external_` only on "fungi" databases. On "all" and "fungi + rep" databases an
`external_` given its true lineage must be **neither TN nor FP** (it is not
"not present in the community" in any useful sense: it is a correct answer),
otherwise option A would count correct answers as errors.

My recommendation is **A**: it matches the reference the manuscript cites,
it makes F1 sensitive to over-classification (which is the failure the
controls exist to detect), and the caveat is already covered by comparing
within a dataset. FN stays a count of truth rows in both options (property 3).

Also in comparpq: `tests/testthat/test-compare_taxo.R:159` asserts
`TN + FN == 0` when every control is assigned, i.e. it pins the current
denominator and would change with the fix. The vignette
`benchmark-da-methods.Rmd:878` already uses the textbook
`Specificity = TN / (TN + FP)` in another context. **Nothing edited in
comparpq**; the proposal is in section 5 (Q-a).

**Implemented 2026-09-22** (ROADMAP 0.6 and 0.7): the per-unit truth in
`R/mock_truth.R` and the scoring in `comparpq::tc_metrics_unit()`, which
replaces the set-membership counting below for the mocks. The chapters still
call `tc_metrics_mock()` until ROADMAP 0.8.

#### Metrics decided (decision 7)

| Code | Definition | Needs | Available on |
|---|---|---|---|
| `NA_real` | share of **real** sequences left NA at the rank | nothing new | every input, including `B.*` |
| `NA_fake` | share of `fake_` left NA (should be 1) | `+fake` | every input |
| `ext_fungi` | share of `external_` assigned Kingdom = Fungi: **a false positive on every database** | `+ext` | every input |
| `ext_correct` | share of `external_` given their true non-fungal lineage | `+ext`, truth parsed from their UNITE headers | "all" databases at every rank; "fungi + rep" (and cut) **at Kingdom only** (decision 15) |
| `ctrl_assigned` | controls given a non-NA value, split `fake` / `ext` (= FP_ctrl) | controls | every input |
| `misassign_seq` | share of sequences assigned to a wrong taxon, per sequence | **per-sequence truth** | `S.iss` and `CV` only (in CV it is `wrong_classifications`); **not computable on the mocks** |
| TNR (redefined) | option A or B of the table above (Q-a) | controls | every input with controls |

Dropped from the list: `edgar_oc` / `edgar_mc` were not retained (ROADMAP M1.2
stays open).

#### Abundance-filtered metrics (decision 11)

The same mock metrics recomputed after discarding the units (ASVs or OTUs)
with **20 reads or fewer** ("> 20" kept, as the developer wrote). Free: a
post-hoc filter on stored columns. Controls are not filtered (their abundance
is a copy convention, not a signal). Not applicable to presence / absence
datasets.

#### Cross-validation and cost metrics (unchanged)

CV (`R/cross_val.R:406-420`): `good_classifications`, `wrong_classifications`,
`prop_NA` split into `NA_classif` / `NA_database`, per rank, mean ± sd over
folds — the only per-sequence error rate on real reference sequences.

Cost (`with_autometric()`): `wall_time_s`, `peak_resident_mb`, `mean_cpu_pct`
per computation, shared by its parameter rows; CO₂eq through
`greenAlgoR::ga_footprint()` (M2.2, open).

---

## 2. The objectives table

"default" = `min_bootstrap = 0.5` for dada2 / sintax, `rel_majority` for blastn
(`filter_default_settings()`, `analysis/00_setup.R:70-78`); lca's default
becomes `lca_cutoff = 1` (today's value). **Bold** metrics answer the question.

| Objective (question) | Inputs | Additions | Databases | Methods | Parameters | Metrics that answer it |
|---|---|---|---|---|---|---|
| **Q0 — Central trade-off**: which method × database × parameter gives the best quality for its cost? | `P.asv`, `T.asv`; `CV` for robustness | `+fake`, `+ext` | all | all 4 | all sweeps | **F1 and TPR / PPV at Genus / Species against `wall_time_s` and `peak_resident_mb`** (Pareto front), **`NA_real`**, **TNR (redefined)**, `ext_fungi` |
| **Q1 — Algorithm and parameters** | `P.asv`, `T.asv` (agreement between the two mocks is itself a result) | `+fake`, `+ext` | all (database as random effect, Q1.7) | all 4 | dada2 / sintax bootstrap; **lca `lca_cutoff`**; blastn vote × **`min_id`** | **PPV and TPR along each sweep**, **`NA_real`**, **`ctrl_assigned`**; F1 per rank; abundance-filtered variants |
| **Q2 — Database source and simplification** | `P.asv`, `T.asv`; `CV` | `+fake`, **`+ext` essential** | **3 sources × 3 simplifications** (no cut, decision 24) | all 4 | default | **paired Δ PPV / TPR / F1** (all → fungi, fungi → fungi + rep), **`ext_fungi`**, **`ext_correct`** (all), `NA_real`, TNR; not ACC / MCC (§1.2); cost per database (measurement to redo) |
| **Q3 — Consensus voting** | `P.asv`, `T.asv` | `+fake`, `+ext` | reuses Q1 / Q2 columns | reuses | default | **F1 of each consensus vs the best single column**, **`NA_real`**, `ctrl_assigned`; cost = sum of member computations |
| **Q4 — ITSx preprocessing** | `P.asv` vs `P.itsx`; `T.asv` vs `T.itsx` | `+fake`, `+ext` (same names and seed in both inputs) | fungi, fungi + rep | all 4 | default; blastn `min_cover` 80 vs 95 on the ITSx input | **paired Δ F1 raw vs ITSx**, **blastn `NA_real`**, `ext_fungi` (length-independent, unlike `fake_`-based TN); share of ASVs undetected by ITSx |
| **Q5 — ASV vs OTU** | `P.asv` vs `P.otu`; `T.asv` vs `T.otu` — **`d_vs` assigned directly** | `+fake`, `+ext` on each input | all | all 4 | default | **`NA_real` per unit**, **`ext_fungi`**, **TNR (redefined)** — all unit-independent; F1 only with a unit-consistent definition (property 3); abundance-filtered variants |
| **D1a — Cross-validation** | `CV` | none (ROADMAP TODO: fake sequences in folds for a TN) | **subset of the databases** (decision 24, Q-l); fungal records queried only (decision 20) | all 4 | default; `remove_tested` TRUE / FALSE | **`good_classifications`, `wrong_classifications` (= `misassign_seq`), `NA_classif`** per rank, standard vs leaked |
| **D1b — In silico** | `S.iss` only (source `Unite_s_all_20250219`) | `+fake`, `+ext` | all | all 4 | default | **`misassign_seq`** (per-sequence truth), **`NA_real`**, F1 |
| **D1c — Biological communities** | `B.*` | **`+fake`, `+ext`** (decided through the scope) | all, **restricted by region** (ITS2 datasets) | all 4 | default | **distribution of `NA_real`, `NA_fake`, `ext_fungi` across studies**, set against the mocks' values for the same method × database; no TP / FP / FN; no ecological reanalysis |

### Notes per objective

**Q0.** Cost is per computation and shared by its parameter rows, so along a
threshold sweep only quality moves. The Pareto front is drawn over method ×
database × input, with parameters as points on the same cost. The lca change
(§1.5) keeps that property for lca.

**Q1.** Raising a threshold trades FP for NA; reading PPV, TPR and `NA_real`
together shows whether it removes errors or only answers.

**Q2.** The fungi filter has one specific risk — non-fungal sequences called
fungal once their relatives are removed — and only `ext_fungi` measures it;
`fungi + rep` is the attempt to fix it. dada2 on the full EUKARYOME already
gives a Q2 result: 0 % genus, because its bootstrap collapses on a broad
reference (ROADMAP "Databases building").

**Q5.** Assigning `d_vs` directly removes the derivation from member ASVs but
not the unit problem: TP / FP are still counted on 196 ASVs vs 143 OTUs and FN
on truth rows. The unit-independent rates carry the answer.

**D1b — leakage, developer's position (decision 5).** The drawn records stay in
the references, but in silico sequencing modifies them, unlike CV, so D1b is
not treated as an upper bound. One measurable check keeps that claim honest:
the **share of `S.iss` ASVs identical to their source record** after DADA2
(denoising corrects most substitution errors). If it is high, D1b behaves like
the leaked CV variant for those ASVs and the discussion should say so.

**D1c.** On real data a higher `NA_real` can be a correct abstention (the taxon
is absent from the database) or a failure; only the controls and the
comparison with the mocks for the same method × database separate the two.

---

## 3. Coverage: which input serves which objective

● main input, ○ secondary.

| Objective | `P.asv` | `P.otu` | `P.itsx` | `T.asv` | `T.otu` | `T.itsx` | `S.iss` | `B.*` | `CV` |
|---|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| Q0 trade-off | ● | | | ● | | | | | ○ |
| Q1 algorithm × parameters | ● | | | ● | | | ○ | | ○ |
| Q2 databases | ● | | | ● | | | | ○ | ● |
| Q3 consensus | ● | | | ● | | | | | |
| Q4 ITSx | ● | | ● | ● | | ● | | | |
| Q5 ASV vs OTU | ● | ● | | ● | ● | | | | |
| D1a cross-validation | | | | | | | | | ● |
| D1b in silico | | | | | | | ● | | |
| D1c biological | ○ | | | ○ | | | | ● | |

---

## 4. Decisions settled 2026-09-17

The developer's answers, verbatim after **>>**, with what each changes.

1. **"UNITE_sh"** — read as UNITE *with* singletons (`Unite_s_all_20250219`),
   so UNITE = `Unite_all_20250219`.
   >> **I Confirm**.

   → §1.3.

2. **`_Fungi_cut`** — drop it from Q2, keep it for the Pauvert ITS1 mock only,
   or keep it as a fourth simplification level?
   >> **I deliberatly remove _Fungi_cut because it is very close to wat _Fungi_rep_cut will do**. High computation price for a low rewards.

   → §1.3: removed from the grid as a level of its own; `_Fungi_rep_cut` is
   defined in decision 16.

3. **`fungi + rep`** — which groups, how many records per group, from which
   release; confirm the exclusion of the 100 `external_` control records.
   >> I confirm the exclusion of the 100 external_ control records. I think we may take 5% of the number of fungal sequences as a good total number with a equal repartition (when possible) accross groups. For eukaryome take from eukaryome_all, for Unite take from Unite_all (therefore without singletons).

   → §1.3: sizes settled; sources revised by decision 14, groups by decision 13.

4. **OTU inputs (Q5)** — assign `d_vs` directly, or keep deriving the OTU
   taxonomy from member ASVs?
   >> We must assign d_vs directly instead of deriving from ASV. Your are write d_vs instead of d_otu is clear.

   **Widened 2026-09-22 (developer), while building ROADMAP 0.8a:** the
   directly assigned OTUs are the input of the grid, as decided above, **and**
   the taxonomy derived from the member ASVs is kept beside them, at no
   computational cost. Q5 then answers a question neither alone can: "cluster
   then assign" against "assign then cluster". In `analysis/01_load_and_clean.qmd`
   and in `res_comp_tax` the two are the `unit` values `OTU` and `OTU_derived`
   (the ASVs being `ASV`). `resolve_otu_taxonomy()` and the `otu_clusters`
   target stay in service for the second one.

   → §1.1, §1.3 (+36 computations per dataset), Q5.

5. **In silico source and leakage.**
   >> In silico source must be Unite_s_all (with singletons), in all cases the idea is that some error will occur during in silico "sequencing" . So it is not a true upper bound, yes at the start only sequences from the references databes will be selected, but in opposition to CV (cross-validation), in silico DO modify the sequences in order to make assignation less easy.

   → §1.1, D1b note (one check proposed).

6. **Scoring of `external_` controls on "all" databases.**
   >> external controls must be taken into account differently regarding the database. With "all" database, no TN is compute but the score `ext_fungi` / `ext_correct` is compute. With _fungi database, it is used for computing TN.

   → §1.2 table; "fungi + rep" in decisions 15 and 22.

7. **Metric set and TNR.**
   >> Yes we must compute NA_fake, NA_real, ext_fungi, ext_correct, ctrl_assigned and misassign_seq. Your are right, we must modify TNR as TN / (TN + FP). Also Verify in the package comparpq.

   → §1.6. Verified in comparpq: the FP of the formula must be the FP **on
   controls**, see Q-a.

8. **Controls on biological data.** Not answered directly; settled by the scope
   decision (12), since `NA_fake` and `ext_fungi` need controls. The abundance
   convention question disappears: every control metric counts sequences, not
   reads.

9. **lca parameters.**
   >> For lca method, we must try to make lca_cutoff variable c(0.8, 0.9, 1), for this, we may need to modify the MiscMetabar::assign_vsearch_lca function in order to use the cutoff after a first computation ? As lca is the second longest method, new computation is a real problem.

   → §1.5: feasible from the `--userout` hit table, one lca rerun to keep it
   (location: decision 18).

10. **blastn `min_id`.**
    >> blastn's min_id must be swept at no compute cost apply ROADMAP S1.3.

    → §1.5; values in decision 17.

11. **Set-membership TP.**
    >> We must take into memory for the discussion that "TP is set membership, not per-sequence identity. An ASV of species A assigned to species B, when B is also in the mock, counts as a TP. The mock carries no per-ASV truth, so misassignment within the mock is invisible." It is true but not a big deal because species in the mock communities are selected to be sufficiently different to be assigned. May be a good approach to mitigate this bias is to compute complementary metrics after discarding ASV/OTUs with a low number of sequences (>20).

    → §1.6 property 1 and abundance-filtered metrics.

12. **Scope of the project.**
    >> Lets clarify the delimitation of the project : we want to test for a very specific steps -> taxonomic_assignation on Fungi AND NOTHING ELSE. So, all part that are taking about ecology are not interesting here. We must remove mia and hubbel distribution and in analyse of real community we will focus only on the distribution of NA_real, NA_fake and ext_fungi across studies.

    → header, §1.1 (`S.mia` dropped), D1b, D1c.


### Second round (answers to the first "Still open" list)

13. **Groups for `fungi + rep`** (was Q-b).
    >> We must only select the most important kingdoms, make me the list of kingdom for
    UNITE and EUKARYOME, but certainly, we will removed `cf.*`, `*.reg*`, `_mitochondrion`, `_Bacteria`, `_Archaea`, `_plastid`, `_pseudogene`, spikes, Eukaryota_kgd_Incertae_sedis.

    → §1.3: lists with counts and allocations; which kingdoms are "most
    important" is settled by decision 21.

14. **UNITE_sh reps** (was Q-b').
    >> Yes we must take the rep from the corresponding _all database, so
    take UNITE_sh_fungi_rep from UNITE_sh_all.

    → §1.3. This **revises decision 3**: each database draws its reps from its
    own "all" release (UNITE_sh from `Unite_s_all_20250219`).

15. **`external_` on `fungi + rep`** (was Q-c).
    >> ext_correct of fungi + rep must only be compute at kingdom level. But 
    in all case, the idea of ext_correct is to see how non-fungal species 
    can be filtered out by the database x assignation methods. 

    → §1.2, §1.6. TN on `fungi + rep`: decision 22.

16. **`_Fungi_rep_cut`** (was Q-d).
    >>_Fungi_rep_cut is the merging of a Fungi_cut (not use in the analysis) and 
    the uncut representative sequences. The main is to see if removing useless 
    nucleotide for a given marker can decrease time/cpu and increase the quality
    of the asignation.  

    → Superseded for the grid by decisions 24 and 25: the cut is a final zoom.

17. **blastn `min_id` values** (was Q-e).
    >> {90, 92, 95, 97}

    → §1.5: 12 blastn rows per computation.

18. **lca cutoff derivation** (was Q-f).
    >> cutoff derivation must lives in MiscMetabar. The ideal is to use
    vsearch's `--lcaout` and then use the output (as we did in assign_blastn previously) 

    → §1.5, with one correction: `--lcaout` cannot give the other cutoffs, the
    `--userout` hit table can.

19. **OTU input** (was Q-g).
    >> No mumu in this project. No ITSx on d_vs

    → §1.3 computation counts.

20. **`fungi + rep` in cross-validation** (was Q-h).
    >> only fungi in cross-validation

    → §1.3, D1a.


### Third round

21. **Kingdoms for `fungi + rep`** (was Q-b).
    >> All kingdoms with ≥ 10 records for Unite and Unite_sh and All kingdoms with ≥ 100 records
    for Eukaryome. However, we must not take ext_sequences from the removed sequences, so
    we must apply this filter when selecting external_ sequences.

    → §1.3: lists and allocations; the current `external_` file breaks the
    rule for 9 records and must be regenerated; which list filters them is
    decision 26.

22. **TN on `fungi + rep`** (was Q-i).
    >> I confirm. 

    → §1.2: `external_` do not enter TN on `fungi + rep`.

23. **lca** (in conversation, answering the `hits_table` proposal of §1.5).
    >> Ok pour le lca.

    → §1.5.

24. **Databases and cross-validation** (in conversation).
    >> Concernant les bases je pense qu'il serait pertinent d'enlever le cut par cutadapt de toute les bases de données. On pourra voir à la fin sur les 2 meilleurs bases (une unite et une eukaryome) si ça apporte quelque choses et uniquement sur Pauvert. On adopte une stratégie ou pour certaines questions on préfère faire des mini-zoom plutôt que multiplier le nombre de calcul total. Pour la validation croisée, c'est pareil on peux sans doute faire une validation croisée sur moins de database.

    → §1.3: no cut level; cross-validation subset in Q-l. Supersedes the grid
    part of decision 16.

25. **Cut level** (was Q-j and Q-k).
    >> No cut for the moment. Forget the idea of cut until we found the few best 
    databases x methods. We will try the cut on thoses combinaison only.

    >> No cut for the moment. Forget the idea of cut until we found the few best 
    databases x methods. We will try the cut on thoses combinaison only. Just 
    remember (with a notes) that there is a question of different length if we 
    cut the "fungi + rep" database. 

    → §1.3: zoom on the best database × method combinations (Pauvert only,
    decision 24), with the length note.

26. **Kingdom list for the `external_` controls** (was Q-m, answered in
    conversation): "Intersection des 3" — only kingdoms retained in UNITE,
    UNITE_sh and EUKARYOME.

    → §1.3: implemented and file regenerated 2026-09-17
    (`config.R::rep_kingdom_min_records`, `rep_kingdom_exclude`,
    `kingdom_aliases`; `make_databases.R::external_control_kingdoms()`;
    tests in `tests/test_make_databases_helpers.R`).

27. **Databases in cross-validation** (was Q-l, answered in conversation):
    decided later, from the results on the mock communities.

    → §1.3; the CV grid stays as it is until then.

---

## 5. Still open

- ~~**Q-a — FP, FN, TN and TNR.**~~ **Settled 2026-09-22** in
  `docs/hleap_2021_metrics.md` §5–6 and implemented as ROADMAP 0.6 and 0.7:
  every unit falls in exactly one cell against its **own** truth, matched to
  the Sanger sequence of a strain (`R/mock_truth.R`); a unit is scored down to
  its `truth_depth` and leaves the matrix below it; TNR is kept, as
  TN / (TN + FP) with every FP, controls and real misassignments alike; the
  `external_` controls enter the matrix on the Fungi databases only. The
  options of §1.6 predate that reading and are superseded by it.

  One departure from Hleap is deliberate and belongs in the Discussion: the
  parameters here are reported on the whole mock, not tuned on one third of it
  (ROADMAP 0.11c).

## 6. Files that now contradict these decisions (not edited)

- ~~`config.R::benchmark_dbs`~~, ~~`pipelines/assign_taxo.R` /
  `R/values_map.R`~~: resolved 2026-09-17 (ROADMAP 0.1–0.5).
- `pipelines/cross_val.R`: runs on `config.R::cv_db_list`, the 7 former
  databases, on purpose until decision 27 is taken.
- ~~`analysis/in_silico_simulation.qmd` (D1b-B branch), `CONTEXT.md` glossary
  entry, `ROADMAP.md` D1b-B section~~: resolved 2026-09-22 (ROADMAP 0.9), the
  miaSim branch is removed and its code lives in the git history.
- ~~`analysis/01_load_and_clean.qmd`~~: resolved 2026-09-22 (ROADMAP 0.8a). It
  now scores the directly assigned OTUs (decision 4); the taxonomy derived from
  the member ASVs is **kept beside it**, on the developer's decision of the
  same day, so Q5 compares "cluster then assign" with "assign then cluster"
  (`unit` = `OTU` against `OTU_derived`).
- ~~`docs/experimental_design.md`~~: resolved 2026-09-22 (ROADMAP 0.10) except
  its §3, whose measured costs were read on the 44-computation grid and are
  kept as an order of magnitude until a production run on the new grid is
  measured. Its header says so.
- ~~`CLAUDE.md`~~: resolved 2026-09-17 and kept current since.
