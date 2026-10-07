# Manuscript framing — proposals to triage

Moved out of `ROADMAP.md` on 2026-09-30 (it was its section 2). These items
shape how the paper is written, not what the pipeline computes, so they are
worked through when the Results are drafted, not night by night. Nothing here
is applied: each item is a proposal for the developer to accept, amend or
refuse. An accepted framing decision is written in `docs/objectives_design.md`
§4 (as the next OD number) and summarised in `HISTORY.md` §3d.

Sources: the conceptual half of [`reviews/critique_2026-09-22.md`](../reviews/critique_2026-09-22.md),
read together with its **addendum** (the developer's framing of 2026-09-22) and
with [`reviews/literature_2026-09-22.md`](../reviews/literature_2026-09-22.md). The `→ §…`
pointers refer to sections of the critique. Ids are kept as they were in the
ROADMAP, since other files cite them.

Tags: `[Priority / facility]`, Priority = Critical → High → Medium → Low,
facility = easy (< 1 day) → moderate (a few days) → hard.

## 1. Framing and reporting (critique 2026-09-22)

- [ ] **C22-P1** `[Critical / easy]` Write the Introduction around 3–4 falsifiable predictions (P1a curation bias, P1b breadth not size, P1c abstention trade-off, P1d trimming no-op) instead of the "best trade-off" shopping question; the trade-off becomes the practical conclusion. The grid already tests all four at zero extra compute. → §1.1 C1
- [x] **C22-P2** `[High / easy]` Pre-specify **one primary endpoint** and label everything else secondary or exploratory in the manuscript; the critique proposes F1 and PPV at Genus and Species on the two mocks. **Supersedes decision 5 and narrows decision 31** — left open on 2026-09-22 (developer), to settle when the first real figures exist; the accepted text goes to `docs/objectives_design.md` §4. → §1.1 C2
  - Settled 2026-10-05 (developer): F1 and PPV at Genus and Species, the four co-primary endpoints (OD 37).
  - Revised 2026-10-06 (developer, OD 43): primary endpoints F1 and PPV at Genus. The paper says it benchmarks taxonomic assignment as a whole, not species assignment; Species F1 and PPV go to a species section (exact-match methods), which expects fewer species named than with the bootstrap classifiers.
- [ ] **C22-P3** `[High / easy]` Held-out confirmation: any *recommendation* is selected on Pauvert and confirmed on Tedersoo, so a parameter set is never recommended from the mock that picked it. 0 compute; narrows decision 31 (which keeps the sweeps untuned on purpose). → §1.1 C3
- [ ] **C22-P4** `[Medium / easy]` Sensitivity run with the truthless units left out of the matrix (the `truth_depth` machinery already does this for ties); report whether the Q2 ranking survives. 42 % of the Tedersoo ASVs carry no truth. → §1.1 C4
- [ ] **C22-P5** `[High / easy]` Promote two findings buried in prose to Results: the **dada2 collapse on references above 65 536 taxonomies**, a 32-bit overflow in dada2 itself, not breadth (`HISTORY.md` §2, "Databases building", corrected 2026-09-30) and the **`_cut` no-op** (`docs/experimental_design.md` §5.5). Demote **Q5 only** — the critique's "demote Q3" is withdrawn by its addendum, Q3 is a primary axis. → §1.2 C5, addendum
- [ ] **C22-P6** `[Critical / moderate]` **Required, not optional** (addendum, because Q3 is a headline claim): analyse Q2 as paired within-source contrasts and write "UNITE vs EUKARYOME", not "databases" (n = 2 providers, the rest is a nested filtering ladder); report the Q3 gain against the **effective** number of independent voters, not the column count. → §1.3 T1
- [ ] **C22-P7** `[Medium / moderate]` One-species-out novelty test: remove every record of one mock species from a reference and score its ASVs. Turns a closed-book benchmark into an open one, answers Orsholm et al. 2026's central question, and doubles as the Q3 stress test. ≤ 1 grid run, or a scoring variant of the CV. → §1.3 T2, `reviews/literature_2026-09-22.md` §4
- [ ] **C22-P8** `[High / easy]` Headline F1/PPV/TPR **on real units only**; the controls report through their own rates (`NA_fake`, `ext_fungi`, `ext_correct`, `ctrl_assigned`). A reporting choice, not a rewrite — the per-unit cells are already separate. Otherwise the absolute metrics depend on `prop_fake = 0.5` and on the control count, which are design constants. → §1.3 T3
- [ ] **C22-P9** `[Medium / easy]` State the Pauvert/Tedersoo confound (ITS region, primer pair, sequencing run and community all differ at once): restrict the claim to "two mocks differing in ITS region and protocol"; only the in silico set or a third mock can separate region from community. → §1.3 T4
- [ ] **C22-P10** `[Medium / easy]` Sensitivity column with strict submitted-name-only scoring: synonym tolerance interacts with the UNITE/EUKARYOME nomenclature difference that Q2 measures. One-line variant of `tc_metrics_unit()`, the truth table already carries `Genus_accepted` / `Species_accepted`. → §1.3 T5
- [ ] **C22-P11** `[High / moderate]` Bootstrap CIs on the **paired** differences (10⁴ resamples of the ~150 scored units, Δ recomputed per resample) and McNemar on the discordant units. Cheap — the confusion cells are per unit — and it turns every Q2 claim into a supported one instead of eye-judged point estimates. → §1.3 T6
- [ ] **C22-P12** `[High / moderate]` D1a: give the cross-validation one explicit estimand in the Methods ("consistency of each method with each database's own labels at fixed reference size"), add the negative controls to the folds (own TODO in `R/cross_val.R`, finally a CV-side TNR), and report the per-database novelty rate beside every CV figure. Only if the cross-validation is run (OD 28, `docs/parked.md`). → §1.3 T7
- [ ] **C22-P13** `[Medium / easy]` Keep M2 secondary: wall time and peak RAM as measured (after the final timing run, ROADMAP M2.4), CO₂eq only with a sensitivity note — if ±2× on the carbon intensity changes no ranking, say exactly that. → §1.3 T8
- [ ] **C22-T9** `[Low / easy]` One sentence in D1c saying what biological datasets can and cannot do: NA rates on real communities show that the mocks are not unrepresentative in their NA behaviour, they cannot support accuracy claims. → §1.3 T9
- [ ] **C22-N1** `[High / easy]` Declaration of interests: the harness and part of the benchmarked software share an author. Declare it, pin the package commits per run (the `session_info` targets already record them), and state that the engines (vsearch, BLAST+, dada2) are independent third-party tools — the wrapper can only flatter its own parameter plumbing, which `tests/test_assign_compute.R` audits against the direct calls. → addendum, `reviews/literature_2026-09-22.md` §4.4

## 2. Discussion notes

- [ ] **0.11c** `[Low / easy]` Discussion notes: mono-nucleotide shuffle (M1.6), set-membership bias of the former metrics made moot by the per-unit truth, parameters reported on the whole mock rather than tuned on a third (hleap §5 point 8).
- [ ] **M1.6** `[Low / easy]` Discussion note on the mono-nucleotide shuffle (Hleap used `esl-shuffle`, di-nucleotide preserving), which likely inflates TNR for the k-mer methods. Same as 0.11c.
- [ ] **ITSx origin** ITSx 1.1.3 calls a quarter of the fungal Pauvert ASVs (short ITS1 amplicons, mostly Hypocreales and Glomerellales) Tracheophyta: its putative origin is no fungal filter on short inputs (SA7, `HISTORY.md` §6.4).
- [ ] **Length cap (OD 35)** Methods: the full EUKARYOME is used minus its 7 871 records over 5 kb, and why: its long mitochondrial genomes drew dada2 to a genus for most shuffled controls; show the paired `fake_`-only test (`analysis/sandbox/check_length_cap_fakes.R`).
- [ ] **Tedersoo read retention** Methods: `ERR16773639` keeps 15.3 % of its read pairs against 53.7 % for `ERR16773638`, because its reverse primer (ITS4ngsUni) is found at the start of only 15 % of its R2 reads; both runs kept (developer, 2026-10-05), scoring being per unit. → `HISTORY.md` §2 D1c
- [ ] **Read filtering (OD 41)** Methods: default `filterAndTrim()` kept, citing Tosadori & Bosch 2026 (Web Ecology 26:47–59), who found no significant difference in reads retained between `maxEE = 2` / `truncQ = 2` and the Rolling et al. 2022 settings (`maxEE = 8` / `truncQ = 8`) on simulated data and on the same Pauvert mock, and "essentially identical results" (their Figs. S1, S2, 2). Our test agrees from the other end: dada2's own default is `maxEE = Inf` (1.40.0; their "default" 2 is the tutorial value), and `maxEE = 8` / `truncQ = 8` changed nothing on the units both settings keep (`analysis/sandbox/test_read_filtering.R`). State our `maxEE = Inf` exactly. (developer, 2026-10-06)
- [ ] **D1b scope** In silico is a positive control only, not a Q2 dataset (OD 29): its reads come from records of `Unite_s_all_20250219`, so any UNITE-vs-EUKARYOME contrast on it measures the leak. Poorly known lineages were excluded and phyla flattened by design when drawing its input — state it.

## 3. Findings flagged "FOR THE PAPER"

Measured and archived in `HISTORY.md` §2 ("Databases building"), each with its
numbers and the script that produced them:

- sintax and one Laccaria clone: a single over-long reference record turns a k-mer classifier into a false-positive machine (SA5, OD 32).
- Without ITSx, blastn is an order of magnitude slower for some markers (`docs/experimental_design.md` §4.1).
- stampa against the common LCA: the result that led to OD 33.
- dada2 above 65 536 distinct taxonomies (full EUKARYOME, every archived release, full `Unite_s_all`): a 32-bit overflow in `assignTaxonomy()` scores taxonomies with the wrong k-mer profile; patched, a full-breadth EUKARYOME gives 84 % correct genera (2026-09-30; supersedes "breadth dilutes").

Candidate added when this file was created (2026-09-30), not flagged by the
developer: sintax was not reproducible run to run without a seed (SA6,
`HISTORY.md` §6.2).
