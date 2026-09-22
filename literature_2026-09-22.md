# Literature positioning — `benchmar_assign_taxo`

Companion to `critique_2026-09-22.md` (critic T5 "why another benchmark?").
Written 2026-09-22 from Crossref / Europe PMC / PubMed records and the papers'
abstracts (search date 2026-09-22). Short by design: one entry per paper, a
comparison matrix, then where this project sits. Full texts were not re-read;
summaries come from abstracts and records, so verify wording before quoting in
the manuscript.

---

## 1. The six papers this project must position against

**Hleap, Littlefair, Steinke, Hebert & Cristescu (2021)** — *Assessment of
current taxonomic assignment strategies for metabarcoding eukaryotes*,
Molecular Ecology Resources. **The direct ancestor.** Four COI animal mock
communities (fish, insects, zooplankton), per-sequence truth curated by BLAST
+ phylogenetic placement, shuffled sequences as true negatives, and one
program per assignment family: RDP-NBC, IDTAXA (k-mer), BLAST/LCA, PROTAX
(probabilistic), HMMUFOTU (phylogenetic). Methods are **tuned on ⅓ of the
mock and scored on ⅔**. Findings: no method dominates; accuracy falls with
reference divergence; over-classification of novel taxa is the dominant error.
This project already reads its paper *and* its code (`docs/hleap_2021_metrics.md`)
and adopts its per-unit scoring — but studies **fungal ITS**, adds the database
axis and reports parameter sweeps instead of tuned optima.

**Orsholm, Zito, Somervuo, …, Roslin & Furneaux (2026)** — *Discovering the
unseen: a performance comparison of taxonomic classification methods for
unknown DNA barcodes*, Methods in Ecology and Evolution 17:2574–2593
(doi:10.1111/2041-210X.70358). **The current state of the art, and the paper
that owns the question our design does not ask.** COI (arthropods) and ITS
(fungi), with train/test sets deliberately separated by taxonomy so that
"novel" taxa occur at every rank; assesses accuracy on observed *and* novel
taxa, **calibration of confidence estimates**, and compute cost. Finds
phylogenetic placement (EPA-ng) best for COI, **composition-based classifiers
(SINTAX, RDP-NBC, IDTAXA) best for fungal ITS**, and short sub-regions ≈
full-length barcodes. Consequences for us: (i) their ITS result is the **foundation this project
builds on** — see §4: because composition-based classifiers win on ITS, the
algorithm-family question is settled and this project deliberately zooms to
four representative methods; (ii) their novel-taxon benchmarking is exactly
critique T2 — cite and differentiate, or add the one-species-out test (P7);
(iii) they do not vary the reference database, nor its curation, have no
negative-control metrics and do not combine classifications — that is our
opening.

**Winand, D'hooge, Van Uffelen, …, Vanneste & De Keersmaecker (2025)** —
*Investigating fungal diversity through metabarcoding…: assessment of ITS1 and
ITS2 Illumina sequencing using multiple defined mock communities with
different classification methods and reference databases*, BMC Genomics 26:729
(doi:10.1186/s12864-025-11917-y). **The closest overlap with Q2.** 37 fungal
defined mock communities (51 species), crossed with ITS1 vs ITS2, **UNITE vs
a curated in-house database (BCCM/IHEM)**, BLAST vs mothur, species vs genus:
56–100 % of species correctly assigned depending on the combination. Also
shows with Sanger data that **ITS1 and ITS2 alone do not discriminate some
genera** — a truth limitation very like our ties (Q2c) and placeholder
epithets. Differentiators for us: 4 algorithm families and 21 parameter rows
vs 2 tools; a designed 3×3 curation ladder (all / Fungi / Fungi+rep) vs one
custom database; per-unit Sanger truth with scored controls. Their weakness is
our weakness reversed: **37 mocks vs our 2** (critique T4).

**Graetz, … (2025)** — *Benchmarking fungal species classification using
Oxford Nanopore long-read ITS metabarcodes*, Fungal Genetics and Biology
(PMID 41429697). Mock of 54 Dikarya species sequenced on ONT; **eight
classifiers** from alignment and k-mer methods to machine learning; sensitivity,
precision and diversity recovery at species level. Headline: classifiers that
**derive similarity thresholds from the provided reference** are more accurate
and better at placing "unknown" species near their true origin. This supports
our Q1 framing (thresholds are the object, not a nuisance) and warns that our
method set misses the ML family (see §2).

**Bayer, Bennett, Nester, …, Rauschert (2025)** — *A Comprehensive Evaluation
of Taxonomic Classifiers in Marine Vertebrate eDNA Studies*, Molecular Ecology
Resources 25:e14107 (doi:10.1111/1755-0998.14107). Nine classifiers on 12S /
16S / COI with **curated reference databases and exclusion-database tests**,
three positive and **two negative control datasets**. MMSeqs2 / Metabuli beat
BLAST; database exclusion tests show some tools are far more prone to false
positives when a species is absent from the reference. Different kingdom and
marker, but methodologically the nearest precedent to our `external_` controls
and `ext_fungi` / `ext_correct` metrics — cite when presenting them.

**Pauvert, Buée, Laval, …, Vacher (2019)** — *Bioinformatics matters: the
accuracy of plant and soil fungal community data is highly dependent on the
metabarcoding pipeline*, Fungal Ecology 41:23–33
(doi:10.1016/j.funeco.2019.03.005). **The source of our Pauvert mock** (189
Dikarya strains, 181 species, its Sanger ITS1 database now our per-unit truth)
and a benchmark of **360 denoising pipelines** — not of taxonomic assignment.
Natural framing sentence: *where Pauvert et al. (2019) stopped — the pipeline
that produces the ASVs — this project starts: the assignment step itself, on
the same mock and its Sanger sequences.* Companion data paper for the second
mock: **Tedersoo et al. (2026)**, MER (doi:10.1111/1755-0998.70189), 103-species
mock across Illumina / PacBio / ONT (its Table S1 is our Tedersoo truth).

## 2. The rest, one line each

*Classifier papers (the four methods in the grid):*
- **Wang, Garrity, Tiedje & Cole (2007)**, RDP naive Bayesian classifier
  (doi:10.1128/AEM.00062-07) — our `dada2::assignTaxonomy` method.
- **Edgar (2016)**, SINTAX (bioRxiv doi:10.1101/074161) — k-mer + bootstrap,
  no training; our `sintax`. Notes over-classification of novel taxa.
- **Allard et al.**, SPINGO (Bioinformatics, doi:10.1093/bioinformatics/btv614)
  — close-range species assignment; a candidate method we do not run.
- **Somervuo et al. (2016/2017)**, PROTAX (Mol Ecol Resour) — probabilistic
  placement of unknowns; the family Hleap tests and we do not.
- **Murali et al. (2018)**, IDTAXA (PeerJ) — training-set aware naive Bayes
  with an explicit "unclassified" decision; excluded here on cost-parity
  grounds (decision 2026-05-19), but see Orsholm 2026 where it is among the
  best on fungal ITS — reviewers may ask for it.

*Earlier method benchmarks (all non-fungal or 16S unless noted):*
- **Lan et al. (2012)**, 16S classifier comparison — the original k-fold /
  leave-one-out framing.
- **Bokulich et al. (2018)**, q2-feature-classifier (Microbiome,
  doi:10.1186/s40168-018-0521-5) — stratified CV, "novel taxon" runs; our CV
  design (D1a) is consciously Bokulich-style.
- **Edgar (2018)**, TAXXI (PeerJ, doi:10.7717/peerj.4652) — identity-based
  splits; argues k-fold overestimates low-rank accuracy. Our D1a accepts that
  flattery by decision (critique T7).
- **Kaehler et al. (2019)** — CNN classifier benchmark (marker-gene
  classification beyond sequence similarity).
- **McIntyre et al. (2017)** — metagenomic classifier benchmark + ensembles;
  out of scope (shotgun), cited for the ensemble lineage of our Q3.

*Fungal-ITS-specific classifiers (2024–2026 wave our method set misses):*
- **Miranda, Azevedo, Ramos, Renard & Piro (2024)**, HiTaC (BMC Bioinformatics,
  doi:10.1186/s12859-024-05839-x) — hierarchical classifier trained on the
  taxonomy tree for fungal ITS; evaluated on TAXXI against 17 classifiers.
- **HFTC (2026)**, Frontiers in Genetics (doi:10.3389/fgene.2025.1650244) —
  hierarchical random forest on fungal ITS; reports MCC 95 % and beats RDP /
  SINTAX / mothur / QIIME2 at species level (on its own closed test set).
- **ClassifyITS (2026)**, bioRxiv doi:10.64898/2026.08.31.748297 (CRAN) —
  alignment-based, **rank-specific identity/e-value cutoffs** taken from
  Tedersoo's best-practice papers (2021, 2022); conservative by design.

*Database-curation papers (the precedent for Q2 / `fungi + rep`):*
- **Gold et al. (2021)**, Mol Ecol Resour (doi:10.1111/1755-0998.13450) —
  curated, custom 12S fish databases sharply improve assignment over
  defaults.
- **Mugnai et al. (2023)**, PeerJ (PMID 36643652) — *"Be positive: customized
  reference databases… balance false taxonomic assignments"*: custom COInr
  variants (e.g. removing a clade) change false-assignment rates. Nearest
  published statement of our Q2 hypothesis (P1a of the critique), on COI.

*Unit-of-analysis paper:*
- **Tosadori & Bosch (2026)**, Web Ecology 26:47–59 — ASVs vs OTUs through
  six pipelines on in silico communities; the Q5 comparison already published,
  which is why Q5 should be framed as confirmation on a mock with per-sequence
  truth.

## 3. Comparison matrix

| Study | Marker / group | Truth | Method axis | Database axis | Parameter axis | Negative controls | Novel taxa | Cost |
|---|---|---|---|---|---|---|---|---|
| Lan 2012; Bokulich 2018; Edgar 2018; Kaehler 2019 | 16S / COI / general | reference labels | 3–17 classifiers | fixed | tuned or defaults | no | Bokulich, Edgar: yes | some |
| Hleap 2021 | COI, animals | **per sequence** (curated) | 5 families | fixed (subsampled) | tuned on ⅓ | shuffled seqs | partly | no |
| Orsholm 2026 | COI + **fungal ITS** | per sequence | ~6 families incl. phylogenetic, PROTAX | fixed | defaults + confidence calibration | no | **yes (design)** | yes |
| Winand 2025 | fungal ITS1/ITS2 | 37 DMCs, Sanger | 2 tools | **2 (UNITE vs curated)** | few | mock only | no | no |
| Graetz 2025 | fungal ITS (ONT) | 54-species mock | 8 incl. ML | fixed | threshold strategies | "unknown" species | partly | no |
| Bayer 2025 | 12S/16S/COI, vertebrates | simulated + mock | 9 classifiers | **curation + exclusion tests** | defaults | **2 negative sets** | yes (exclusion) | no |
| Mugnai 2023; Gold 2021 | COI / 12S | assignments compared | 1–4 | **custom variants** | defaults | no | no | no |
| Miranda 2024; HFTC 2026; ClassifyITS 2026 | fungal ITS | TAXXI / own sets | 1 (vs 5–17) | fixed | rank-specific cutoffs | no | ClassifyITS: partly | some |
| Pauvert 2019 | fungal ITS1 | 189-strain mock | — (360 *pipelines*) | fixed | pipelines | mock only | no | no |
| **This project** | fungal ITS1 + full ITS | **per unit** (Sanger match, ties, GBIF names) | 4 families (dada2, sintax, lca, blastn) | **3 × 3 curation ladder** (2 providers) | **21 rows, swept not tuned** | **shuffled + 100 external, scored per rank** | only via CV | **yes (wall/RAM/CO₂)** |

## 4. Where `benchmar_assign_taxo` sits

**The framing (developer, 2026-09-22).** Orsholm et al. (2026) settles the
algorithm-family question for fungal ITS: composition-based classifiers
(SINTAX, RDP-NBC, IDTAXA) beat phylogenetic placement and generic similarity
search there, because ITS is too variable to align reliably but rich in
short-motif (k-mer) signal. There is nothing to gain from re-running a
17-classifier shoot-out on ITS. **This project takes that result as its
starting point**: it zooms the method axis to four representatives and puts
the complexity where the literature has none — the composition of the
reference database, and the voting across independent classifications.

The moves:

1. **Zoom on methods (deliberate narrowing, not a gap).** `blastn`, `sintax`,
   `lca`, `dada2`: the family Orsholm hails for ITS (dada2's RDP naive Bayes
   and SINTAX's k-mer bootstrap are composition-based) plus two
   similarity-search methods (vsearch LCA, blastn + vote). One nuance to write
   carefully: **LCA is search-based like blastn** (vsearch alignment, then LCA
   of the top hits), not composition-based — the honest contrast is 2 × 2:
   composition vs similarity search × answering with the nearest name vs
   abstaining behind a threshold. blastn is then not the odd one out but the
   abstaining contrast, and the Q1 sweeps measure what Orsholm left as a black
   box: how each family trades abstention for error along its threshold
   (Graetz 2025's "reference-informed thresholds win", made quantitative).
2. **Database composition as the first complexity layer.**
   `[EUKARYOME | UNITE (± singletons)] × [all | Fungi | Fungi + rep]` — a
   designed curation ladder on both providers, scored per unit with negative
   controls (`ext_fungi`, `ext_correct`). Gold 2021 and Mugnai 2023 showed
   curation matters on COI/12S with 1–2 databases and default parameters;
   Winand 2025 showed it on fungal mocks with 2 databases and 2 tools. Nobody
   crosses the ladder with four method families and per-sequence truth.
3. **Voting across classifications as the second complexity layer (Q3).**
   Does a vote over several method × database classifications beat the best
   single one? McIntyre 2017 did ensembles on shotgun metagenomes; nobody does
   it on fungal ITS with per-unit scoring and scored controls. The
   hypothesis worth stating up front: **a vote helps through error
   independence, so it should gain from algorithm-family and provider
   diversity, not from column count** — testable at zero compute by running
   `resolve_taxo_conflict()` on column subsets (composition-only vs mixed
   families, one provider vs both). Beware the redundant-voter trap (critique
   T1): the nested UNITE pair and two search-based columns count shared
   evidence twice.
4. **One language, one open-source ecosystem (developer, 2026-09-22).** The
   four methods and the voting strategies all live in the **pqverse
   metapackage**: `MiscMetabar` wraps the four assigners
   (`add_new_taxonomy_pq(method = "dada2" | "sintax" | "lca" | "blastn")`,
   dada2's `assignTaxonomy` being an Import of MiscMetabar), and `comparpq`
   provides both the consensus strategies (`resolve_taxo_conflict()`: rel.
   majority, preference, unanimity, consensus, abs. majority) and the
   per-unit scoring (`tc_metrics_unit()`). Every cell of the grid — assignment,
   voting, scoring — is therefore reproducible in **R with one dependency
   stack**, through uniform arguments (`min_bootstrap`, `lca_cutoff`,
   `vote_algorithm`, `min_id`, `strict`, `nb_agree_threshold`), which is what
   makes the 21-row sweeps and the subset votes cheap. No benchmark in §1–2
   can say this: they mix stand-alone binaries, QIIME2 plugins and Python
   packages. **Caveat to state in the paper** (self-implementation): the
   benchmarking harness and part of the benchmarked software share an author,
   so (i) declare it, (ii) pin the exact package commits (the project already
   records them per store via the `session_info` target), and (iii) note that
   the underlying engines are independent third-party tools (vsearch, BLAST+,
   dada2), so the wrapper cannot flatter them — it can only flatter its own
   parameter plumbing, which `tests/test_assign_compute.R` checks against the
   direct calls.

**What no published benchmark does together** (the claims to make):

1. **Database curation as a designed factor** with the over-assignment to
   Fungi measured per rank (`ext_fungi`, `ext_correct`).
2. **Negative controls as scored units in the confusion matrix** (Hleap's
   shuffles, plus outgroups), one cell per sequence per rank as in Hleap's own
   code.
3. **Voting across method × database classifications** tested against the best
   single column, with the diversity-of-voters hypothesis above.
4. **Parameter sensitivity as the result** — 21 rows derived from one
   computation per method × database × input, instead of tuned optima (Hleap)
   or defaults (most of the literature).
5. **Cost and CO₂ per computation** as a first-class axis — Orsholm 2026 is
   the only other recent benchmark to report resource use.
6. **The whole grid in one language and one dependency stack** (R / pqverse:
   MiscMetabar for the four assigners, comparpq for voting and scoring), with
   the package commits pinned per run.

**Overlaps to handle explicitly in the manuscript:**

- *Orsholm et al. 2026* — now the **foundation, not a competitor**: cite as
  the reason for the four-method zoom and as what the project does *not*
  re-claim (composition-based winning on ITS is their result). What remains
  open in their paper — database composition, controls, voting, cost — is this
  paper's scope.
- *Winand et al. 2025* for Q2 — cite first, and state the delta (9 databases
  in a curation ladder × 4 methods × per-unit truth vs 2 databases × 2 tools
  × 37 mocks). Their 37 mocks are better replication than our 2; do not claim
  generality we do not have (critique T4).
- *Tosadori & Bosch 2026* for Q5 — frame as replication with per-sequence
  truth and abundance-aware reading, or move to supplementary.
- *McIntyre et al. 2017* for Q3 — the ensemble precedent is shotgun
  metagenomes; state that Q3 is the metabarcoding, per-sequence-truth version.

**Gaps the literature confirms** (back to `critique_2026-09-22.md`):

- **Novel taxa** (T2/P7): Orsholm 2026 shows this is now the field's central
  question; a closed-book mock benchmark will be asked for it. The
  one-species-out test is the cheapest answer — and it doubles as the perfect
  Q3 stress test (voting should help most when every single method is wrong in
  a different direction).
- **Method set**: with the zoom framing this is no longer a gap but a choice —
  which must be stated as such ("we take Orsholm et al.'s answer for the
  algorithm families and study what it leaves open"). HiTaC, HFTC and
  ClassifyITS (2024–26) then belong to the family question, i.e. to the
  discussion, not the grid.
- **Confidence calibration**: Orsholm scores *calibration* of support values,
  not just accuracy. Our threshold sweeps already give (threshold, error,
  abstention) curves at zero cost — one calibration figure answers their
  question for the abstaining methods.

**Suggested positioning paragraph** (for the Introduction):

> For fungal ITS, the choice of classification algorithm family is settled:
> composition-based classifiers outperform similarity search and phylogenetic
> placement (Orsholm et al. 2026), because ITS is too variable to align
> reliably but rich in k-mer signal. What remains undecided is what the
> practitioner actually controls once the software is chosen: what the
> reference database is made of, and whether combining several classifications
> helps. Zooming on four representative methods — two composition-based
> (dada2's RDP naive Bayes, SINTAX) and two thresholded similarity searches
> (vsearch LCA, blastn + vote) — we cross a designed curation ladder of the
> two main ITS references (EUKARYOME and UNITE, whole / Fungi-only / Fungi +
> non-fungal representatives) with voting strategies, on two fungal mock
> communities with per-sequence Sanger truth. Per rank and with scored negative
> controls, we measure how reference curation and voting trade accuracy against
> abstention and against falsely naming non-targets, at what compute and CO₂
> cost. All four methods, the voting strategies and the scoring are implemented
> in a single open-source R ecosystem (the pqverse metapackage), so the whole
> grid is reproducible in one language and one dependency stack.
