# Diagnostics: cross-validation design and database effects

Measurements made on 2026-09-15 and 2026-09-16 on the production stores, kept
for the manuscript (Methods and Discussion). Each section states what was
measured, on which data, and what it implies. The scripts are in
`analysis/sandbox/` and are runnable from the project root with `Rscript`.

Related: `ROADMAP.md` (items S1.3, S6.9, S8.6, D1a, Q5),
`critique_fable/07_results_validity.md` (defects B21–B24), `CONTEXT.md`
(glossary entries *CV reference pool*, *trimmed CV queries*).

## 1. The cross-validation reference is the drawn pool, not the database

`cross_val()` draws `cv_oversample × cv_max_seq` records, keeps the pool up to
its `cv_max_seq`-th usable query, and then **replaces the reference by that
pool** (`dna <- dna[match(kept$pool, names(dna))]`). This is the documented
behaviour of `max_seq` ("size of the random subsample of the database"), but it
means a CV target never searches the database it names.

Production settings (`cv_max_seq = 5000`, `cv_oversample = 2`,
`cv_fold_number = cv_fold_tested = 10`), script
`analysis/sandbox/check_cv_training_size.R`:

| database | records | drawn | kept as reference | share of the database |
|---|---|---|---|---|
| `Unite_all_20250219` | 156 820 | 10 000 | 5 413 | 3.45 % |
| `EUK_ITS_v2.1` | 1 877 003 | 10 000 | 6 284 | 0.33 % |
| `EUK_ITS_v2.1_Fungi` | 1 076 870 | 10 000 | 5 295 | 0.49 % |

The pool grows until 5 000 of its records carry the ITS2 site, so its size
follows that site rate (80.9 % for `EUK_ITS_v2.1`, 93.8–95.4 % for the Fungi
files), never the size of the database. In the standard variant the tested
fold (≈ 500 records) is removed from that pool.

**Implication for the manuscript.** The D1a database axis compares the
taxonomic composition of seven databases sampled at a common reference size;
it does not compare the databases as the mock-community results do. Any
statement of the form "database X performs better in cross-validation" is a
statement about composition at ≈ 5 400 records.

## 2. Congeners are missing from the pool by construction

Same draw as above, genus taken from the sintax header
(`analysis/sandbox/check_blastn_vote.R` uses the same helper):

| database | queries whose genus is present elsewhere in the kept pool | … in the full database | distinct genera in the pool |
|---|---|---|---|
| `Unite_all_20250219` | 69.7 % | 96.2 % | 2 388 |
| `EUK_ITS_v2.1_Fungi` | 87.6 % | 100 % | 1 098 |

So about 30 % of the UNITE queries have **no** congener in their training part:
no method can classify them correctly at genus, and the ceiling of the
standard variant is set by the draw, not by the methods.

### Singleton genera: the floor that survives a full reference

Once the reference is the whole database and only the tested queries are
removed (design decided on 2026-09-16), a query is still unclassifiable at
genus when **its genus has a single record in the database**. That floor is a
property of each database and must be stated when comparing them:

| database | records | usable genus in the header | distinct genera | genera with a single record | records whose genus is a singleton |
|---|---|---|---|---|---|
| `Unite_s_all_20250219` | 266 589 | 67.6 % | 18 394 | 41.5 % | **4.2 %** |
| `Unite_s_all_20250219_Fungi` | 168 030 | 61.9 % | 5 954 | 31.5 % | 1.8 % |
| `Unite_all_20250219` | 156 820 | 73.0 % | 13 538 | 41.6 % | **4.9 %** |
| `Unite_all_20250219_Fungi` | 102 137 | 69.0 % | 5 356 | 34.5 % | 2.6 % |
| `EUK_ITS_v2.1` | 1 877 003 | 52.6 % | 26 636 | 31.4 % | 0.8 % |
| `EUK_ITS_v2.1_Fungi` | 1 076 870 | 59.4 % | 6 800 | 14.0 % | 0.1 % |
| `EUK_ITS_v2.1_Fungi_cut` | 1 075 023 | 59.3 % | 6 790 | 13.9 % | 0.1 % |

"Usable genus" excludes `unclassified`, `*_Incertae_sedis` and empty values.

Two things to write in the Methods, because they push in opposite directions:

- **UNITE loses more queries to singleton genera than EUKARYOME** (4.2–4.9 %
  of its records against 0.1–0.8 %), so its genus NA rate is structurally
  higher, whatever the method. A UNITE-versus-EUKARYOME gap at genus of a few
  points is explained by the databases, not by the classifiers.
- **EUKARYOME carries fewer usable genus labels in the first place** (52.6 % of
  the full release against 73.0 % for `Unite_all_20250219`), so its truth table
  is thinner: queries whose own header has no genus are excluded from the
  comparison rather than counted as failures.

No stratified draw is used (decision 2026-09-16): the queries stay a uniform
random sample, and these two rates are reported instead.

## 3. Why blastn abstains and lca is wrong instead

blastn leaves 87–95 % of the standard CV queries without a genus, while lca
(same vsearch alignments) leaves 1–5 % but gives a wrong genus to 37–60 % of
them. 200 trimmed queries of `Unite_all_20250219`, query records removed,
filters of the pipeline (`min_cover = 80`, bit score ≥ 50, e-value ≤ 1e-30);
scripts `analysis/sandbox/check_blastn_cv.R`,
`analysis/sandbox/check_blastn_vote.R`,
`analysis/sandbox/check_blastn_refsize.R`:

| reference | identity ≥ 90 | identity ≥ 95 (the `assign_blastn()` default) |
|---|---|---|
| full database (156 620 records) | 82.5 % assigned | 58.5 % assigned |
| random 5 413-record subsample | 26.5 % | 9.0 % |
| the CV pool, pipeline path, one fold | — | 2 % |

Best-hit statistics against the full database: median identity 96.3 %, median
query cover 100 %, 85 % of the queries keep a hit passing every filter at
identity 90. Reference size is therefore the dominant cause; the residual gap
between 9 % and 2 % is still unexplained and is logged as open in B24.

The methods differ in what they do when no close relative is available: blastn
abstains (a filter), lca and sintax return the nearest available taxon (an
answer, often wrong). That difference is a result worth stating in the
Discussion, but the magnitude measured here is an artefact of the reference
size.

Three candidate causes were tested and cleared, each by a direct measurement
(`analysis/sandbox/check_blast_pq_table.R`,
`analysis/sandbox/check_blastn_names.R`):

- `assign_blastn()`: against the full database it gives a genus to 29 of 50
  queries (58 %) and its `taxa_names` match the phyloseq taxa 30 of 30.
- `blast_pq()`: asks BLAST for `qcovs`, so "Query cover" is the query cover
  (median 100 % for best hits); the returned table has one row per hit and the
  vote branch reduces it to one row per query.
- `resolve_vector_ranks(method = "rel_majority")`: assigns on a plurality
  (`A,A,B,C,D` → `A`) and returns NA only on a tie.

Secondary: the CV passes `min_cover` but leaves `min_id` at the
`assign_blastn()` default of 95, which costs about 24 points even against the
full database (82.5 % → 58.5 %).

## 4. dada2 on the full EUKARYOME release: bootstrap collapse, not format

On the mock community, `dada2__EUK_ITS_v2.1___0.5` gives a genus to 0 % of the
195 ASVs (Kingdom 10.8 %), against 95.4 % on `EUK_ITS_v2.1_Fungi`. Neither the
file format nor the organelle records explain it.

- **Format identical.** The dada2 files of `EUK_ITS_v2.1`, `_Fungi`,
  `_Fungi_cut` and `Unite_s_all_20250219` all have 7 fields and a trailing
  semicolon on 100 % of their records.
- **dada2 does assign.** In the raw stored computation
  `compute_dada2__EUK_ITS_v2.1` (before any threshold or cleaning): 0 % NA,
  100 % of the ASVs with a Kingdom value, 84.7 % with a real (non
  `unclassified`) Genus.
- **The bootstrap collapses**, and the benchmark thresholds (0.4, 0.5, 0.6)
  then wipe everything:

| rank | median bootstrap, `EUK_ITS_v2.1` | median bootstrap, `EUK_ITS_v2.1_Fungi` |
|---|---|---|
| Kingdom | 31 | 100 |
| Genus | 17 | 68 |

- **The driver is taxonomic breadth**: 336 distinct kingdom values and 26 637
  genera in the full release, against 1 and 6 801 in `_Fungi`; 47.4 % of its
  records are `unclassified` at Genus (40.6 % in `_Fungi`). With many
  competing non-fungal taxa, no candidate wins enough bootstrap replicates.
- **Consistent with a controlled test** (`analysis/sandbox/check_dada2_long_refs.R`,
  195 mock ASVs, `minBoot = 50`, reference = 10 000 Fungi records of
  `mini_EUK_ITS_v2.1_Fungi`): 96.4 % genus alone; 50.8 % after adding the 7 871
  records > 5 kb; 9.7 % after adding 7 871 random non-Fungi records ≤ 5 kb.
  Breadth hurts more than length. The lengths of the added short records were
  not controlled, so the two factors are not fully separated.

Organelle content of `EUK_ITS_v2.1`, for the record: 2 010 `_mitochondrion`
records (0.11 %, median 30 050 bp) and 118 `_plastid` (median 3 682 bp), i.e.
3.2 % of the nucleotides and 1 248 of the 1 309 records > 10 kb. Decision
2026-09-16: keep the release as it is.

**Implication for the manuscript.** dada2 is the only method whose bootstrap
makes it abstain on a broad reference; sintax (97.4 %), blastn (99.5 %) and
lca (66.2 %) still assign on the same file. This is a Q2 result, not a
pipeline defect.

## 5. sintax cross-validation compared shuffled rows (fixed)

`cross_val()` compares the assignment table to the truth table by position,
and `assign_sintax(behavior = "return_matrix")` returns the rows in the order
of the vsearch `--tabbedout` file, which with `--threads` > 1 is the
completion order (`analysis/sandbox/check_sintax_order.R`):

| `mini_Unite_all_20250219`, 2000 leaked queries | genus good, by position | genus good, joined on `taxa_names` |
|---|---|---|
| `nproc = 1` | 1.000 | 1.000 |
| `nproc = 4` | 0.237 | 1.000 |

Fixed on 2026-09-15 by `cv_align_rows()` (`R/cv_queries.R`, tested in
`tests/test_cv_queries.R`), which orders the sintax values, the sintax
bootstraps and the lca rows on `taxa_names(fake_pq)`; `cross_val()` also stops
when a method returns another number of rows than queries. Confirmed on the
production run: sintax leaked genus 0.999–1.000 on the seven databases,
against 0.49–0.54 before. Every sintax CV result produced before that date is
invalid (B23).

## 6. ASV versus OTU (Q5): what the comparison does and does not measure

On the production cache, at Genus, mean F1 is 0.782 over the 132 single-method
ASV columns against 0.743 for the same columns on OTUs, and no column does
better on OTUs (0 of 132). This is **not** a clean test of post-clustering:
`tc_metrics_mock_vec()` counts taxa, so moving from 195 ASVs to 143 OTUs
removes the true positives carried by the non-representative ASVs *and* the
false positives they carried, while the truth table keeps its rows. A
unit-invariant metric, or a comparison restricted to the 143 representative
sequences, is needed before answering Q5.

## 7. Cost of the alternatives for the CV reference (estimate)

Per-computation times measured in `assign_taxo` against the **full** databases
(392 ASVs, `assign_threads_dada2 = 3`, one thread for the others). Reference
reading and indexing dominate, so a CV fold of ≈ 500 queries costs about the
same as one of these computations:

| database | dada2 | lca | sintax | blastn | sum |
|---|---|---|---|---|---|
| `EUK_ITS_v2.1` | 33.4 | 31.3 | 9.0 | 6.2 | 79.9 min |
| `EUK_ITS_v2.1_Fungi` | 10.4 | 6.0 | 3.2 | 4.7 | 24.3 min |
| `EUK_ITS_v2.1_Fungi_cut` | 19.6 | 4.8 | 10.4 | 5.1 | 39.9 min |
| `Unite_all_20250219` | 7.3 | 2.4 | 1.1 | 0.9 | 11.7 min |
| `Unite_all_20250219_Fungi` | 5.4 | 1.0 | 0.5 | 0.8 | 7.7 min |
| `Unite_s_all_20250219` | 9.3 | 3.5 | 1.9 | 1.2 | 15.9 min |
| `Unite_s_all_20250219_Fungi` | 5.1 | 1.3 | 1.0 | 1.0 | 8.4 min |

One pass over the 28 method × database cells is 188 min. **The 63 h figure this
paragraph used to quote is withdrawn**: it assumed ten folds and was derived
from mock assignments of 392 ASVs, before the design changed to three folds and
a full reference.

What is measured, per fold of 500 queries against the full `EUK_ITS_v2.1`:
blastn 1 min 48, sintax 2 min 26, lca 7 min 05. What is **not** measured is the
one that matters: **dada2 did not finish a single fold in 8 minutes** before
being killed. With 14 dada2 targets × 3 folds, dada2 could dominate the run by
an order of magnitude, and no honest total can be quoted until one dada2 fold
completes on a full EUKARYOME reference. The standard variant cannot share an
index between folds, since each fold removes its own records.

A genus-preserving pool (draw the queries, then keep every record sharing
their genus) would stay near the current cost only if the added congeners are
few; the pool size that rule produces has not been measured yet.

### Memory, not wall time, is the binding constraint for dada2

Measured on 2026-09-16 with `/usr/bin/time -v`, one `cross_val()` call,
`reduce_reference = FALSE`, 500 queries, one fold:

| method × database | reference records | peak memory | wall time | note |
|---|---|---|---|---|
| dada2 × `EUK_ITS_v2.1` | 1 877 003 | **≥ 45.0 GB** | killed at 7 min 58 | fold unfinished, machine at 60/62 GB and 22.9 GB of swap |
| dada2 × `Unite_s_all_20250219` | 266 589 | **23.5 GB** | 4 min 32 | completed; genus standard 0.335 good, 0.138 wrong, 52.7 % NA |
| sintax × `EUK_ITS_v2.1` | 1 877 003 | **5.9 GB** | 2 min 26 | completed; genus standard 0.790 good, 0.054 wrong, 15.6 % NA |

The contrast between the two methods is structural: dada2 loads the reference
and its k-mer profiles into R, while vsearch streams the database, so sintax
costs eight times less memory on the *same* file. It also shows what the
reference size was hiding: sintax reaches **0.79** good genus against the full
EUKARYOME, where the drawn pool gave 0.26–0.52 (§ 1).

| blastn × `EUK_ITS_v2.1` | 1 877 003 | **3.5 GB** | 1 min 48 | completed; genus standard 0.737 good, 0.072 wrong, 19.2 % NA |
| lca × `EUK_ITS_v2.1` | 1 877 003 | **6.2 GB** | 7 min 05 | completed; genus standard 0.754 good, 0.150 wrong, 9.6 % NA |

Taken together, these four runs are a preview of the production cross-validation
on the hardest database, over the same 500 queries: sintax 0.790 good genus,
lca 0.754, blastn 0.737. The behaviour that the reduced reference exaggerated
is still there, in a readable form: **lca abstains least (9.6 % NA) and errs
most (15.0 % wrong), blastn does the opposite** (19.2 % NA, 7.2 % wrong), and
sintax sits between the two on both counts.

**Memory budget for the production run.** `pipelines/cross_val.R` now runs
`dada2_ctrl` (1 worker) next to `fast_ctrl` (`config.R::cv_n_workers_fast`).
The worst case is one dada2 target on `EUK_ITS_v2.1` beside the fast ones on
the same database:

| configuration | worst case | verdict |
|---|---|---|
| 3 fast workers | 45 + 6.2 (lca) + 5.9 (sintax) + 3.5 (blastn) = **60.6 GB** | fits on paper, 1.4 GB left for the system |
| 2 fast workers | 45 + 6.2 + 5.9 = **57.1 GB** | safe — **adopted 2026-09-16** |

Every figure is now measured rather than projected. Two caveats on the
arithmetic itself. First, `dada2_ctrl` holds **one** worker for 14 dada2
targets, so a dada2 target is resident almost continuously during the run, not
occasionally: the worst case above is the normal case. Second, the budget
assumes a fresh worker per target, which `seconds_idle = 60` does not
guarantee — a crew worker that has just finished an EUKARYOME target may still
hold its memory when the next one starts, since R returns freed memory to the
system lazily. Watch `free -g` during the first dada2 targets rather than
trusting the sum.

The dada2 figure is a **floor**, not a peak: that run was killed before it
finished its fold, so the real maximum is higher than 45 GB. That alone argued
for two fast workers, and the developer settled on them on 2026-09-16:
`config.R::cv_n_workers_fast = 2`, a constant dedicated to this project and
read by `pipelines/cross_val.R` only, rather than lowering `n_workers`, which
`pipelines/assign_taxo.R` shares and which stays at 3.

The machine has 62 GB. During that run it reached 60 GB used and 22.9 GB of
swap, so a single such target saturates it. `pipelines/cross_val.R` dispatches
`n_workers = 3` targets on one crew controller, which would put three of them
in memory at once. `pipelines/assign_taxo.R` already solves this for the mock
by giving dada2 its own one-worker controller (`dada2_ctrl`) next to a
`fast_ctrl` for sintax / lca / blastn; the CV pipeline has no such split.

**Memory is not linear in the number of records.** A first extrapolation from
the EUKARYOME point (≈ 24 kB per record) predicted ≈ 6 GB for
`Unite_s_all_20250219` (266 589 records); the measurement reached **23.5 GB
after 4 minutes** on that database, four times the prediction, so the
extrapolation was wrong and has been removed. dada2 holds the reference
sequences *and* their k-mer profiles, and the cost per record grows with
sequence length and with the number of distinct genera, not with the record
count alone. Consequence: every database must be measured, and no CV
configuration should be chosen from a projected figure.

## 7b. Implementation of the 2026-09-16 decision

Decided by the developer on 2026-09-16: Bokulich design, **removing only the
tested queries** (not a tenth of the database per fold), **no stratified
draw**, **3 folds**. The cross-validation is read as a comparison of
algorithms × databases × metrics against each other, not as an absolute
confidence value; its closeness to a flattering leave-one-out is accepted and
stated in the Methods.

Code: `cv_reference_records()` (`R/cv_queries.R`, tested in
`tests/test_cv_queries.R`) and the `reduce_reference` argument of
`cross_val()`, which defaults to TRUE so nothing changes silently for older
callers. `config.R` sets `cv_reduce_reference <- FALSE` and
`cv_fold_number <- 3L`. The queries are drawn exactly as before — same pool,
same seed, same trimming — so the only difference with the run of 2026-09-15
is the reference.

Validated on `cross_val_mini` (2026-09-16, 56 targets, 0 error, mean 6.7 s per
target). Genus, mean over the seven mini databases:

| method | good, standard | NA, standard | good, leaked |
|---|---|---|---|
| lca | 0.754 | 1.2 % | 1.000 |
| sintax | 0.726 | 13.9 % | 1.000 |
| dada2 | 0.655 | 22.9 % | 0.965 |
| blastn | 0.310 | 66.3 % | 0.971 |

The run of 2026-09-15, same targets with the reduced reference, gave 0.17–0.85
(lca), 0.15–0.70 (sintax), 0.08–0.20 (dada2) and 0–0.03 (blastn) per database.
blastn still leaves 66 % NA here because a `mini_*` file is itself only 10 000
records; on the production databases the measurements of § 3 predict roughly
40 % NA at identity 95.

## 8. How the cited benchmarks built their reference (literature, 2026-09-16)

Checked online on 2026-09-16, because none of these papers is in
`docs/references/`. **No cited study reduces the reference to a random
subsample of the database, which is what `cross_val()` does today** (§ 1).
Four designs exist, and they answer different questions.

| Study | Design | Reference used to classify a query |
|---|---|---|
| Bokulich et al. 2018 (q2-feature-classifier, tax-credit) | stratified 10-fold cross-validation of the database | the **whole** database minus the test fold |
| Bokulich et al. 2018, "novel taxon" evaluation | same folds, plus exclusion by taxonomy | the whole database minus **every** sequence sharing the query's taxon at rank *L* |
| Murali et al. 2018 (IDTAXA) | leave-one-out | the whole database minus **the one** tested sequence |
| Edgar 2018 (TAXXI) | cross-validation by identity | train/test split such that **no pair across sets exceeds a given identity**, results reported per identity bin |
| Hleap et al. 2021 | mock communities, no k-fold on the database | full assembled reference; a **taxon-aware** subsample (one representative per species, within the lowest rank shared by the mock) only for the methods needing a global alignment |

Points worth keeping for the Methods and Discussion:

- **Bokulich stratifies by taxonomy**: taxonomic labels with fewer than 10
  instances are amalgamated so every stratum is large enough. The consequence
  is the opposite of ours (§ 2): their folds are built so that a query's taxon
  stays represented, ours drops the congeners of ~30 % of UNITE queries by
  chance.
- **Edgar's motivation for TAXXI** is explicitly that leave-one-out flatters
  classifiers, because near-identical relatives stay in the reference. Our
  leaked variant is the extreme of that critique (the query's own record is
  present), and our standard variant is *harsher* than leave-one-out, since the
  draw also removes unrelated congeners.
- **Hleap et al. subsample the reference on purpose**, but by taxonomy (one
  representative per species inside the relevant clade), and they justify it as
  "a realistic use scenario" for large references with alignment-based methods.
  That precedent supports a taxon-aware pool, not a random one.
- **Hleap et al. also add shuffled sequences as true negatives**, the same
  device as our `fake_` / `external_` controls, and they report F1 and MCC,
  which is where decision 5 of the ROADMAP comes from.
- Their "pseudo 3-fold cross-validation" is a grid search for hyperparameters
  on the mock communities, not a reference-construction scheme; it should not
  be cited as a CV design.

Consequence for D1a: the current design is closest to nothing published. The
three defensible options are Bokulich (whole database minus the fold, the
costly one, § 7), Hleap (taxon-aware subsample, cheap, needs the pool rule of
§ 2), or Edgar (identity-binned split, which answers "how well does a method
do at a known distance from the reference" and would replace the standard /
leaked pair by a distance axis).
