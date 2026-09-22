# Experimental design: what is tested, against what, at what cost

> **Updated 2026-09-22 (ROADMAP 0.10).** The counts of §1, §2 and §6, the
> dataset table of §4 and the two `_cut` implications of §5.2 and §5.3 now
> describe the design of 2026-09-17: 9 databases × 3 inputs, 96 computations
> and 504 rows per dataset, no trimmed database in the grid, per-unit truth for
> both mocks, and no miaSim branch.
>
> **One thing is still stale, and it is the important one: the measured costs
> of §3 were read on the 44-computation grid** (2026-09-16). They are the only
> figures anyone has, so they are kept as an order of magnitude — dada2
> dominates, a whole release peaks above 50 GB — but the totals now describe
> less than half the work, and the per-database split was never usable (see
> §3). Remeasure after the first production run on the new grid.

Written 2026-09-16, at the developer's request, to settle the implications of
the benchmark's design **before** launching large computations. It is a design
document, not an implementation note: it records what each objective actually
measures, what it costs, and which questions are still open.

The benchmark's central question, which every decision below serves:

> Which **method × database (including database modifications) × parameter**
> combinations give the best trade-off between assignment quality, compute
> time and compute power?

Companion documents: `docs/reference_databases.md` (how each reference is
built), `docs/diagnostics_cv_and_databases.md`
(cross-validation reference design and memory measurements).

## 1. The framing finding: objectives are not independent cost units

It is tempting to give each objective its own dataset and script. On the mock
community that would save nothing, and splitting the objectives into separate
`{targets}` projects would cost **more**, because they share their
computations:

- The 96 computations are 4 methods × 9 databases on the ASVs (36), the same
  on the OTUs (36) and 4 methods × 6 Fungi databases on the ITSx input (24).
- **Q2 requires all nine databases by definition.** Q1, Q3 and Q5 re-read
  exactly the same columns of `d_all_taxo` / `d_all_taxo_otu` and recompute
  nothing.
- The 21 parameter rows per database and input are derived from those 96
  computations (`R/assign_compute.R`, ROADMAP S8.2 and 0.4), so Q1.3's
  bootstrap sweep, Q1.4's vote sweep, the lca `lca_cutoff` sweep and the
  blastn `min_id` sweep are **free**.

The axis that really multiplies cost is the **number of datasets**: each new
one replays the whole 96-computation grid. With the developer's decisions of
2026-09-16 (full grid on every dataset, ITSx everywhere), total cost is
**linear in the number of datasets**, and that count is the only remaining
lever.

## 2. Objective → grid map

| Objective | What varies | Truth | Own computations | Marginal cost |
|---|---|---|---:|---|
| **Q1** algorithm × parameters | 21 parameter rows | mock | 0 | **none** (re-reads Q2) |
| **Q2** database × simplification | 9 databases | mock | 36 | **the ASV bill** |
| **Q3** consensus voting | existing columns | mock | 0 | **none** (post-hoc) |
| **Q4** ITSx preprocessing | preprocessed input | mock | 24 | **+67 %** of the ASV bill |
| **Q5** ASV → OTU | OTU input `d_vs` | mock | 36 | **+100 %** of the ASV bill |
| **D1a** cross-validation | 7 DB × 4 methods × 2 variants | the database itself | 56 targets | 47 min, **≥ 45 GB** for dada2 |
| **M1** performance metrics | — | — | 0 | none |
| **M2** compute metrics | — | — | 0 | none (`with_autometric()` already wraps every target) |

**Only Q2, Q4, Q5 and D1a cost anything.** Q1, Q3, M1 and M2 are re-readings of
material already on disk. Q5 changed side on 2026-09-17: the OTUs are now
assigned in their own right (objectives_design decision 4) instead of being
derived from the ASV assignments, which is why it carries 36 computations.

Verify this claim with `build_values_map(dbs = db_list)`: 504 rows, 96 distinct
`compute_name` (36 ASV + 36 OTU + 24 ITSx).

## 3. Measured cost (production store, not extrapolated)

Read from `store_assign_taxo::benchmark_costs` on 2026-09-16, over the 44
distinct computations:

| Method | n | Target time | Peak RAM |
|---|---:|---:|---:|
| dada2 | 11 | **1.14 h (53 %)** | **53.1 GB** |
| sintax | 11 | 0.60 h | 1.3 GB |
| blastn | 11 | 0.22 h | 1.3 GB |
| lca | 11 | 0.19 h | 1.2 GB |
| **Total** | 44 | **2.15 h** (1 h 34 wall) | |

dada2 is the whole problem: half the time and forty times the memory of any
other method.

**Measurement defect — do not quote the per-database split yet.**
`EUK_ITS_v2.1_Fungi` comes out at 0.06 h against **1.19 h** for
`EUK_ITS_v2.1_Fungi_cut`, which holds the same records *shortened*. That is
physically impossible: the autometric logs mix the runs of 2026-09-14 and
2026-09-15, since the `_cut` computations were rebuilt on their own and the
aggregation keeps the newest file per phase (`R/autometric_helpers.R`). The
**total is reliable, the per-database split is not**. This is what ROADMAP S1.1
already announces (regenerate the M2 figures after a clean production run).
Remeasure on the next complete run **before** producing any Q2.5 or M2 figure.

The same caveat retired the estimates of the former `docs/compute_budget.md`,
which were extrapolated on 2026-09-11 from a 2025-02 sequential run at ±50 %;
that page was deleted on 2026-09-22 and what it still said about the structure
of the cost is in section 8.

## 4. Datasets

Status verified 2026-09-16. **Region is the structuring constraint**, not a
detail: it decides which databases a dataset can legitimately be assigned
against (ROADMAP D1c, second caveat).

| Dataset | Region | Truth | Status | Missing |
|---|---|---|---|---|
| Mock Pauvert | ITS1 | 189 strain rows / 181 species, per unit since 2026-09-22 (`R/mock_truth.R`) | complete | — |
| Mock Tedersoo Illumina | full ITS | 103 species, Table S1 Sanger reads; per unit too | `d_asv` built (238 ASVs) | the assignment run |
| Mock SRR30413326 (D1d) | to confirm | 16 taxa, on disk | fastq absent | `fasterq-dump` |
| In silico D1b-A | ITS1 | input fasta lineages | fastq generated (3 seeds) | `_insilico` projects |
| Tedersoo long reads | ITS1 + ITS2 + full | mock known | NextITS not run | NextITS pipeline |
| GloSED (~10 sites) | full ITS, PacBio | none | phyloseq published | loader + **count the OTUs** |
| GSSP air | **ITS2** | none | CSV published | loader + `_ITS2_cut` |
| Korhonen wood ×2 | **ITS2** | none | `RepSeq` published | loader + `_ITS2_cut` |
| Karlsson leaf 2017 | **unknown** | none | fasta published | **read the article's Methods** |

## 5. The five implications to settle

None of these is an implementation detail: each changes what the benchmark
measures.

### 5.1 "ITSx everywhere" requires a per-dataset `itsx_region`

`itsx_region` is a single global constant (`config.R:178`, read once by
`pipelines/assign_taxo.R:181`). An ITS2 dataset extracted as ITS1 would have
essentially no overlap with its queries, and a full-ITS dataset extracted as
ITS1 would throw ITS2 away.

Measured on the Tedersoo ASVs (ITSx 1.1.3, all organism profiles, 210 of 238
ASVs with both ITS1 and ITS2 detected): **193 bp of flank at the median** —
154 bp before ITS1, 39 bp after ITS2 — for a full-ITS span of 555 bp inside
749 bp ASVs. ITS1 alone is 191 bp at the median, ITS2 202 bp.

Implementation note for later: `--save_regions` accepts only
`{SSU,ITS1,5.8S,ITS2,LSU,all,none}` — there is **no `full` value** — but the 
all seems to be the one we want in that case.


### 5.2 The ITS2 datasets require `_ITS2_cut` variants — *moot for the grid since 2026-09-17*

Written when trimmed references were part of the grid. They no longer are
(objectives_design decisions 16, 24, 25): the nine databases are untrimmed, and
trimming is tested at the end, as a zoom on the best combinations, on the
Pauvert mock only (ROADMAP 0.11b). An ITS2 dataset therefore runs the same nine
databases as every other. **What the section still says correctly** is what an
ITS2 variant would cost if the zoom is ever extended to one:
`config.R::primer_min_overlap` is derived from `rev_primer_sequences`, and the
`_cut` derivation also reads `cut_discard_untrimmed` and `cut_min_length`, so it
needs its own consistent set of primer-derived constants, not a partial
override.

### 5.3 `_cut` is dropped for Tedersoo — *fixed 2026-09-17*

Developer decision of 2026-09-16, implemented as ROADMAP 0.1b: `config.R` no
longer substitutes `EUK_ITS_v2.1_Fungi_cut_full_ITS` into `db_list`,
`itsx_db_list` and `db_meta$db` when `assign_taxo_tedersoo_illumina` runs; the
`bio_datasets$cut_suffix` column stays, holding NA. `tar_manifest()` on that
project now gives the same 96 computations and 504 rows as `assign_taxo`.

Why it mattered: the grid carried two databases that are 96.9 % identical, and
Q2 would have displayed an empty comparison that reads like a result.

### 5.4 The query-side cost is unprofiled

The D1a measurements show that dada2's memory is dominated by the **reference**
(23.5 GB on `Unite_s_all_20250219`, ≥ 45 GB on `EUK_ITS_v2.1`). They say
**nothing** about the number of queries. The mock and Tedersoo both carry ~400
sequences after control injection; GloSED at 10 sites carries an **unknown**
number of OTUs.

Measure one subsample before sizing that run. Do not extrapolate from the mock.

### 5.5 The `_cut` no-op is a Q2 result, not a bug

Finding of 2026-09-16, to write into the manuscript:

> A trimmed reference is a **no-op** when the amplicon's primer sites lie in
> flanks the reference does not carry.

Evidence, measured on `EUK_ITS_v2.1_Fungi` (1 076 870 records):

| Primer pair | 5' site present | 3' site present | Median length | Records changed |
|---|---:|---:|---|---:|
| ITS1F / ITS2 (Pauvert) | 14 914 (1.4 %) | **736 014 (68 %)** | 523 → **211** | massive |
| ITS9mun / ITS4ngsUni (Tedersoo) | 23 174 (2.2 %) | **44 (0.004 %)** | 523 → **521** | **32 858 of 1 076 389 (3.1 %)** |

`_Fungi_cut` works because its 3' anchor (in the 5.8S) is present in two thirds
of the records. Tedersoo's pair has neither anchor, so 96.9 % of the records
came through unchanged and the derived file is a 1.4 GB duplicate of its
source. The corollary matters for the manuscript: **for a full-ITS amplicon, an
ITS reference is already at the right region — there is nothing to trim.**

## 6. Decisions recorded (2026-09-16)

| Decision | Consequence |
|---|---|
| Full grid on every new dataset (44 computations then, 96 since 2026-09-17) | Maximum comparability between datasets; cost linear in dataset count |
| ITSx on every dataset | +24 computations per dataset; **required 5.1**, done (ROADMAP 0.5) |
| No `_cut` arm for the Tedersoo mocks | **Required 5.3**, done (ROADMAP 0.1b) |
| Three variants wanted (full ITS, ITS1, ITS2), ITS1/ITS2 deferred | ITS1/ITS2 sub-region references cannot come from a primer pair; they need ITSx extraction of the reference |
| Biological communities serve to compare NA counts against real data | No ecological reanalysis; abundances not needed |

## 7. What must not happen

- Do **not** run `assign_taxo_tedersoo_illumina` before 5.3 is fixed: it would
  fill a store with a duplicated database arm.
- Do **not** publish Q2.5 or M2 figures from the current `benchmark_costs`:
  the per-database split is corrupted (section 3).
- Do **not** re-derive the `_cut` databases or rerun the CV merely because
  cutadapt was upgraded to 4.9 — that premise was tested and not reproduced
  (ROADMAP D1c).

## 8. Where the cost goes, and the levers (from the former `compute_budget.md`)

Kept from `docs/compute_budget.md`, deleted on 2026-09-22 because its estimates
(2026-09-11, before any production run on the current databases) were replaced
by the measurements of section 3.

**Shape of the cost.** One computation per method × database × input; every
parameter row is derived from it in the main process, so a threshold, an
`lca_cutoff` or a `min_id` costs nothing (`R/assign_compute.R`). The cost is
therefore the number of computations (96 per dataset with the 2026-09-17 grid)
times the cost of one, which grows with the size of the reference and is
dominated by dada2: 53 % of the target time and 53.1 GB of peak memory against
1.2–1.3 GB for the other three methods.

**Levers already applied** (ROADMAP S8.1–S8.5): one computation per method
(L1–L4), threads for every method (L5, `config.R::assign_threads_dada2`,
`assign_threads_fast`, `cv_threads`), the ITSx input restricted to the
Fungi-filtered databases (L7), cross-validation folds and subsample set to
their publication values (L6).

**Levers still available**, if a run has to be made smaller: fewer
cross-validation folds or a smaller `cv_max_seq`; a narrower cross-validation
grid (`config.R::cv_db_list`, which objectives_design decision 27 leaves open);
and dropping the ITSx input on the databases where Q4 turns out to be flat.
