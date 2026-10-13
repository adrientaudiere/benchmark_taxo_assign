# Results — draft

Outline of the Results, following the framing decisions of
`docs/objectives_design.md` §4 (OD 37, 43–51). Numbers are not
copied here while the chapters can still change: each subsection points to
the chapter section that computes them. Prose is added subsection by
subsection as `docs/manuscript_framing.md` is worked through.

Three uses, each with its endpoint (OD 51): filtering groups (PPV and
`NA_real` from Kingdom to Order, with the control rates), filling the
taxonomy for traits (F1 and PPV at Genus on the mock ASVs, OD 43), target
species (each mock species as the target: detected, false detections).
Genus F1 and PPV are scored on the real units only, the control rates
beside them (OD 49); every other figure is secondary. Every score is read with the units without truth set
aside, the stricter Q2b reading given beside it (OD 44). Each paired
difference carries a 95 % bootstrap interval (OD 50). Every reading,
the practical recommendation (one per use) included, is made on the three
mocks at once (OD 48, 51).

## 1. The four predictions (OD 45)

- P1a curation bias — chapter 03 (to read, ROADMAP §4).
- P1b breadth, not size — chapter 03 (to read).
- P1c abstention trade-off — not supported: SA1 (ROADMAP §2), chapter 03.
- P1d trimming no-op — holds on Genus F1, not on the controls: chapter 03
  Q2.6.

## 2. Methods and parameters (Q1)

Chapter 02 summary; held under both truth rules (chapter 02, "Under Q2b").

## 3. UNITE against EUKARYOME, and the filtering ladder (Q2, OD 46)

Paired within-source contrasts, n = 2 providers: chapter 03 Q2.1, Q2.2;
both truth rules: chapter 03, "Under Q2b".

## 4. Trimming the references (final zoom, OD 47)

Chapter 03 Q2.6.

## 5. Consensus (Q3, OD 46)

Gain against the best single column (chapter 04) and against the effective
number of voters (to compute, ROADMAP §4).

## 6. ITS region (Q4)

Chapter 07; the order holds under both truth rules.

## 7. dada2 above 65 536 taxonomies (OD 47)

`HISTORY.md` §2 "Databases building" (2026-09-30), `patches/`.

## 8. Secondary results

- Q5, ASVs or OTUs (OD 47): the OTU gain depends on how the units without
  truth are scored; chapter 08 and its "Under Q2b" section.
- Species by exact match: chapter 09.
- Costs (M2.4): chapter 05.

## 9. Recommendation per use (OD 48, 51)

One recommendation per use, read on Pauvert, Tedersoo and Rolling at once;
where the mocks disagree, say so:

- filtering groups: upper-rank PPV and NA, control rates (to compute,
  ROADMAP §4);
- traits: Genus F1 and PPV (chapters 02, 03);
- target species: per-species detection and false detections (to compute,
  ROADMAP §4), beside the pooled reading of chapter 09.
