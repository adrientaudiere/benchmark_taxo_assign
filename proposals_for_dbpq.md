# Functions to consider adding to dbpq

The benchmark project keeps these local because dbpq does not (yet) cover
them. Each is a thin operation that fits dbpq's scope ("download, format,
summarize, and modify FASTA reference databases") and would let
`make_databases.R` and `R/create_fake_pq_from_refseq.R` collapse to pure dbpq calls.

Ordered roughly by usefulness × ease of porting. **Four live proposals**
(§1, §2, §4, §6), two withdrawn (§3, §5) and the dbpq defects of §0.

**Rule (2026-09-22, ROADMAP C22-S2): do not upstream dead code.** A proposal is
live only while its local call site is reachable from `derive_all_variants()`.
The numbering is kept stable — a withdrawn item keeps its number and says why,
because other documents cite these sections.

---

## 0. Defects found while rebuilding the databases (2026-09-11)

These are bugs in dbpq, worked around locally in `make_databases.R`.

**`format2sintax()` / `format2dada2()` drop the UNITE kingdom.** The `unite`
branch of `.parse_tax_header()` splits the header on `;` and takes the first
field as the identifier. UNITE general-release headers are
`Name|Acc|SH|type|k__Kingdom;p__…`, so the identifier becomes
`…|k__Kingdom` and the kingdom is lost: the sintax output starts at `p:` and
the positional dada2 output starts at the phylum. Local workaround:
`general_headers_fixed()` rewrites `|k__` into `;k__` before conversion.
Suggested fix: in the `unite` branch, split identifier and taxonomy at the
first `k__` (or at the last `|` before it) rather than at the first `;`, and
add a test with a real general-release header.

**EUKARYOME name qualifiers pass through `format2sintax()` /
`format2dada2()`.** EUKARYOME general headers carry `g__Lactarius(Fungi)`,
`g__(Candida)` (sometimes `g__(Candida]`), `f__Gonostomatidae(Sporadotrichida)`,
`g__Mortierella.s.str` or `f__Xxx.nom.prov` (v2.1: 11 413 ITS and 5 586 SSU
headers with parentheses). dbpq copies them into both outputs, where they
break the consumers: vsearch prints `g:Lactarius(Fungi)(0.97)`, whose
bootstrap `MiscMetabar::assign_sintax()` reads as NA (so `min_bootstrap`
never filters that genus), and dada2 or lca return names that never match a
curated taxonomy. This belongs to dbpq (the database is malformed for the
sintax and dada2 formats), not to the assignment functions. Local workaround:
`general_headers_fixed()` (sed rules, tested on real headers). Suggested
addition: a qualifier clean-up in the `eukaryome` branch of
`.parse_tax_header()` (keep the name, unwrap `(Name)`, drop `.s.str` and
`.nom.prov`), plus a warning when a converted rank value still contains `(`,
`)`, `[` or `]`.

**`filter_db()` matches the pattern anywhere in the header.** Filtering on
`"Fungi"` keeps `f:Fungiidae` corals, `k:cf.Fungi` and names such as
`fungiformis` (965 non-Fungi records in EUKARYOME ITS v2). Local workaround:
`derive_kingdom_only()` builds a kingdom-anchored pattern per header format.
Suggested addition: `filter_db(..., rank = "k", value = "Fungi")`, parsing the
header with `.parse_tax_header()` so the match is on the rank value.

**`download_file(timeout = Inf)` warns on every download.** `options(timeout =
Inf)` is coerced to an integer by `utils::download.file()`, which emits
"NAs introduits lors de la conversion automatique en 'integer'" and leaves the
effective timeout undefined. Suggested fix: translate `Inf` into a large
finite value (e.g. `.Machine$integer.max`) before setting the option.

**EUKARYOME v2.1 archives are nested.** `General_EUK_<marker>_v2.1.zip`
contains a `.7z` archive, not a FASTA, so `download_eukaryome_db()` returns a
file that still needs two extraction steps. Local workaround in
`make_databases.R::extract_single_fasta()`. Suggested addition: an `extract`
argument mirroring `download_unite_db()`.

**`download_unite_db(doi = …)` saves an HTML page.** UNITE DOIs resolve to a
PlutoF landing page. Fixed in dbpq on 2026-09-11 with the new `url` and
`extract` arguments; the `doi` argument could query
`https://api.plutof.ut.ee/v1/public/dois/?identifier=<doi>` to find the
archive URL automatically.

---

## 1. `filter_db(..., invert = TRUE)`

Extend the existing [`filter_db()`](../pqverse_pkg/dbpq/R/modify.R) so the
caller can keep records that **do not** match the pattern. The current
implementation only supports inclusion (`grep -i 'pattern' -A 1`); adding an
`invert` arg would let it run `grep -v` instead.

**Local call site** — `make_databases.R::derive_no_pattern()` (currently
shells out to `grep -vFf - … < <(grep -A1 'pattern' …)`).

**Suggested signature**
```r
filter_db(ref_fasta, pattern, output = NULL,
          invert = FALSE,
          force_two_lines_per_seq = TRUE,
          keep_temporary_files = FALSE)
```

**Implementation sketch** — `grep -v 'Fungi' -A 1` does not work for inversion
because `-A 1` only applies to *matched* lines, not the unmatched ones you want
to keep. `grep -vFf` (exclude exact lines from a prior grep pass) works only
for two-line FASTA and breaks silently for multi-line sequences.

The correct approach (same as `make_databases.R::derive_no_pattern()` after the
2026-05-26 fix): normalize to two-line format first, then use awk to exclude
matching records:

```bash
cat input.fasta | sed ':a;N;/>/!s/\n//;ta;P;D' > tmp_2line.fasta
awk '/^>/{keep=!/pattern/} keep{print}' tmp_2line.fasta > output.fasta
```

The `force_two_lines_per_seq` normalisation step already in `filter_db` covers
this — `invert = TRUE` just needs to flip the awk condition.

---

## 2. `subset_balanced_db(ref_fasta, output, n, by_rank, seed)`

Sample `n` records from a fasta, picking one record per `by_rank` value
first and filling the rest randomly. Used to build the
`fake_ref_asv_100.fasta` non-Fungi control set.

**Local call site** — `make_databases.R::derive_fake_ref()`.

**Suggested signature**
```r
subset_balanced_db(ref_fasta, output, n = 100,
                   by_rank = "Phylum",
                   input_format = "auto",
                   seed = NULL)
```

**Implementation sketch** — the current local version regexes `p__…;c__` to
extract phylum from UNITE-format headers. A dbpq version should reuse
`detect_tax_format()` + the `.parse_tax_header()` machinery already present
in `dbpq/R/format.R` so it works on any supported format. Output is a fasta
written via `Biostrings::writeXStringSet()`.

---

## 3. `cluster_db(ref_fasta, output, identity)` — **WITHDRAWN 2026-09-22**

**Its local call site is dead code.** `make_databases.R::derive_clustered()` is
defined and never called: the only other occurrences of its name in the project
are its own `stop()` message and the comment at `make_databases.R:31` saying the
orchestration does not call it. Clustered databases left the grid with the
2026-09-17 design (`config.R::benchmark_dbs` has no clustered level), so there is
no local user to port. Verified by grep on 2026-09-22.

Withdrawn rather than deleted because the idea is sound and cheap if a clustered
reference ever returns to the design: dbpq already exports `find_vsearch()` and
`is_vsearch_installed()`, so only the wrapper below would be missing. **Do not
propose it upstream in this state** — write the caller first.

**Former local call site** — `make_databases.R::derive_clustered()`, **deleted
on 2026-09-22** (ROADMAP C22-S2). The body is in the git history and the sketch
below is enough to rewrite it.

**Suggested signature**
```r
cluster_db(ref_fasta, output, identity = 0.99,
           threads = 1, verbose = TRUE)
```

**Implementation sketch** —
```r
status <- system2(
  find_vsearch(),
  args = c("--cluster_fast", ref_fasta,
           "--id", format(identity, nsmall = 2),
           "--centroids", output,
           "--threads", threads)
)
```
plus the standard "check vsearch on PATH" guard and an output-already-exists
short circuit.

---

## 4. `head_db(ref_fasta, output, n)`

Take the first `n` *records* from a fasta — not the first `n` lines, which
silently truncates the last record when sequences span multiple lines. The
local `derive_mini()` now counts records with awk
(`awk '/^>/{count++} count > n {exit} {print}'`, tested in
`tests/test_make_databases_helpers.R`), so a dbpq port can start from that.

**Local call site** — `make_databases.R::derive_mini()`.

**Suggested signature**
```r
head_db(ref_fasta, output, n = 1000)
```

**Implementation sketch** — `Biostrings::readDNAStringSet(ref_fasta)[seq_len(n)]`
then `writeXStringSet()`. Order of magnitude slower than `head` on huge
files; for the smoke-test use case this is fine.

---

## 5. `remove_parens_db(ref_fasta, output)` — **WITHDRAWN 2026-09-22**

Would have stripped parenthesized synonyms from headers, i.e.
`sed 's/([^)]*)//g'`, for EUKARYOME taxonomy strings like
`Genus_name (synonym_name)`.

**Withdrawn for two independent reasons.** (i) **Superseded on 2026-09-11** by
the EUKARYOME qualifier item of §0, which is the correct fix: a blind
`s/([^)]*)//g` empties `g__(Candida)` and misses `g__(Candida]`, `.s.str` and
`.nom.prov`, whereas `general_headers_fixed()` handles all four. (ii) **Its local
call site is dead code**: `make_databases.R::derive_no_parens()` is defined,
never called and never tested (verified by grep on 2026-09-22).

What should reach dbpq is the §0 qualifier clean-up in the `eukaryome` branch of
`.parse_tax_header()`, not this function.

**Former local call site** — `make_databases.R::derive_no_parens()`, **deleted
on 2026-09-22** (ROADMAP C22-S2).

---

## 6. `create_fake_pq_from_refseq(refseq, taxonomy_in_names = TRUE, taxa_ranks)`

Build a degenerate `phyloseq` object from a reference fasta: 1-column
`otu_table` of zeros, `refseq` slot with the sequences, and (optionally) a
`tax_table` reconstructed from the sintax-style header. Used by
`cross_val()` to wrap each test-fold's sequences in a phyloseq so they can
be passed to `assign_*` functions.

**Local call site** — `R/create_fake_pq_from_refseq.R::create_fake_pq_from_refseq()`, called
inside `R/cross_val.R`.

**Suggested signature**
```r
create_fake_pq_from_refseq(references_sequences,
                           taxonomy_in_names = TRUE,
                           taxa_ranks = c("Kingdom","Phylum","Class",
                                          "Order","Family","Genus","Species"))
```

This sits awkwardly between dbpq and `MiscMetabar` — it produces a phyloseq
object (MiscMetabar territory) but operates on a fasta (dbpq territory). My
preference is dbpq because the inputs are databases; the result happens to
be a phyloseq so the assignment helpers can consume it.

---

## Not proposed (kept local on purpose)

- `R/combine_taxo_assignments.R::combine_taxo_assignments()` — extracts
  suffixed columns from N assignment phyloseqs and cbinds them onto a base
  tax_table. This is a benchmark-specific consensus step, not a generic
  phyloseq verb; better placed in `comparpq` if it ever leaves this project.
- `R/cross_val.R::cross_val()` / `cross_val_param()` — k-fold cross-
  validation of an assignment method. The body wires together
  `MiscMetabar::assign_*()` + `dbpq::format2dada2()` + dada2's
  `assignTaxonomy()`. It is more of an evaluation harness than a database
  operation; arguably belongs in `comparpq` (which already houses
  `tc_metrics_mock` and the resolve/compare verbs).
