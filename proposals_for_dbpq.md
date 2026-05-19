# Functions to consider adding to dbpq

The benchmark project keeps these local because dbpq does not (yet) cover
them. Each is a thin operation that fits dbpq's scope ("download, format,
summarize, and modify FASTA reference databases") and would let
`make_databases.R` and `R/functions.R` collapse to pure dbpq calls.

Ordered roughly by usefulness × ease of porting.

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

**Implementation sketch** — add `if (invert) grep_args <- "-v"` to the shell
pipeline; everything else stays. The harder case is "give me the non-Fungi
records as full two-line entries" because `grep -v 'Fungi' -A 1` does not do
what you want — that flag only affects matched lines. The current bash
trick (extract Fungi block, then exclude *those exact lines* from the input)
is what `filter_db(invert = TRUE)` would need to reproduce.

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

## 3. `cluster_db(ref_fasta, output, identity)`

vsearch `--cluster_fast` wrapper that returns the centroid fasta. dbpq
already exports `find_vsearch()` and `is_vsearch_installed()` for the
machinery; only the wrapper is missing.

**Local call site** — `make_databases.R::derive_clustered()`.

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
local `derive_mini()` uses `head -n 10000` and tolerates the truncation
because the inputs are all two-line-per-record.

**Local call site** — `make_databases.R::derive_mini()`.

**Suggested signature**
```r
head_db(ref_fasta, output, n = 1000)
```

**Implementation sketch** — `Biostrings::readDNAStringSet(ref_fasta)[seq_len(n)]`
then `writeXStringSet()`. Order of magnitude slower than `head` on huge
files; for the smoke-test use case this is fine.

---

## 5. `remove_parens_db(ref_fasta, output)`

Strip parenthesized synonyms from headers, i.e. `sed 's/([^)]*)//g'`. The
EUKARYOME v1.9.3 release ships taxonomy strings like
`Genus_name (synonym_name)` that confuse downstream parsers.

**Local call site** — `make_databases.R::derive_no_parens()`.

**Suggested signature**
```r
remove_parens_db(ref_fasta, output)
```

This is narrow but it is the kind of database-specific quirk dbpq's
download/format functions already absorb for other releases — better to
keep the fix close to the format-detection logic.

---

## 6. `create_fake_pq_from_refseq(refseq, taxonomy_in_names = TRUE, taxa_ranks)`

Build a degenerate `phyloseq` object from a reference fasta: 1-column
`otu_table` of zeros, `refseq` slot with the sequences, and (optionally) a
`tax_table` reconstructed from the sintax-style header. Used by
`cross_val()` to wrap each test-fold's sequences in a phyloseq so they can
be passed to `assign_*` functions.

**Local call site** — `R/functions.R::create_fake_pq_from_refseq()`, called
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
