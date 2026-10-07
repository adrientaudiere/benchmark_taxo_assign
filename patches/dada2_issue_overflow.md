**Title:** assignTaxonomy(): 32-bit overflow in get_best_genus() when the reference holds more than 65,536 distinct taxonomies

**dada2 version:** 1.40.0 (Bioconductor 3.23, git 52cd78e); the same line is on `master`.

### Summary

When a reference has more than 65,536 distinct taxonomy strings, `assignTaxonomy()` scores most of them with the wrong k-mer profile. Assignments degrade silently: bootstraps collapse, and the reported lineage is drawn at random among unrelated taxa. No error or warning is raised. Several widely used references are above this limit. The full EUKARYOME ITS release v2.1 has 206,399 distinct taxonomies (DADA2 file from eukaryome.org). The UNITE general release for all eukaryotes, with singletons (2025-02-19), has 88,989.

### Cause

In `src/taxonomy.cpp`, `get_best_genus()`:

```cpp
int get_best_genus(int *karray, float *out_logp, unsigned int arraylen, unsigned int n_kmers, unsigned int ngenus, float *lgk_probability) {
  ...
  int kmer, g, max_g = -1;
  ...
  for(g=0;g<ngenus;g++) {
    lgk_v = &lgk_probability[g*n_kmers];
```

`g` is an `int` and `n_kmers` an `unsigned int` (4^8 = 65,536). The offset `g*n_kmers` is therefore computed in 32 bits, and it wraps once `g >= 65,536` (65,536 × 65,536 = 2^32).

The matrix itself is allocated and filled correctly with `size_t` indices in `C_assign_taxonomy2()` (lines 239, 250, 266). Only this read is wrong. Taxonomies `g`, `g + 65536`, `g + 131072`, … are all scored with the profile of taxonomy `g mod 65536`. They tie exactly, and the tie is broken at random.

### Minimal reproducible example (~17 GB RAM)

```r
library(dada2)
set.seed(1)
random_seq <- function(n, len = 300) {
  vapply(seq_len(n), \(i) paste(sample(c("A", "C", "G", "T"), len, replace = TRUE), collapse = ""), character(1))
}
target <- random_seq(1)
ref <- Biostrings::DNAStringSet(c(random_seq(65536), target))
names(ref) <- c(sprintf("Dummy;D%05d;", 1:65536), "Real;Target;")
ref_file <- tempfile(fileext = ".fasta")
Biostrings::writeXStringSet(ref, ref_file)

res <- assignTaxonomy(target, ref_file, taxLevels = c("Kingdom", "Genus"),
                      minBoot = 0, outputBootstraps = TRUE)
data.frame(taxonomy = res$tax[1, ], bootstrap = res$boot[1, ])
```

The query is identical to the 65,537th reference, so the expected answer is `Real / Target` with bootstrap 100. With dada2 1.40.0 we get a random `Dummy` taxonomy instead: `Dummy / D11528` in our run, with a genus bootstrap of 2. With the one-line fix below, we get `Real / Target` with bootstrap 100.

### Effect on real references

We measured this on 195 ASVs of a fungal mock community (ITS1, Pauvert et al. 2019), `minBoot = 50`:

| reference | distinct taxonomies | ASVs given the correct genus, 1.40.0 → patched |
|---|---|---|
| UNITE all eukaryotes 2025-02-19, without singletons | 64,976 (below the limit) | 90 % (unchanged) |
| UNITE all eukaryotes 2025-02-19, with singletons | 88,989 | 54 % → 90 % |
| EUKARYOME ITS v1.9 (DADA2 file from eukaryome.org) | 142,793 | 24 % → 84 % |
| EUKARYOME ITS v2.1 | 206,399 | 0 % (median Kingdom bootstrap 28) → not rerun yet |

Run time and memory are unchanged. Because the aliased taxonomies tie exactly, the effect also makes results non-reproducible: the tie-breaking RNG is seeded from `std::random_device`, which `set.seed()` does not reach.

### Proposed fix

```diff
--- a/src/taxonomy.cpp
+++ b/src/taxonomy.cpp
@@ -82,7 +82,7 @@
   std::uniform_real_distribution<> cunif(0.0, 1.0);
   
   for(g=0;g<ngenus;g++) {
-    lgk_v = &lgk_probability[g*n_kmers];
+    lgk_v = &lgk_probability[(size_t) g * n_kmers];
     logp = 0.0;
```

Declaring `g` and `n_kmers` as `size_t` in `get_best_genus()` would work as well.

Possibly related: #1032. It reports mostly NA beyond phylum with a large custom training set (SILVA clustered at 99 %), appearing after v1.12, when exact ties between best genera started to be broken at random. That is what this overflow produces: the aliased taxonomies tie exactly, and since the random tie-break their bootstraps collapse. The issue does not give the number of distinct taxonomies, so I cannot confirm it.

I am happy to open a pull request with the one-line fix above if that helps.
