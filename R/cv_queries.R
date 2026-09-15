# Cross-validation queries trimmed to the amplicon region (ROADMAP D1a,
# decided 2026-09-15): the test sequences are cut with cutadapt like the reads
# of an Illumina run, the training part keeps the full-length records.
#
# Requires Biostrings and the dbpq checkout (load_pqverse("dbpq")).

# Trim the records of `dna` (a named DNAStringSet) to the region between the
# primers with dbpq::cutadapt_rm_primers_db(): the forward primer is trimmed
# when found, the reverse-primer site (its reverse complement) is required, so
# records without a complete site are dropped. `prelude` activates cutadapt
# (NULL: the dbpq default). Returns the trimmed records, names unchanged.
trim_cv_queries <- function(
  dna,
  primer_fw,
  primer_rev,
  min_overlap = nchar(primer_rev),
  nproc = 1,
  prelude = NULL
) {
  input <- tempfile(pattern = "cv_queries_", fileext = ".fasta")
  output <- tempfile(pattern = "cv_queries_trimmed_", fileext = ".fasta")
  on.exit(unlink(c(input, output)))
  Biostrings::writeXStringSet(dna, input)
  args <- list(
    ref_fasta = input,
    output = output,
    primer_fw = primer_fw,
    primer_rev = primer_rev,
    rc_primer_rev = TRUE,
    primer_rev_required = TRUE,
    discard_untrimmed = TRUE,
    min_overlap = min_overlap,
    nproc = nproc,
    verbose = FALSE
  )
  if (!is.null(prelude)) {
    args$args_before_cutadapt <- prelude
  }
  do.call(dbpq::cutadapt_rm_primers_db, args)
  Biostrings::readDNAStringSet(output)
}

# Identifiers of the records of a sintax-format fasta (the text before
# ";tax="). dada2-format headers are the taxonomy only, so cross_val() names
# the records of a dada2-format file with the identifiers of the sintax-format
# file of the same database (ROADMAP D1a, B22). `widths`, when given, are the
# sequence lengths of the paired file: they must match record by record.
record_ids <- function(fasta, widths = NULL) {
  index <- Biostrings::fasta.index(fasta)
  if (
    !is.null(widths) &&
      !identical(as.integer(index$seqlength), as.integer(widths))
  ) {
    stop(
      fasta,
      " does not hold the same records in the same order as the paired file",
      " (sequence lengths differ)."
    )
  }
  ids <- sub(";tax=.*$", "", index$desc)
  if (anyDuplicated(ids) > 0) {
    stop("Duplicated record identifiers in ", fasta, ".")
  }
  ids
}

# Subsample of a cross-validation. `pool_names` are the drawn records in their
# random order, `query_names` the records usable as queries (primer site
# found). The pool is kept up to its `max_seq`-th query, so every database gets
# max_seq queries whatever its share of records without primer site; the other
# records of the kept pool are training records only. Returns the kept pool
# and query names, in drawn order.
cv_select_queries <- function(pool_names, query_names, max_seq = NULL) {
  is_query <- pool_names %in% query_names
  last <- length(pool_names)
  if (!is.null(max_seq) && sum(is_query) > max_seq) {
    last <- which(cumsum(is_query) == max_seq)[1]
  } else if (!is.null(max_seq) && sum(is_query) < max_seq) {
    message(
      "cv_select_queries(): ",
      sum(is_query),
      " queries only for max_seq = ",
      max_seq,
      "; increase `oversample` in cross_val()."
    )
  }
  pool <- pool_names[seq_len(last)]
  list(pool = pool, queries = pool[is_query[seq_len(last)]])
}
