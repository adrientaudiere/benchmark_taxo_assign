# From the raw exact hits of the species project (pipelines/species.R,
# R/species_exact.R) to one species call per query, and to tax_table columns
# that chapter 01's scoring can read (analysis/09_species_exact.qmd, OD 43).
#
# The rule (developer, 2026-10-07, OD 43): a query is named a species when
# every exact hit that names a species gives the same binomial; hits naming
# no species do not block the call; several binomials give no species.
# assignSpecies() follows it by construction: its reference holds only the
# records that name a species.
#
# Requires dplyr, tibble and tidyr.

# Status of a query against one database for one method:
#   "one species"       every named hit gives the same binomial: the call;
#   "several species"   the named hits disagree: no call;
#   "only unnamed hits" hits, none naming a species: no call
#                       (vsearch and blastn; assignSpecies cannot see them);
#   "no hit"            no hit at that match level.
species_status_levels <- c(
  "one species",
  "several species",
  "only unnamed hits",
  "no hit"
)

# Match levels of the near-match hits (R/species_exact.R, near_hit_columns()),
# each with the tag of its scored column: "exact" (no difference over the
# query), "1 difference" (at most one), "90 % covered" (100 % identity over at
# least 90 % of the query, the rest beyond an end of the record: the record
# stops early; developer 2026-10-07).
match_tags <- c(
  "exact" = "exact",
  "1 difference" = "diff1",
  "90 % covered" = "cover90"
)

# The hits of `hits` (with `diffs` and `exact_cover`) at one match level.
hits_at <- function(hits, match) {
  switch(
    match,
    "exact" = hits[hits$diffs == 0, ],
    "1 difference" = hits[hits$diffs <= 1, ],
    "90 % covered" = hits[hits$exact_cover, ],
    stop("Unknown match level: ", match)
  )
}

# Calls of one hit table with a `binomial` column (vsearch, blastn): one row
# per (db, query) that has at least one exact hit.
calls_from_hits <- function(hits, method) {
  hits |>
    dplyr::distinct(.data$db, .data$query, .data$binomial) |>
    dplyr::group_by(.data$db, .data$query) |>
    dplyr::summarise(
      n_binomials = dplyr::n_distinct(.data$binomial, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      hits |>
        dplyr::filter(!is.na(.data$binomial)) |>
        dplyr::distinct(.data$db, .data$query, .data$binomial) |>
        dplyr::group_by(.data$db, .data$query) |>
        dplyr::filter(dplyr::n() == 1) |>
        dplyr::ungroup(),
      by = c("db", "query")
    ) |>
    dplyr::mutate(
      method = method,
      status = dplyr::case_when(
        .data$n_binomials == 1 ~ "one species",
        .data$n_binomials > 1 ~ "several species",
        TRUE ~ "only unnamed hits"
      ),
      call = dplyr::if_else(
        .data$status == "one species",
        .data$binomial,
        NA_character_
      )
    ) |>
    dplyr::select("db", "query", "method", "status", "call", "n_binomials")
}

# Calls of dada2::assignSpecies() (assign_species_exact()): Genus is NA when
# the hits span several genera, and the epithets of several species are
# joined by "/".
calls_from_assign_species <- function(assigned) {
  assigned |>
    dplyr::filter(!is.na(.data$genus) | !is.na(.data$epithets)) |>
    dplyr::mutate(
      method = "assignSpecies",
      single = !is.na(.data$genus) &
        !is.na(.data$epithets) &
        !grepl("/", .data$epithets, fixed = TRUE),
      status = dplyr::if_else(.data$single, "one species", "several species"),
      call = dplyr::if_else(
        .data$single,
        paste0(.data$genus, "_", .data$epithets),
        NA_character_
      ),
      n_binomials = dplyr::if_else(.data$single, 1L, NA_integer_)
    ) |>
    dplyr::select("db", "query", "method", "status", "call", "n_binomials")
}

# Every (db, query, method) of the three methods with its status and call at
# exact match; `queries` is the species_queries$seqs names (all queries
# searched), so that the queries without an exact hit get the status "no hit".
# `blastn_hits` are the exact hits (hits_at(<blastn near hits>, "exact")).
species_calls <- function(assigned, vsearch_hits, blastn_hits, queries) {
  found <- dplyr::bind_rows(
    calls_from_assign_species(assigned),
    calls_from_hits(vsearch_hits, "vsearch"),
    calls_from_hits(blastn_hits, "blastn")
  )
  grid <- tidyr::expand_grid(
    db = unique(c(assigned$db, vsearch_hits$db, blastn_hits$db)),
    query = queries,
    method = c("assignSpecies", "vsearch", "blastn")
  )
  grid |>
    dplyr::left_join(found, by = c("db", "query", "method")) |>
    dplyr::mutate(
      status = factor(
        dplyr::coalesce(.data$status, "no hit"),
        levels = species_status_levels
      )
    )
}

# Every (db, query, method, match) of vsearch and blastn at the three match
# levels, from their near hits, with its status and call (same rule as at
# exact match).
near_calls <- function(vsearch_near, blastn_near, queries) {
  found <- dplyr::bind_rows(lapply(names(match_tags), function(m) {
    dplyr::bind_rows(
      calls_from_hits(hits_at(vsearch_near, m), "vsearch"),
      calls_from_hits(hits_at(blastn_near, m), "blastn")
    ) |>
      dplyr::mutate(match = m)
  }))
  tidyr::expand_grid(
    db = unique(c(vsearch_near$db, blastn_near$db)),
    query = queries,
    method = c("vsearch", "blastn"),
    match = names(match_tags)
  ) |>
    dplyr::left_join(found, by = c("db", "query", "method", "match")) |>
    dplyr::mutate(
      status = factor(
        dplyr::coalesce(.data$status, "no hit"),
        levels = species_status_levels
      ),
      match = factor(.data$match, levels = names(match_tags))
    )
}

# Column name of one species assignment in the scored object, read by
# assignment_db() (R/score_assignments.R) as database `db`; `match` is a
# level of match_tags.
exact_column <- function(input, method, db, match = "exact") {
  paste0(input, "_", method, "__", db, "___", unname(match_tags[match]))
}

# `physeq` with one "Gen_sp_<column>" tax_table column per (method, db) of
# `calls` (per (method, db, match) when `calls` has a `match` column) for the
# taxa of one input, plus an all-NA "no_rank" column, and the
# ranks_df that score_assignments() takes for them: tc_metrics_unit() needs
# every rank of truth_ranks (it places a unit's truth depth among them), so
# the six upper ranks point to "no_rank" and only the Species rows of the
# result are meaningful. `map` is species_queries$map restricted to the input.
add_exact_columns <- function(physeq, calls, map, input) {
  if (!"match" %in% names(calls)) {
    calls$match <- "exact"
  }
  calls$match <- as.character(calls$match)
  cols <- calls |>
    dplyr::distinct(.data$method, .data$db, .data$match) |>
    dplyr::mutate(
      column = exact_column(input, .data$method, .data$db, .data$match)
    )
  per_taxon <- calls |>
    dplyr::filter(!is.na(.data$call)) |>
    dplyr::inner_join(map, by = "query", relationship = "many-to-many") |>
    dplyr::mutate(
      column = exact_column(input, .data$method, .data$db, .data$match)
    )
  taxa <- phyloseq::taxa_names(physeq)
  added <- matrix(
    NA_character_,
    nrow = length(taxa),
    ncol = nrow(cols),
    dimnames = list(taxa, paste0("Gen_sp_", cols$column))
  )
  rows <- match(per_taxon$taxon, taxa)
  keep <- !is.na(rows)
  added[cbind(
    rows[keep],
    match(paste0("Gen_sp_", per_taxon$column[keep]), colnames(added))
  )] <- per_taxon$call[keep]
  physeq@tax_table <- phyloseq::tax_table(cbind(
    methods::as(phyloseq::tax_table(physeq), "matrix"),
    added,
    no_rank = NA_character_
  ))
  ranks_df <- as.data.frame(stats::setNames(
    lapply(cols$column, \(cl) c(rep("no_rank", 6), paste0("Gen_sp_", cl))),
    cols$column
  ))
  list(physeq = physeq, ranks_df = ranks_df)
}
