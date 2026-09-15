# Reference databases: sources, derivation, updates

How every reference FASTA used by the benchmark is obtained and built, and
what to change when a new UNITE or EUKARYOME release comes out. Code:
`make_databases.R`. Configuration: the "Reference databases" block of
`config.R`. Tests: `tests/test_make_databases_helpers.R`,
`tests/test_values_map.R`.

## 1. Overview

```
provider release (UNITE .tgz, EUKARYOME .zip)
        │  step 1  download_reference_sources()
        │          dbpq::download_unite_db(url, extract = TRUE)
        │          dbpq::download_eukaryome_db(url) + unzip
        ▼
data/data_raw/refseq/sources/<source>.fasta          general FASTA, k__/p__ headers
data/data_raw/refseq/sources/manifest.csv            provenance (url, doi, md5, date)
        │  step 2  derive_all_variants()
        ├─► dada2_format/<source>.fasta              derive_dada2()   (dada2 methods)
        ├─► sintax_format/<source>.fasta             derive_sintax()  (sintax, lca, blastn)
        │        ├─► <source>_Fungi.fasta            derive_kingdom_only(), kingdom == Fungi
        │        │        └─► <source>_Fungi_cut.fasta   derive_cutadapted(), config primers, untrimmed records kept
        │        └─► mini_<db>.fasta                 derive_mini(), first 10 000 records
        └─► sources/<fake_ref_source>_wo_Fungi.fasta ─► data/data_raw/fake_ref/fake_ref_asv_100.fasta
```

Both steps are idempotent: an existing output is kept and a message says so.
Pass `force = TRUE` to rebuild.

```r
source("make_databases.R")
download_reference_sources()   # network: ~250 MB of archives, ~2.7 GB of FASTA
derive_all_variants()          # long: format conversion of 1.9 M ITS records, cutadapt
```

`make.R` runs both before the pipelines. Sourcing `make.R` also runs the
three production projects (hours); to build a single project, source
`R/run_project.R` instead (§5).

## 2. Sources (config.R::reference_sources)

| source | provider | release | content | DOI / page |
|---|---|---|---|---|
| `Unite_s_all_20250219` | UNITE | 19.02.2025 | general FASTA release for eukaryotes **2**: all eukaryotes, dynamic, **with** singletons (`sh_general_release_dynamic_s_all_19.02.2025.fasta`) | 10.15156/BIO/3301232 |
| `Unite_all_20250219` | UNITE | 19.02.2025 | general FASTA release for eukaryotes: all eukaryotes, dynamic, **without** singletons (`sh_general_release_dynamic_all_19.02.2025.fasta`) | 10.15156/BIO/3301231 |
| `EUK_ITS_v2.1` | EUKARYOME | 2.1 (2026-09-04) | general FASTA, ITS | https://eukaryome.org/generalfasta/ |

`EUK_SSU_v2.1` (EUKARYOME general FASTA, SSU) was a source until 2026-09-15:
ITS1 is not part of the 18S gene, so the SSU release is no longer benchmarked
(ROADMAP decision 7 revised, S6.9). Its files stay on disk.

The exact download URLs are in `config.R`. The UNITE archives also contain a
`_dev` FASTA (sequences with flanking regions); it is ignored.

Before 2026-09-11 the benchmark used `data/data_raw/refseq/Unite.fasta`
(byte-identical to `Unite_s_all_20250219`), `Unite_RefS.fasta` (identical to
`Unite_all_20250219`, despite its name) and EUKARYOME v2.0. Those files are
deleted on 2026-09-11 (see §7).

## 3. Benchmarked databases (config.R::benchmark_dbs)

| db | source | simplification |
|---|---|---|
| `Unite_s_all_20250219` | `Unite_s_all_20250219` | full |
| `Unite_s_all_20250219_Fungi` | `Unite_s_all_20250219` | Fungi |
| `Unite_all_20250219` | `Unite_all_20250219` | full |
| `Unite_all_20250219_Fungi` | `Unite_all_20250219` | Fungi |
| `EUK_ITS_v2.1` | `EUK_ITS_v2.1` | full |
| `EUK_ITS_v2.1_Fungi` | `EUK_ITS_v2.1` | Fungi |
| `EUK_ITS_v2.1_Fungi_cut` | `EUK_ITS_v2.1` | Fungi+cut |

The `Fungi+cut` variant is trimmed by `dbpq::cutadapt_rm_primers_db()` at
ITS1F when found and at the reverse complement of ITS2 when found, with a
minimum overlap of `config.R::primer_min_overlap` (20 bp). Records without
primer sites are kept whole (`config.R::cut_discard_untrimmed = FALSE`): about
94 % of the Fungi records are trimmed, to a median of about 210 bp. Records
shorter than `config.R::cut_min_length` (50 bp) after trimming are dropped
(cutadapt `-m`): ITS2-only records that start at the 5.8S site are left empty
or nearly (317 empty and 1530 shorter than 50 bp in `EUK_ITS_v2.1_Fungi_cut`). Before
2026-09-15 the reverse primer was not reverse-complemented and untrimmed
records were discarded, which kept 4.8 % of the records (ROADMAP S6.9).

`db_list` and `db_meta` (used by `R/values_map.R`, the pipelines and the
chapters) are computed from this table. The seed taxonomy of the DADA2
pipeline (`seed_taxonomy_db`), the release feeding the negative controls
(`fake_ref_source`) and the `preference` consensus database
(`preference_db`) are set right below it.

Names: `<source>` carries the release; `_Fungi` and `_Fungi_cut` mark the
simplification; `mini_` prefixes the smoke-test subsets. Target names are
`<method>__<db>___<parameters>`, so a new release produces new targets and
new store objects, and results of two releases can be compared.

## 4. Header formats at each stage

| stage | example header |
|---|---|
| UNITE general | `>Fusarium_sp\|KF1\|SH1.10FU\|refs\|k__Fungi;p__Ascomycota;…;s__Fusarium_sp` |
| EUKARYOME general | `>EUK0000001;k__Fungi;p__Ascomycota;…;s__unclassified` |
| dada2 format | `>Fungi;Ascomycota;…;Fusarium_sp;` |
| sintax format | `>Fusarium_sp\|KF1\|SH1.10FU\|refs;tax=k:Fungi,p:Ascomycota,…,s:Fusarium_sp` |

Check with `head -1` after any rebuild: the kingdom must be the first rank in
both converted formats.

## 5. Moving to a new release

1. **Find the download URL.**
   - UNITE: https://unite.ut.ee/repository.php lists the general FASTA
     releases and their DOIs. The DOIs resolve to an HTML page, so ask the
     public PlutoF API for the archive URL:
     ```sh
     curl -s "https://api.plutof.ut.ee/v1/public/dois/?identifier=10.15156/BIO/3301232" \
       | grep -o '"url":"[^"]*"'
     ```
     "eukaryotes 2" / "Fungi 2" are the releases with singletons.
   - EUKARYOME: copy the `General_EUK_<marker>_v<version>.zip` links from
     https://eukaryome.org/generalfasta/.
2. **Edit `config.R`.** Add or replace rows of `reference_sources` with a new
   versioned `source` name, then update `benchmark_dbs`, and if needed
   `seed_taxonomy_db`, `fake_ref_source` and `preference_db`. Keep the old
   rows to benchmark both releases side by side.
3. **Check the configuration:** `Rscript tests/test_values_map.R`.
4. **Build:** `download_reference_sources()` then `derive_all_variants()`.
   Inspect `sources/manifest.csv` and the first header of each new file (§4).
5. **Smoke test.** Do not `source("make.R")` for this: it runs the
   production projects too. `R/run_project.R` only defines `run_project()`:
   ```r
   source("R/run_project.R")
   run_project("assign_taxo_mini")   # 12 methods x length(db_list) targets
   run_project("cross_val_mini")     # 2 folds, 200 sequences (config.R, mini_db)
   ```
   Check `targets::tar_errored()` with the same `TAR_PROJECT`, or follow the
   run with `targets::tar_poll()` in a second session.
6. **Production:** `run_project("dada2")`, `run_project("assign_taxo")`,
   `run_project("cross_val")` (or `source("make.R")`, which runs step 4 and
   these three in order), then `quarto render analysis`. Restore the
   publication values of the CV knobs in `config.R` first.
   `tar_prune()` in `run_project()` removes the targets of databases no longer
   in `benchmark_dbs`.

Nothing else hard-codes a database name: the chapters, the pipelines and the
tests read `config.R`.

## 6. Pitfalls handled by the code

- **UNITE separator.** In UNITE general headers a `|` precedes `k__Kingdom`.
  `dbpq::format2sintax()` and `format2dada2()` then take `…|k__Kingdom` as the
  identifier and drop the kingdom, shifting every rank by one in both output
  formats. `derive_sintax()` and `derive_dada2()` rewrite `|k__` into `;k__`
  first (`general_headers_fixed()`; ROADMAP S6.2).
- **EUKARYOME name qualifiers.** EUKARYOME disambiguates homonyms with a
  kingdom or an order (`g__Lactarius(Fungi)`,
  `f__Gonostomatidae(Sporadotrichida)`), brackets some genera
  (`g__(Candida)`, sometimes `g__(Candida]`) and adds `.s.str`, `.s.str.`,
  `.nom.prov` or `.nom.provis`. v2.1 has 11 413 ITS and 5 586 SSU headers
  with parentheses (v2.0 already had 9 641 in ITS). Left as is, the sintax
  prediction `g:Lactarius(Fungi)(0.97)` holds two `(`:
  `MiscMetabar::assign_sintax()` reads its bootstrap as NA, so
  `min_bootstrap` never filters that genus; dada2 and lca return
  `Lactarius(Fungi)`, which never matches the mock (`Lactarius`,
  `Cryptococcus`, `Mortierella` and `Peziza` are affected).
  `general_headers_fixed()` keeps the name only (`g__Lactarius`,
  `g__Candida`, `g__Mortierella`) before the conversion. Placeholders such as
  `Dothideales.gen05` or `.fam.incertae.sedis` are kept (ROADMAP S7.8).
- **Rebuilding a file in place invalidates exactly its targets.** Each
  reference fasta is a file target of `pipelines/assign_taxo.R` and
  `pipelines/cross_val.R` (`ref_dada2__<db>`, `ref_sintax__<db>`; column
  `ref_file` of `values_map`; ROADMAP S7.9). After
  `derive_all_variants(force = TRUE)` under unchanged names, the next
  `tar_make()` reruns the assignment and CV targets that read a changed file,
  and only those.
- **Fungi filter.** A plain `Fungi` pattern also keeps `f:Fungiidae` (corals),
  `k:cf.Fungi` or species names such as `fungiformis`: 965 non-Fungi records
  in the former EUKARYOME ITS v2 Fungi file, 22 in SSU, 19 in UNITE.
  `derive_kingdom_only()` anchors the pattern on the kingdom field of each
  format (`kingdom_pattern()`).
- **Nested EUKARYOME archives.** `General_EUK_<marker>_v2.1.zip` contains a
  `.7z` archive, which contains the FASTA (`General_EUK_SSU_v2.1.fasta`:
  867 MB, 558 632 records, against 353 797 in v2.0). `extract_single_fasta()`
  unzips, then extracts the `.7z` with the `7z` command-line tool (p7zip),
  which must be on `PATH`.
- **Fake reference kept by idempotence.** `derive_all_variants()` keeps an
  existing `data/data_raw/fake_ref/fake_ref_asv_100.fasta` like any other
  output, even when `fake_ref_source` or the non-Fungi pool changed. Rebuild
  it explicitly:
  ```r
  source("make_databases.R")
  derive_fake_ref(
    file.path(sources_dir, paste0(fake_ref_source, "_wo_Fungi.fasta")),
    here(fake_ref_fasta), seed = targets_seed, force = TRUE
  )
  ```
  `pipelines/assign_taxo.R` tracks it as the file target `fake_ref_file`, so
  the next `tar_make()` rebuilds `d_asv_for_assignation` and every assignment.
  Regenerated on 2026-09-11 from `Unite_s_all_20250219_wo_Fungi.fasta`
  (100 records, 73 phyla, no Fungi; 80 sequences shared with the May 2026
  file, which came from the legacy `Fungi` pattern).
- **UNITE DOIs** resolve to an HTML page; `dbpq::download_unite_db()` takes a
  direct `url` instead (added to dbpq on 2026-09-11).
- **Long steps.** Converting the EUKARYOME ITS release (≈1.6 M records) with
  dbpq takes tens of minutes per format; primer trimming needs cutadapt in the
  conda env of `config.R::cutadapt_conda_prelude`.

## 7. Legacy files

Deleted on 2026-09-11 (40 files, 14.3 GB): the unversioned
`data/data_raw/refseq/{Unite,Unite_RefS,Euk_ITS_v2,Euk_SSU_v2,Euk_ITS_v2_1}.fasta`
and the unversioned files of `dada2_format/` and `sintax_format/`
(`Unite*.fasta`, `UNITE_RefS.fasta`, `EUK_ITS_v2*.fasta`, `EUK_SSU_v2*.fasta`
and their `mini_*`). Their UNITE sintax headers lost the kingdom and their
`_Fungi` variants contained non-Fungi records. The EUKARYOME v2.0 releases
are no longer on disk; UNITE 19.02.2025 and EUKARYOME v2.1 are in `sources/`.
