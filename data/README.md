# data/

Only `data/data_raw/metadata/*.csv` is tracked in git; everything else is
regenerated (`make_databases.R`, the pipelines) or downloaded. The reference
databases are documented in detail in `docs/reference_databases.md`.

## data_raw/

| Folder / file | Content | Provenance | Used by |
|---|---|---|---|
| `metadata/sam_data.csv` | Sample metadata of the mock community (3 samples `ITSEnz_Dik{A,B,C}`) | Pauvert et al. 2019 | `pipelines/dada2.R` (`sam_data_matching_names`) |
| `metadata/taxo_mock.csv` | Truth table of the mock community (Kingdom → Species, `MockStrain`) | Pauvert et al. 2019 | `pipelines/assign_taxo.R` (`taxo_mock`), `tc_metrics_mock()` |
| `metadata/taxo_mock_SRR30413326.csv` | Truth table of a second fungal mock (16 taxa: Russula, Lactifluus, Sydowia, …) | SRA run SRR30413326 | planned dataset D1d (ROADMAP) |
| `rawseq/*.fastq` | Paired-end ITS1 reads of the mock community | Pauvert et al. 2019 | `pipelines/dada2.R` |
| `rawseq/SRR30413326/` | Empty; destination of the SRR30413326 fastqs (`fasterq-dump SRR30413326`) | SRA | D1d |
| `refseq/sources/<source>.fasta` | General FASTA release of each row of `config.R::reference_sources` (`Unite_s_all_20250219`, `Unite_all_20250219`, `EUK_ITS_v2.1`; `EUK_SSU_v2.1` until 2026-09-15, its files kept) | `make_databases.R::download_reference_sources()` via `dbpq::download_unite_db()` / `dbpq::download_eukaryome_db()` | step 2 of `make_databases.R` |
| `refseq/sources/manifest.csv`, `<source>.provenance.csv` | One row per source: provider, release, DOI, URL, original file name, bytes, md5, download date | idem | provenance for the manuscript |
| `refseq/sources/<fake_ref_source>_wo_Fungi.fasta` | Non-Fungi records of the release used for the negative controls | `derive_no_pattern()` | `derive_fake_ref()` |
| `refseq/dada2_format/<db>.fasta`, `mini_<db>.fasta` | Every source and benchmarked database in dada2 header format (`>Fungi;Ascomycota;…;`), plus the first 10 000 records | `derive_all_variants()` | `dada2` method, `pipelines/dada2.R` seed taxonomy |
| `refseq/sintax_format/<db>.fasta`, `mini_<db>.fasta` | Same in sintax header format (`>ID;tax=k:Fungi,p:…`) | `derive_all_variants()` | `sintax`, `lca`, `blastn` methods |
| `fake_ref/fake_ref_asv_100.fasta` | 100 non-Fungi UNITE records, one per phylum then random fill (seed `targets_seed`) | `make_databases.R::derive_fake_ref()` | `pipelines/assign_taxo.R` (`add_external_seq_pq`, TN denominator) |
| `mock_hleap2021/Mocks/` | Five mock communities of Hleap et al. 2021 (`*_realized.fa` + taxonomy `.txt`; 32–387 sequences each). **Fish** (Chordata) 12S/COI mocks, not fungal ITS. | Hleap et al. 2021 supplementary data | nothing yet; decision pending (ROADMAP D1d) |

### Legacy reference files (deleted 2026-09-11)

40 files, 14.3 GB (ROADMAP S7.5): `refseq/{Unite,Unite_RefS,Euk_ITS_v2,Euk_SSU_v2,Euk_ITS_v2_1}.fasta` and the unversioned files of `refseq/{dada2,sintax}_format/` (`Unite`, `Unite_Fungi`, `Unite_wo_fungi`, `UNITE_RefS`, `EUK_ITS_v2*`, `EUK_SSU_v2*` and their `mini_*`). `Unite.fasta` and `Unite_RefS.fasta` were byte-identical to `sources/Unite_s_all_20250219.fasta` and `sources/Unite_all_20250219.fasta` (the `RefS` name was misleading), `Euk_ITS_v2_1.fasta` to `sources/EUK_ITS_v2.1.fasta`; the EUKARYOME v2.0 releases are no longer on disk. The derived files lost the UNITE kingdom in sintax headers and let non-Fungi records into the `_Fungi` files (ROADMAP S6.2, S7.2).

## data_intermediate/ (~160 MB)

Written by `pipelines/dada2.R`: `seq_wo_primers/` (cutadapt output + json
reports), `filterAndTrim_fwd/`, `filterAndTrim_rev/`. Safe to delete; the
pipeline recreates them.

## data_final/

| Path | Content |
|---|---|
| `autometric/assign_taxo/<full_name>__<timestamp>.txt` | One {autometric} log per assignment target and per run, written inside the crew worker (`R/autometric_helpers.R`). Aggregated by the `benchmark_costs` target (newest file per target). Empty files belong to targets shorter than the 1 s sampling interval. |
| `autometric/dada2/<phase>__<timestamp>.txt` | Same for the DADA2 phases; aggregated by `benchmark_costs_dada2`. |
| `autometric_log_assign_taxo.txt`, `autometric_log_dada2.txt` | Legacy single-file logs (2025-02 to 2026-05). No longer written; kept only for `analysis/sandbox/autometric_regex_costs.qmd`. |
