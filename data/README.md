# data/

Only `data/data_raw/metadata/*.csv` is tracked in git; everything else is
regenerated (`make_databases.R`, the pipelines) or downloaded. The reference
databases are documented in detail in `docs/reference_databases.md`.

## data_raw/

| Folder / file | Content | Provenance | Used by |
|---|---|---|---|
| `metadata/sam_data.csv` | Sample metadata of the mock community (3 samples `ITSEnz_Dik{A,B,C}`) | Pauvert et al. 2019 | `pipelines/dada2.R` (`sam_data_matching_names`) |
| `metadata/taxo_mock.csv` | Truth table of the mock community (Kingdom → Species, `MockStrain`) | Pauvert et al. 2019 | `pipelines/assign_taxo.R` (`taxo_mock`), `tc_metrics_mock()` |
| `metadata/pauvert2019_table_s1_mock_sanger.csv` | Table S1 of Pauvert et al. 2019: the 189 mock strains with their ITS Sanger sequence (`Sanger_Sequence`, median 542 bp), taxonomy, `Unique_ITS1_Sanger_Sequence` (175 Yes) and `Present_In_Raw_MiSeq_Data` (160 Yes); `;`-separated. Downloaded 2026-09-17 from `https://ars.els-cdn.com/content/image/1-s2.0-S1754504818302800-mmc2.csv`, the file named in the ADDENDUM of the authors' code deposit (doi:10.15454/VKTWKR) | Pauvert et al. 2019, Fungal Ecology, supplementary mmc2 | planned per-ASV truth (`docs/hleap_2021_metrics.md` §4) |
| `metadata/tedersoo_mock_table_s1_sanger.csv` | Table S1 of the Tedersoo mock (doi:10.1111/1755-0998.70189): 103 specimens (`Specimen`, `Identification` = binomial only), `Sanger_read` (upper-cased, 466–953 bp) and read counts per platform × pipeline × PCR cycles. Written from the xlsx below | Tedersoo et al., supplementary Table S1 | planned per-ASV truth (`docs/hleap_2021_metrics.md` §6) |
| `metadata/tedersoo_mock_table_s1.xlsx` | The original supplementary file (`men70189-sup-0002-tables1.xlsx`, md5 `d9a4fcb8b295534c43ad586902f7bd8a`), downloaded by hand on 2026-09-17 (Wiley refuses scripted downloads). Not tracked by git | Tedersoo et al. | source of the csv above |
| `metadata/pauvert_gna_lineage.csv` | GBIF lineage (Kingdom → Genus) and accepted name of the 181 Pauvert species names, from the GNA verifier (`R/mock_truth.R::refresh_gna_lineage_cache()`, 2026-09-22). Read by the per-unit truth builder, which never calls the network | GNA verifier, GBIF backbone | `R/mock_truth.R` |
| `metadata/tedersoo_illumina_gna_lineage.csv` | Same, for the 103 Tedersoo names (their Table S1 gives a binomial only) | GNA verifier, GBIF backbone | `R/mock_truth.R` |
| `metadata/taxo_mock_SRR30413326.csv` | Truth table of a second fungal mock (16 taxa: Russula, Lactifluus, Sydowia, …) | SRA run SRR30413326 | planned dataset D1d (ROADMAP) |
| `metadata/sam_data_tedersoo_illumina.csv` | Sample metadata of the Tedersoo Illumina mock (2 runs, columns `Sample_names` / `platform` / `community`; the second variable exists because `add_shuffle_seq_pq()` needs at least two) | written from the ENA run report of PRJEB108994 | `pipelines/dada2_bio.R` (`sam_data_matching_names`) |
| `rawseq/*.fastq` | Paired-end ITS1 reads of the mock community | Pauvert et al. 2019 | `pipelines/dada2.R` |
| `rawseq/SRR30413326/` | Empty; destination of the SRR30413326 fastqs (`fasterq-dump SRR30413326`) | SRA | D1d |
| `rawseq_tedersoo_illumina/*.fastq.gz` | Paired-end full-length ITS reads of the Tedersoo Illumina mock: `ERR16773638` (`ILLUvs__MOCK`) and `ERR16773639` (`ILLUvs__MOCK_2`), iSeq 100, ~495 bp per read, 19.7 MB. `ena_md5.txt` holds the ENA checksums, all four verified on download (2026-09-16). Fetch over **HTTPS** (`https://ftp.sra.ebi.ac.uk/vol1/fastq/…`): EBI's FTP endpoint fails on these paths. | ENA project PRJEB108994, Tedersoo et al. 2026 | `pipelines/dada2_bio.R` (project `dada2_tedersoo_illumina`) |
| `rawseq_*/` | One folder per `config.R::bio_datasets` entry (D1c). Ignored by git as a glob, since the table is meant to grow | ENA / SRA | `pipelines/dada2_bio.R` |
| `refseq/sources/<source>.fasta` | General FASTA release of each row of `config.R::reference_sources` (`Unite_s_all_20250219`, `Unite_all_20250219`, `EUK_ITS_v2.1`; `EUK_SSU_v2.1` until 2026-09-15, its files kept) | `make_databases.R::download_reference_sources()` via `dbpq::download_unite_db()` / `dbpq::download_eukaryome_db()` | step 2 of `make_databases.R` |
| `refseq/sources/manifest.csv`, `<source>.provenance.csv` | One row per source: provider, release, DOI, URL, original file name, bytes, md5, download date | idem | provenance for the manuscript |
| `refseq/sources/<fake_ref_source>_wo_Fungi.fasta` | Non-Fungi records of the release used for the negative controls | `derive_no_pattern()` | `derive_fake_ref()` |
| `refseq/dada2_format/<db>.fasta`, `mini_<db>.fasta` | Every source and benchmarked database in dada2 header format (`>Fungi;Ascomycota;…;`), plus the first 10 000 records | `derive_all_variants()` | `dada2` method, `pipelines/dada2.R` seed taxonomy |
| `refseq/sintax_format/<db>.fasta`, `mini_<db>.fasta` | Same in sintax header format (`>ID;tax=k:Fungi,p:…`) | `derive_all_variants()` | `sintax`, `lca`, `blastn` methods |
| `fake_ref/fake_ref_asv_100.fasta` | 100 non-Fungi UNITE records from the 13 kingdoms retained in every release (`config.R::rep_kingdom_min_records`, since 2026-09-17), one per phylum then random fill (seed `targets_seed`) | `make_databases.R::derive_fake_ref()` | `pipelines/assign_taxo.R` (`add_external_seq_pq`, TN denominator) |

### Legacy reference files (deleted 2026-09-11)

40 files, 14.3 GB (ROADMAP S7.5): `refseq/{Unite,Unite_RefS,Euk_ITS_v2,Euk_SSU_v2,Euk_ITS_v2_1}.fasta` and the unversioned files of `refseq/{dada2,sintax}_format/` (`Unite`, `Unite_Fungi`, `Unite_wo_fungi`, `UNITE_RefS`, `EUK_ITS_v2*`, `EUK_SSU_v2*` and their `mini_*`). `Unite.fasta` and `Unite_RefS.fasta` were byte-identical to `sources/Unite_s_all_20250219.fasta` and `sources/Unite_all_20250219.fasta` (the `RefS` name was misleading), `Euk_ITS_v2_1.fasta` to `sources/EUK_ITS_v2.1.fasta`; the EUKARYOME v2.0 releases are no longer on disk. The derived files lost the UNITE kingdom in sintax headers and let non-Fungi records into the `_Fungi` files (ROADMAP S6.2, S7.2).

## data_intermediate/ (~160 MB)

Written by `pipelines/dada2.R`: `seq_wo_primers/` (cutadapt output + json
reports), `filterAndTrim_fwd/`, `filterAndTrim_rev/`. Safe to delete; the
pipeline recreates them.

`pipelines/dada2_bio.R` writes the same three folders suffixed with the dataset
name (`seq_wo_primers_tedersoo_illumina/`, `filterAndTrim_fwd_tedersoo_illumina/`,
`filterAndTrim_rev_tedersoo_illumina/`), so a biological dataset never
overwrites the mock-community intermediates.

## data_final/

| Path | Content |
|---|---|
| `autometric/assign_taxo/<full_name>__<timestamp>.txt` | One {autometric} log per assignment target and per run, written inside the crew worker (`R/autometric_helpers.R`). Aggregated by the `benchmark_costs` target (newest file per target). Empty files belong to targets shorter than the 1 s sampling interval. |
| `autometric/dada2/<phase>__<timestamp>.txt` | Same for the DADA2 phases; aggregated by `benchmark_costs_dada2`. |
| `autometric/dada2_<dataset>/<phase>__<timestamp>.txt` | Same again for each D1c biological dataset (`pipelines/dada2_bio.R`), in a folder of its own so a dataset run never overwrites the mock-community costs. |
| `autometric_log_assign_taxo.txt`, `autometric_log_dada2.txt` | Legacy single-file logs (2025-02 to 2026-05). No longer written; kept only for `analysis/sandbox/autometric_regex_costs.qmd`. |
