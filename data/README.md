# data/

Only `data/data_raw/metadata/*.csv` is tracked in git; everything else is
regenerated (`make_databases.R`, the pipelines) or downloaded. Sizes as of
2026-09-10.

## data_raw/ (~13 GB)

| Folder / file | Content | Provenance | Used by |
|---|---|---|---|
| `metadata/sam_data.csv` | Sample metadata of the mock community (3 samples `ITSEnz_Dik{A,B,C}`) | Pauvert et al. 2019 | `pipelines/dada2.R` (`sam_data_matching_names`) |
| `metadata/taxo_mock.csv` | Truth table of the mock community (Kingdom → Species, `MockStrain`) | Pauvert et al. 2019 | `pipelines/assign_taxo.R` (`taxo_mock`), `tc_metrics_mock()` |
| `metadata/taxo_mock_SRR30413326.csv` | Truth table of a second fungal mock (16 taxa: Russula, Lactifluus, Sydowia, …) | SRA run SRR30413326 | planned dataset D1d (ROADMAP) |
| `rawseq/*.fastq` | Paired-end ITS1 reads of the mock community | Pauvert et al. 2019 | `pipelines/dada2.R` |
| `rawseq/SRR30413326/` | Empty; destination of the SRR30413326 fastqs (`fasterq-dump SRR30413326`) | SRA | D1d |
| `refseq/Unite.fasta`, `refseq/Unite_RefS.fasta` | UNITE (all eukaryotes; `_RefS` = RefSeq-only singletons variant, side question in ROADMAP) | UNITE download, dada2/`k__` header format | `make_databases.R` source file |
| `refseq/Euk_ITS_v2.fasta`, `refseq/Euk_SSU_v2.fasta` | EUKARYOME v2 ITS and SSU | EUKARYOME download, `k__` header format | `make_databases.R` source files |
| `refseq/dada2_format/` | Every derived variant in dada2 header format (`Unite`, `Unite_Fungi`, `EUK_ITS_v2[_Fungi[_cut]]`, `EUK_SSU_v2[_Fungi[_cut]]`, `Unite_wo_fungi`, `mini_*`) | `make_databases.R::derive_all_variants()` | `dada2` method (`db_list` in `config.R`) |
| `refseq/sintax_format/` | Same variants in sintax `tax=` header format | idem | `sintax`, `lca`, `blastn` methods |
| `fake_ref/fake_ref_asv_100.fasta` | 100 non-Fungi UNITE records, one per phylum then random fill (seed `targets_seed`) | `make_databases.R::derive_fake_ref()` | `pipelines/assign_taxo.R` (`add_external_seq_pq`, TN denominator) |
| `mock_hleap2021/Mocks/` | Five mock communities of Hleap et al. 2021 (`*_realized.fa` + taxonomy `.txt`; 32–387 sequences each). **Fish** (Chordata) 12S/COI mocks, not fungal ITS. | Hleap et al. 2021 supplementary data | nothing yet; decision pending (ROADMAP D1d) |

## data_intermediate/ (~160 MB)

Written by `pipelines/dada2.R`: `seq_wo_primers/` (cutadapt output + json
reports), `filterAndTrim_fwd/`, `filterAndTrim_rev/`. Safe to delete; the
pipeline recreates them.

## data_final/ (~6 MB)

| Path | Content |
|---|---|
| `autometric/assign_taxo/<full_name>__<timestamp>.txt` | One {autometric} log per assignment target and per run, written inside the crew worker (`R/autometric_helpers.R`). Aggregated by the `benchmark_costs` target (newest file per target). |
| `autometric/dada2/<phase>__<timestamp>.txt` | Same for the DADA2 phases (`cutadapt`, `filtered`, `err_*`, `dd*`, `tax_tab`); aggregated by `benchmark_costs_dada2`. |
| `autometric_log_assign_taxo.txt`, `autometric_log_dada2.txt` | Legacy single-file logs (2025-02 to 2026-05). No longer written; kept only for `analysis/sandbox/autometric_regex_costs.qmd`. Delete once the new logs exist. |
