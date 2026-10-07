# ETAT

Dated figures of the project, newest first. Written by `tools/etat.R`
(`R/etat_helpers.R`), never by hand: README, CONTEXT, AGENTS and
ROADMAP state the decisions and point here for the counts.

## 2026-10-07

### Targets projects

|project                       | targets| built| outdated| errored| orphans| computations| assignment_rows|
|:-----------------------------|-------:|-----:|--------:|-------:|-------:|------------:|---------------:|
|dada2                         |      32|    32|       27|       0|       0|            0|               0|
|assign_taxo                   |     485|   484|      466|       0|     615|           72|             378|
|assign_taxo_mini              |     485|   485|        0|       0|     615|           72|             378|
|cross_val                     |      72|    72|       62|       0|       0|            0|               0|
|cross_val_mini                |      72|    72|        0|       0|       0|            0|               0|
|dada2_tedersoo_illumina       |      30|    30|       29|       0|       0|            0|               0|
|assign_taxo_tedersoo_illumina |    1399|  1398|     1376|       0|       0|          216|            1134|
|dada2_insilico                |      30|    30|       29|       0|       0|            0|               0|
|assign_taxo_insilico          |    1399|  1398|     1376|       0|     615|          216|            1134|
|dada2_rolling_mock            |      30|    30|       25|       0|       0|            0|               0|
|assign_taxo_rolling_mock      |     483|   482|      455|       0|       0|           72|             378|
|species                       |      52|    52|       34|       0|       0|            0|               0|
|species_mini                  |      12|    12|        4|       0|       0|            0|               0|
|dada2_gssp_air                |       3|     3|        3|       0|       0|            0|               0|
|assign_taxo_gssp_air          |     254|   253|      239|       0|     117|           36|             189|
|dada2_korhonen_wood           |       3|     3|        3|       0|       0|            0|               0|
|assign_taxo_korhonen_wood     |     254|   253|      239|       0|     385|           36|             189|

### Assignment inputs (units)

|project                       |input    | real_units| controls|
|:-----------------------------|:--------|----------:|--------:|
|assign_taxo                   |ITS1_asv |        195|      197|
|assign_taxo                   |ITS1_otu |        150|      175|
|assign_taxo_mini              |ITS1_asv |        195|      197|
|assign_taxo_mini              |ITS1_otu |        150|      175|
|assign_taxo_tedersoo_illumina |ITS2_asv |        199|      199|
|assign_taxo_tedersoo_illumina |full_asv |        237|      218|
|assign_taxo_tedersoo_illumina |ITS1_asv |        213|      206|
|assign_taxo_tedersoo_illumina |ITS2_otu |        176|      188|
|assign_taxo_tedersoo_illumina |full_otu |        180|      190|
|assign_taxo_tedersoo_illumina |ITS1_otu |        175|      187|
|assign_taxo_insilico          |ITS2_asv |         68|      134|
|assign_taxo_insilico          |full_asv |         71|      135|
|assign_taxo_insilico          |ITS1_asv |         71|      135|
|assign_taxo_insilico          |ITS2_otu |         68|      134|
|assign_taxo_insilico          |full_otu |         68|      134|
|assign_taxo_insilico          |ITS1_otu |         70|      135|
|assign_taxo_rolling_mock      |ITS1_asv |        109|      154|
|assign_taxo_rolling_mock      |ITS1_otu |         49|      124|
|assign_taxo_gssp_air          |ITS2_asv |       2572|     1386|
|assign_taxo_korhonen_wood     |ITS2_asv |        875|      537|

### Tests

17 test files, 95 `test_that()` blocks, 382 `expect_*()` calls (counted in the sources).

