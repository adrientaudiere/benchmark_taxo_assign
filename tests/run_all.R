# Runs every tests/test_*.R file. From the project root:
#   Rscript tests/run_all.R
# Each file is self-contained (loads what it needs), so they can also be run
# one by one with Rscript tests/<file>.R.

library("testthat")
library("here")

here::i_am("tests/run_all.R")
res <- testthat::test_dir(here("tests"), reporter = "summary", stop_on_failure = FALSE)
df <- as.data.frame(res)
cat(sprintf(
  "\n%d test files, %d expectations, %d failed, %d errors\n",
  length(unique(df$file)), sum(df$nb), sum(df$failed), sum(df$error)
))
if (sum(df$failed) + sum(df$error) > 0) {
  quit(status = 1)
}
