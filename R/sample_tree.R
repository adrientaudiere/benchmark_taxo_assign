# Memory of a process and of all its descendants, sampled from outside it
# (code review §1.3). Run as a script by with_autometric()
# (R/autometric_helpers.R):
#
#   Rscript R/sample_tree.R <pid> <csv> <seconds>
#
# Appends "epoch,rss_mb,n_procs" every <seconds> for <pid> and its
# descendants, the sampler itself excluded (it is one of them), until <pid>
# exits or the sampler is killed. autometric samples only the R process,
# which sits idle while vsearch or BLAST+ run as its children, so their memory
# was never measured. RSS is summed over the processes, so pages they share
# count once per process: an upper bound on the tree's footprint.

args <- commandArgs(trailingOnly = TRUE)
root <- ps::ps_handle(as.integer(args[[1]]))
path <- args[[2]]
seconds <- as.numeric(args[[3]])
self <- Sys.getpid()

cat("epoch,rss_mb,n_procs\n", file = path)
while (ps::ps_is_running(root)) {
  kids <- tryCatch(
    ps::ps_children(root, recursive = TRUE),
    error = function(e) list()
  )
  kids <- Filter(function(h) ps::ps_pid(h) != self, kids)
  procs <- c(list(root), kids)
  rss <- vapply(
    procs,
    function(h) {
      tryCatch(ps::ps_memory_info(h)[["rss"]], error = function(e) 0)
    },
    numeric(1)
  )
  cat(
    sprintf(
      "%.3f,%.1f,%d\n",
      as.numeric(Sys.time()),
      sum(rss) / 1e6,
      length(procs)
    ),
    file = path,
    append = TRUE
  )
  Sys.sleep(seconds)
}
