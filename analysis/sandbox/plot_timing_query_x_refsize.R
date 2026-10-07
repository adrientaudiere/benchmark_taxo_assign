# Figure of Archives/sandbox_pre_od34/timing_query_x_refsize.R: elapsed time of one
# computation against the number of queries, one line per reference size, one
# panel per method, with the points measured by the pipeline runs (same
# database, 206 to 1412 units, other computations running at the same time).
#   Rscript analysis/sandbox/plot_timing_query_x_refsize.R

library(ggplot2)

grid <- read.csv("data/data_final/timing/timing_query_x_refsize.csv")
pipe <- read.csv(
  "data/data_final/timing/timing_query_x_refsize_pipeline_points.csv"
)
grid$reference <- factor(
  paste0(format(grid$n_ref, big.mark = " "), " records"),
  levels = paste0(
    format(sort(unique(grid$n_ref), decreasing = TRUE), big.mark = " "),
    " records"
  )
)
pipe$dataset <- sub("^pipeline run \\((.*)\\), concurrent$", "\\1", pipe$source)

fig <- ggplot(grid, aes(n_query, elapsed_s)) +
  geom_line(aes(colour = reference)) +
  geom_point(aes(colour = reference), size = 1.5) +
  geom_point(
    data = pipe,
    aes(shape = dataset),
    size = 2.5,
    colour = "black"
  ) +
  facet_wrap(~method, nrow = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "number of queries (log)",
    y = "elapsed time of one computation, s (log)",
    colour = "Unite_all_20250219_Fungi",
    shape = "pipeline run,\nwhole database"
  ) +
  theme_bw()

dir.create("figures", showWarnings = FALSE)
ggsave(
  "figures/timing_query_x_refsize.png",
  fig,
  width = 12,
  height = 4.5,
  dpi = 150
)
