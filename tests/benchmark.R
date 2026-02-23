library(fastqnsn)
library(robustbase)
library(microbenchmark)
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

# Configuration
set.seed(42)
iterations <- 30

sample_sizes <- unique(c(seq(10, 100, by=10),
                         seq(200, 1000, by=100),
                         seq(2000, 10000, by=1000),
                         seq(20000, 100000, by=10000),
                         seq(200000, 1000000, by=100000)))

cat(sprintf("Total sample sizes: %d\n", length(sample_sizes)))

results_list <- list()

for (n in sample_sizes) {
  cat(sprintf("[%s] Benchmarking n = %d...\n", Sys.time(), n))
  x <- rnorm(n)

  bm <- microbenchmark(
    sn_robustbase = robustbase::Sn(x),
    sn_fastqnsn = fastqnsn::sn(x),
    qn_robustbase = robustbase::Qn(x),
    qn_fastqnsn = fastqnsn::qn(x),
    times = iterations
  )

  sum_bm <- summary(bm, unit = "ms")
  sum_bm$n <- n
  results_list[[as.character(n)]] <- sum_bm
}

all_results <- bind_rows(results_list)

# Generate the plot
plot_data <- all_results %>%
  mutate(
    Estimator = ifelse(grepl("^sn_", expr), "Sn", "Qn"),
    Package = ifelse(grepl("_fastqnsn$", expr), "fastqnsn", "robustbase"),
    Group = paste(Estimator, Package, sep = " - ")
  )

p <- ggplot(plot_data, aes(x = n, y = median, color = Group, linetype = Estimator)) +
  geom_line(linewidth = 1) +
  geom_point(size = 1.5) +
  scale_x_log10(labels = scales::label_number(scale_cut = cut_short_scale())) +
  scale_y_log10() +
  labs(
    title = "Performance Comparison: fastqnsn vs robustbase",
    subtitle = paste("Median execution time over", iterations, "iterations"),
    x = "Sample Size (n)",
    y = "Median Time (ms)",
    color = "Implementation",
    linetype = "Estimator"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  ) +
  scale_color_brewer(palette = "Set1")

dir.create("man/figures", showWarnings = FALSE, recursive = TRUE)
ggsave("man/figures/benchmark.png", p, width = 10, height = 7)
cat("\nBenchmark plot saved to man/figures/benchmark.png\n")

# Prepare summary
summary_pivot <- plot_data %>%
  group_by(n, Estimator, Package) %>%
  summarize(median = median(median), .groups = "drop") %>%
  pivot_wider(names_from = Package, values_from = median) %>%
  mutate(speedup = robustbase / fastqnsn)

write.csv(summary_pivot, "tests/benchmark_results.csv", row.names = FALSE)
