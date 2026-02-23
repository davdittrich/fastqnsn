library(fastqnsn)
library(robustbase)
library(microbenchmark)
library(ggplot2)
library(tidyr)
library(dplyr)

sample_sizes <- c(
  seq(10, 100, by = 10),
  seq(200, 1000, by = 100),
  seq(2000, 10000, by = 1000),
  seq(20000, 100000, by = 10000),
  seq(200000, 1000000, by = 100000)
)

results_list <- list()

for (n in sample_sizes) {
  cat(sprintf("Benchmarking n = %d...\n", n))
  set.seed(42)
  x <- rnorm(n)

  iter <- if (n >= 200000) 5 else 20

  bm <- microbenchmark(
    fastqnsn_sn = fastqnsn::sn(x),
    fastqnsn_qn = fastqnsn::qn(x),
    robustbase_sn = robustbase::Sn(x),
    robustbase_qn = robustbase::Qn(x),
    times = iter,
    unit = "us"
  )

  bm_summary <- summary(bm, unit = "us") %>%
    mutate(n = n)

  results_list[[as.character(n)]] <- bm_summary
}

results <- do.call(rbind, results_list)

# Save raw results
saveRDS(results, "benchmark_results.rds")

# Prepare data for plotting
plot_data <- results %>%
  select(expr, median, n) %>%
  mutate(Implementation = ifelse(grepl("fastqnsn", expr), "fastqnsn", "robustbase"),
         Estimator = ifelse(grepl("_sn", expr), "Sn", "Qn"))

# Create plot
dir.create("man/figures", recursive = TRUE, showWarnings = FALSE)

p <- ggplot(plot_data, aes(x = n, y = median, color = Implementation, linetype = Estimator)) +
  geom_line(linewidth = 1) +
  geom_point() +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  labs(title = "Execution Time Comparison: fastqnsn vs robustbase",
       subtitle = "Median of 30 iterations per sample size",
       x = "Sample Size (n)",
       y = "Median Execution Time (microseconds)",
       color = "Implementation",
       linetype = "Estimator") +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("man/figures/benchmark.png", p, width = 10, height = 6)

# Generate a clean summary table
summary_table <- plot_data %>%
  pivot_wider(id_cols = n, names_from = c(Implementation, Estimator), values_from = median) %>%
  mutate(Speedup_Sn = robustbase_Sn / fastqnsn_Sn,
         Speedup_Qn = robustbase_Qn / fastqnsn_Qn)

write.csv(summary_table, "benchmark_summary.csv", row.names = FALSE)

cat("Benchmark complete. Results saved to benchmark_results.rds and man/figures/benchmark.png\n")
