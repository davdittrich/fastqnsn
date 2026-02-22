library(fastqnsn)
library(robustbase)
library(microbenchmark)
library(dplyr)
library(tidyr)

ns <- c(10, 10000, 1000000)
results <- data.frame(N = integer(), Estimator = character(), Pkg = character(), Time_us = numeric())

cat("\n--- Benchmarking fastqnsn vs robustbase ---\n")

for (n in ns) {
  cat(sprintf("Benchmarking N = %d...\n", n))
  set.seed(42)
  x <- rnorm(n)
  
  iters <- if (n <= 100) 100 else if (n <= 10000) 20 else 5
  
  # Sn
  res_sn_ref <- microbenchmark(val = robustbase::Sn(x), times = iters)
  res_sn_new <- microbenchmark(val = fastqnsn::sn(x), times = iters)
  
  results <- rbind(results, data.frame(N = n, Estimator = "Sn", Pkg = "robustbase", Time_us = mean(res_sn_ref$time) / 1000))
  results <- rbind(results, data.frame(N = n, Estimator = "Sn", Pkg = "fastqnsn", Time_us = mean(res_sn_new$time) / 1000))
  
  # Qn
  res_qn_ref <- microbenchmark(val = robustbase::Qn(x), times = iters)
  res_qn_new <- microbenchmark(val = fastqnsn::qn(x), times = iters)
  
  results <- rbind(results, data.frame(N = n, Estimator = "Qn", Pkg = "robustbase", Time_us = mean(res_qn_ref$time) / 1000))
  results <- rbind(results, data.frame(N = n, Estimator = "Qn", Pkg = "fastqnsn", Time_us = mean(res_qn_new$time) / 1000))
}

summary_tab <- results %>% 
  pivot_wider(names_from = Pkg, values_from = Time_us) %>%
  mutate(Speedup = robustbase / fastqnsn) %>%
  arrange(Estimator, N)

print(summary_tab)
write.csv(summary_tab, "bench_ref_summary.csv", row.names = FALSE)
