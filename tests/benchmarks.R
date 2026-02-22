args <- commandArgs(trailingOnly = TRUE)
version_name <- if (length(args) > 0) args[1] else "unknown"

library(fastqnsn)
library(microbenchmark)

ns <- c(10, 100, 1000, 10000, 100000)
results <- data.frame(Version = character(), N = integer(), Estimator = character(), Time_us = numeric())

cat(sprintf("\n--- Starting Benchmarks for %s ---\n", version_name))

for (n in ns) {
  cat(sprintf("Benchmarking N = %d...\n", n))
  set.seed(42)
  x <- rnorm(n)
  
  # Number of iterations scales down for larger n to save time
  iters <- if (n <= 1000) 50 else if (n <= 10000) 20 else 5
  
  # Sn
  res_sn <- microbenchmark(val = fastqnsn::sn(x), times = iters)
  mean_sn_us <- mean(res_sn$time) / 1000
  results <- rbind(results, data.frame(
    Version = version_name, N = n, Estimator = "Sn", Time_us = mean_sn_us
  ))
  
  # Qn
  res_qn <- microbenchmark(val = fastqnsn::qn(x), times = iters)
  mean_qn_us <- mean(res_qn$time) / 1000
  results <- rbind(results, data.frame(
    Version = version_name, N = n, Estimator = "Qn", Time_us = mean_qn_us
  ))
}

out_file <- sprintf("bench_%s.csv", version_name)
write.csv(results, out_file, row.names = FALSE)
cat(sprintf("Wrote results to %s\n", out_file))
