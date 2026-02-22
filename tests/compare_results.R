library(dplyr)
library(tidyr)

new_res <- read.csv("bench_new_zig.csv")
old_res <- read.csv("bench_old_local.csv")

all_res <- rbind(new_res, old_res)

summary_tab <- all_res %>% 
  pivot_wider(names_from = Version, values_from = Time_us) %>%
  mutate(Speedup = old_local / new_zig) %>%
  arrange(Estimator, N)

cat("| Estimator | Sample Size (N) | Old Version (µs) | New Zig Version (µs) | Speedup |\n")
cat("| :--- | :--- | :--- | :--- | :--- |\n")
for (i in seq_len(nrow(summary_tab))) {
  cat(sprintf("| %s | %d | %.2f | %.2f | **%.2fx** |\n", 
              summary_tab$Estimator[i],
              summary_tab$N[i],
              summary_tab$old_local[i],
              summary_tab$new_zig[i],
              summary_tab$Speedup[i]))
}
