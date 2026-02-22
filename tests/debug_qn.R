library(fastqnsn)
library(robustbase)

set.seed(42)
n <- 150
x <- rnorm(n)

# Reference
h <- floor(n/2) + 1
k <- h * (h - 1) / 2
constant_qn <- 2.21914446598508
# Factors for n=150
predict_qn_factor <- function(i) if (i %% 2 == 1) 1 - 1.594/i + 3.22/i^2 else 1 - 3.672/i + 11.087/i^2
factor <- predict_qn_factor(n)
ref_val_qn <- robustbase::Qn(x, constant = constant_qn * factor)

# New
new_val_qn <- fastqnsn::qn(x)

cat("n:", n, "\n")
cat("k:", k, "\n")
cat("Factor:", factor, "\n")
cat("Constant:", constant_qn, "\n")
cat("Combined Constant:", constant_qn * factor, "\n")
cat("Reference Qn:", ref_val_qn, "\n")
cat("New Qn:", new_val_qn, "\n")
cat("Difference:", new_val_qn - ref_val_qn, "\n")

# Get raw differences for comparison
all_diffs <- sort(dist(x))
raw_ref <- all_diffs[k]
raw_new <- new_val_qn / (constant_qn * factor)

cat("Reference Raw:", raw_ref, "\n")
cat("New Raw:", raw_new, "\n")
cat("Raw Difference:", raw_new - raw_ref, "\n")
