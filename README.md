# fastqnsn

`fastqnsn` is a high-performance R package for computing the **Rousseeuw-Croux $Q_n$ and $S_n$** robust scale estimators. It is designed to provide peak performance across all sample sizes while maintaining superior numerical accuracy compared to legacy implementations.

## Key Features
- **Hybrid Architecture:**
  - **Micro-Scale ($n < 100$):** High-speed scalar $O(n^2)$ engine with zero threading overhead.
  - **Macro-Scale ($n \ge 100$):** Multi-threaded $O(n \log n)$ implementation using **RcppParallel (Intel TBB)**.
- **Superior Accuracy:** 
  - Implements corrected $D_\infty = 2.21914446598508$ (fixing the legacy typo $2.2219$).
  - Uses modern finite-sample bias corrections from **Akinshin (2022)**.
- **Memory Optimized:** Allocation-free $S_n$ workers and efficient vector reuse in $Q_n$ refinement loops.
- **Robustness:** Built-in `std::isfinite` checks and 64-bit rank calculations to handle large datasets safely.

## Installation
```R
# Requires Rcpp, RcppParallel, and a C++ compiler
remotes::install_local("fastqnsn")
```

## Usage
```R
library(fastqnsn)
x <- rnorm(10000)

scale_sn <- sn(x)
scale_qn <- qn(x)
```

## Benchmarks
Results from comparison with `robustbase` (median execution time):

| Sample Size ($n$) | Estimator | `robustbase` | `fastqnsn` | Speedup |
| :--- | :--- | :--- | :--- | :--- |
| **10** | $S_n$ | 4.10 µs | 1.97 µs | **~2.1x** |
| **10** | $Q_n$ | 8.83 µs | 1.76 µs | **~5.0x** |
| **10,000** | $S_n$ | 916 µs | 580 µs | **~37%** |
| **10,000** | $Q_n$ | 5,798 µs | 4,250 µs | **~26%** |

*Note: `fastqnsn` provides bit-identical results to `robustbase` when matching consistency constants are provided.*

## Authors
**Dennis Alexis Valin Dittrich** (ORCID: 0000-0002-4438-8276)  
Email: [davd@economicscience.net](mailto:davd@economicscience.net)

## References
- Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the Median Absolute Deviation. *JASA*.
- Akinshin, A. (2022). Finite-sample Rousseeuw-Croux scale estimators. *arXiv:2209.12268*.
- Johnson, D. B., & Mizoguchi, T. (1978). Selecting the Kth element in X + Y. *SIAM J. Comput.*
