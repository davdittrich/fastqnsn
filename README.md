# fastqnsn

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18727053.svg)](https://doi.org/10.5281/zenodo.18727053)

`fastqnsn` is a high-performance R package for computing the **Rousseeuw-Croux $Q_n$ and $S_n$** robust scale estimators. It is designed to provide peak performance across all sample sizes while maintaining absolute bit-identical correctness compared to `robustbase`.

## Key Features
- **Hybrid Architecture:**
  - **Serial Path ($n \le 2000$):** High-speed deterministic Zig kernels with zero threading overhead for small and medium data.
  - **Parallel Path ($n > 2000$):** Multi-threaded $O(n \log n)$ implementation using **RcppParallel (Intel TBB)**.
- **Deterministic Kernels:**
  - **$S_n$:** Implements the $O(n \log n)$ Shamos (1976) overall-median algorithm in Zig.
  - **$Q_n$:** Optimized Johnson-Mizoguchi (1978) selector with branchless counting.
- **Superior Accuracy:** 
  - Implements corrected $D_\infty = 2.21914446598508$ (fixing the legacy typo $2.2219$).
  - Uses modern finite-sample bias corrections from **Akinshin (2022)**.
- **Memory Optimized:** Allocation-free $S_n$ workers and efficient vector reuse in $Q_n$ refinement loops.
- **Robustness:** Built-in `std::isfinite` checks and 64-bit rank calculations for massive datasets.

## Installation
```R
# Requires Rcpp, RcppParallel, and a C++ compiler / Zig
# install.packages("remotes")
remotes::install_github("davdittrich/fastqnsn")
```

## Usage
```R
library(fastqnsn)
x <- rnorm(10000)

scale_sn <- sn(x)
scale_qn <- qn(x)
```

## Benchmarks
Results from comparison with `robustbase` (median execution time on $N=1,000,000$):

| Estimator | `robustbase` | `fastqnsn` | Speedup |
| :--- | :--- | :--- | :--- |
| **$S_n$** | 134.4 ms | 30.7 ms | **4.37x** |
| **$Q_n$** | 895.8 ms | 523.4 ms | **1.71x** |

Small-sample results ($n=10$):

| Estimator | `robustbase` | `fastqnsn` | Speedup |
| :--- | :--- | :--- | :--- |
| **$S_n$** | 4.6 µs | 2.0 µs | **2.3x** |
| **$Q_n$** | 10.0 µs | 2.3 µs | **4.3x** |

*Note: `fastqnsn` provides bit-identical results to `robustbase` when matching consistency constants are used.*

## Authors
**Dennis Alexis Valin Dittrich** (ORCID: 0000-0002-4438-8276)  

## References
- Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the Median Absolute Deviation. *JASA*.
- Akinshin, A. (2022). Finite-sample Rousseeuw-Croux scale estimators. *arXiv:2209.12268*.
- Johnson, D. B., & Mizoguchi, T. (1978). Selecting the Kth element in X + Y. *SIAM J. Comput.*
