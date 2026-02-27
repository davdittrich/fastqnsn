# fastqnsn

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18727053.svg)](https://doi.org/10.5281/zenodo.18727053)

`fastqnsn` is a high-performance R package for computing the **Rousseeuw-Croux $Q_n$ and $S_n$** robust scale estimators. It delivers consistent speedups over `robustbase` across all sample sizes from $N=10$ to $N=10^8$, with cache-aware algorithm dispatch that self-tunes to the target CPU architecture at install time.

## Key Features

- **Cache-Aware Hybrid Architecture:** Six threshold parameters are derived from the CPU's L2 cache size (detected at install time via `sysctl`/`getconf`), controlling algorithm dispatch across three regimes:
  - **Micro-Scale ($N \le 2048$ for $Q_n$):** Ultra-fast $O(n^2)$ exact brute-force kernel. Working set sized to fit within L2 cache.
  - **Mid-Scale (serial $O(n \log n)$):** Johnson-Mizoguchi iterative algorithm for $Q_n$; sweep-based algorithm for $S_n$. Parallelization thresholds ($S_n$: 12288, $Q_n$: 8192) are tuned to avoid premature thread spawning overhead.
  - **Macro-Scale:** Parallelized counting and refinement via **RcppParallel (Intel TBB)**.
- **Floyd-Rivest Selection:** Replaces `std::nth_element` throughout, achieving ~30% fewer comparisons.
- **Arena Memory Allocation:** Single contiguous allocation for all working arrays in both $Q_n$ and $S_n$.
- **Three-Tier Sorting:** `std::sort` for $N \le 256$, Boost Spreadsort for medium $N$, TBB `parallel_sort` for large $N$ (float threshold: 6144, integer threshold: 8192).
- **Superior Accuracy:**
  - Corrected $D_\infty = 2.21914446598508$ (fixing the legacy approximation $2.2219$).
  - Modern finite-sample bias corrections from **Akinshin (2022)**.
  - `(float)` truncation matching robustbase precision semantics.

## Installation

```R
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

Validated across 61 sample sizes from $N=10$ to $N=1{,}000{,}000$, both $S_n$ and $Q_n$ estimators, on double and integer data. `fastqnsn` is faster than `robustbase` at **every** sample size tested (30 iterations per measurement, `microbenchmark`).

### Speedup over robustbase

![Speedup over robustbase](man/figures/validation_speedup.png)

### Absolute Timing

![Absolute Timing](man/figures/validation_timing.png)

### Summary Statistics (v1.1.0 Dynamic)

| Estimator | Min Speedup | Median Speedup | Max Speedup | At $N$ |
|:---------:|:-----------:|:--------------:|:-----------:|:------:|
| $S_n$ | **2.21x** | **4.33x** | **9.89x** | 2,097,152 |
| $Q_n$ | **1.74x** | **4.11x** | **5.84x** | 8 |

### Speedup at Key Sample Sizes (double precision)

| $N$ | $S_n$ Speedup | $Q_n$ Speedup | New in v1.1.0 |
|----:|:-------------:|:-------------:|:--------------|
| 10 | 2.69x | 5.78x | Local Config Caching |
| 64 | 2.63x | 4.31x | **Stack Fast-Path** |
| 128 | 2.45x | 2.21x | **Stack Fast-Path** |
| 1,024 | 2.15x | 1.89x | Optimized Sort Threshold |
| 16,384 | 4.43x | 3.48x | HW-Aware Parallelism |
| 1,048,576 | 6.44x | 4.89x | TBB Parallel Selection |

### Extreme Scale ($10^8$ Frontier)

Rigorous testing up to $N=10^8$ confirms `fastqnsn` safely calculates robust scales on Big Data where legacy implementations struggle with memory pressure and severe performance bottlenecks.

| Sample Size ($N$) | Estimator | `robustbase` | `fastqnsn` (Dynamic) | Speedup |
| :---: | :---: | :--- | :--- | :---: |
| **$10^6$** | $S_n$ | 0.067 s | **0.012 s** | **~5.5x** |
| | $Q_n$ | 0.430 s | **0.095 s** | **~4.5x** |
| **$10^7$** | $S_n$ | 1.310 s | **0.216 s** | **~6.1x** |
| | $Q_n$ | 22.8 s* | **6.51 s** | **~3.5x** |
| **$10^8$** | $S_n$ | 16.1 s* | **2.65 s** | **~6.1x** |
| | $Q_n$ | 94.5 s* | **30.17 s** | **~3.1x** |

*\*Extrapolated or from legacy benchmarks where robustbase limits were exceeded.*

## Runtime Hardware Tuning

`fastqnsn` v1.1.0-dynamic introduces true **Runtime Hardware Discovery**. Unlike static versions, it detects your CPU's L1/L2 cache topology and core count at execution time to optimize every calculation.

| Feature | Dynamic Logic | Benefit |
|:----------|:---|:---------|
| **Stack Fast-Path** | Activated for $N \le 128$ | Zero heap-allocation latency |
| **Qn Brute-Force** | Sized to 50% of available L2 | Optimal SIMD throughput |
| **Sn Parallelism** | Threshold = $L2 / sizeof(double)$ | Amortizes thread-spawning cost |
| **Lookup Caching** | Singleton bypass in hot loops | 21% reduction in branch overhead |

No manual configuration is required. The package self-calibrates to provide the fastest possible estimator on any machine, from laptops to high-core servers.

### Cross-Platform Cache Detection

| Platform | Detection Method | Fallback |
|:---------|:-----------------|:---------|
| macOS (Apple Silicon) | `sysctl -n hw.perflevel1.l2cachesize` (E-core L2) | 4 MB |
| macOS (Intel) | `sysctl -n hw.l2cachesize` | 4 MB |
| Linux | `getconf LEVEL2_CACHE_SIZE` | 4 MB |
| Windows | Static default | 4 MB |

*Note: `fastqnsn` uses updated consistency constants and finite-sample bias corrections from Akinshin (2022).*

## Authors

**Dennis Alexis Valin Dittrich** (ORCID: 0000-0002-4438-8276)

## References

- Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the Median Absolute Deviation. *JASA*.
- Akinshin, A. (2022). Finite-sample Rousseeuw-Croux scale estimators. *arXiv:2209.12268*.
- Johnson, D. B., & Mizoguchi, T. (1978). Selecting the Kth element in X + Y. *SIAM J. Comput.*
- Floyd, R. W., & Rivest, R. L. (1975). Expected time bounds for selection. *CACM*.
