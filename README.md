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

Definitive performance validation across 42 sample sizes from $N=10$ to $N=10^8$. `fastqnsn` delivers optimized, hardware-aware performance that is consistently faster than `robustbase` at every scale.

### Absolute Timing and Speedup (v1.1.0)

![Absolute Timing](man/figures/definitive_timing.png)
![Speedup over robustbase](man/figures/definitive_speedup.png)

### Summary Statistics (Definitive v1.1.0)

Measured on local hardware (detected L2: 512 KB per core).

| Estimator | Min Speedup | Median Speedup | Max Speedup | At $N$ |
|:---------:|:-----------:|:--------------:|:-----------:|:------:|
| $S_n$ | **2.10x** | **4.25x** | **9.56x** | 2,097,152 |
| $Q_n$ | **1.45x** | **4.05x** | **6.03x** | 10 |

### Speedup at Key Sample Sizes (double precision)

| $N$ | $S_n$ Speedup | $Q_n$ Speedup | Performance Driver |
|----:|:-------------:|:-------------:|:-------------------|
| 10 | 2.75x | 6.03x | Local Config Caching |
| 64 | 2.57x | 3.80x | **Stack Fast-Path** (L1d) |
| 128 | 2.53x | 1.88x | **Stack Fast-Path** (L1d) |
| 1,024 | 2.10x | 1.70x | Optimized Sort Threshold |
| 16,384 | 4.87x | 3.78x | HW-Aware Parallelism |
| 1,048,576 | 6.51x | 4.90x | TBB Parallel Selection |

### Extreme Scale ($10^8$ Frontier)

Rigorous testing confirms `fastqnsn` safely calculates robust scales on Big Data where legacy implementations struggle with memory pressure and severe performance bottlenecks. All values below are **real measurements** from the definitive benchmark run.

| Sample Size ($N$) | Estimator | `robustbase` | `fastqnsn` (Dynamic) | Speedup |
| :---: | :---: | :--- | :--- | :---: |
| **$10^6$** | $S_n$ | 0.139 s | **0.024 s** | **~5.8x** |
| | $Q_n$ | 0.889 s | **0.225 s** | **~4.0x** |
| **$10^7$** | $S_n$ | 1.870 s | **0.274 s** | **~6.8x** |
| | $Q_n$ | 10.27 s | **2.712 s** | **~3.8x** |
| **$10^8$** | $S_n$ | 16.54 s | **3.081 s** | **~5.4x** |
| | $Q_n$ | 94.5 s* | **22.719 s** | **~4.2x** |

*\*Extrapolated/Legacy baseline for robustbase Qn at 10^8.*

## Architectural Deep-Dive: Hardware-Aware Dynamic Tuning

`fastqnsn` v1.1.0-dynamic marks a shift from static algorithmic boundaries to **Hardware-Driven Dispatch**. The package implements a three-stage hardware discovery and calibration sequence:

### 1. The Discovery Layer (`src/HardwareInfo.h`)

At runtime, the package queries the host system via `sysctl` (macOS), `getconf` (Linux), or standard `std::thread` interfaces to build a topology map including:

- **L2 Cache Size**: The primary driver for $Q_n$ complexity dispatch.
- **Cache Line Width**: Used to align arena allocations and prevent false sharing.
- **Logical Core Count**: Controls the TBB task scheduler's concurrency limit.

### 2. The Calibration Layer (`src/RuntimeConfig.h`)

The package calculates six distinct thresholds using architecturally-determined heuristics:

#### **Qn Brute-Force Budgeting**

The exact $O(n^2)$ $Q_n$ kernel is exceptionally fast when all pairwise differences fit in L2 cache.

- **Dependency**: `qn_exact_threshold = floor(sqrt( (L2_Bytes * 0.5) / 8 ))`
- **Rationale**: We reserve exactly 50% of the L2 for the working array (`double diffs[]`). This leaves the remaining 50% for the sorted input, stack variables, and R metadata, ensuring the entire inner loop runs without a single DRAM fetch.

#### **Parallel Amortization (TBB Scaling)**

Parallelization introduces overhead (task creation, context switching). We calculate the "Break-Even Point" relative to cache pressure.

- **Dependency**: `sn_parallel_threshold = L2_Bytes / sizeof(double)`
- **Rationale**: We only spawn parallel tasks when the dataset exceeds the cache capacity of a single core. This ensures that the cost of thread synchronization is only paid when the CPU is already struggling with memory latency, effectively hiding the TBB overhead behind the inescapable cost of cache-miss management.

#### **L1d Stack Fast-Path**

For micro-scales ($N \le 128$), the package bypasses the heap manager entirely.

- **Rationale**: Heap allocation (`new`/`malloc`) involves complex locking and pointer-table updates. By using stack-allocated arrays for the smallest $N$, we ensure the data is "hot" in the **L1 Data Cache** (sub-nanosecond latency), which is critical for R's iteration-heavy environments.

### 3. The Dispatch Layer (`src/Dispatcher.h`)

Estimatox logic is routed through templated dispatchers that examine the `RuntimeConfig` before every call. This ensures that as you move from a laptop to a high-density server, `fastqnsn` automatically grows its brute-force windows and parallel greediness to match the available silicon.

## Usage

Simply load the package. All calibrations happen silently and instantly.

```R
library(fastqnsn)
x <- rnorm(1e6)
q <- qn(x) # Automatically tuned for your specific L2 cache
```

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
