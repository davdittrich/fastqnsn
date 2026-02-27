#ifndef FASTQNSN_THRESHOLDS_H
#define FASTQNSN_THRESHOLDS_H

#include <cstddef>

// ============================================================================
// Cache-aware threshold constants for fastqnsn
//
// These thresholds control algorithm dispatch decisions:
//   - Qn: brute-force O(n^2) vs Johnson-Mizoguchi iterative
//   - Sn: stack-allocated vs heap-allocated working arrays
//   - Sn/Qn: serial vs parallel (RcppParallel/TBB)
//   - Sort: std::sort vs Boost spreadsort vs TBB parallel_sort
//
// Values are architecture-aligned (powers of 2 or 3*2^(k-1)) and derived
// from the target CPU's L2 cache size, detected at package install time
// via Makevars. The compile-time constant FASTQNSN_L2_CACHE_BYTES is set
// by Makevars; if unavailable, a conservative 4MB default is used.
//
// Empirical validation: ~/fastqnsn_benchmarks/ (Phase 2A-2D, Phase 3)
// ============================================================================

// --- Compile-time cache parameters (set by Makevars, with safe defaults) ---

#ifndef FASTQNSN_L2_CACHE_BYTES
#define FASTQNSN_L2_CACHE_BYTES 4194304 // 4 MB conservative default
#endif

#ifndef FASTQNSN_CACHE_LINE_BYTES
#ifdef __aarch64__
#define FASTQNSN_CACHE_LINE_BYTES 128 // ARM / Apple Silicon
#else
#define FASTQNSN_CACHE_LINE_BYTES 64 // x86-64
#endif
#endif

namespace fastqnsn {

// --- Helper: constexpr integer square root (Newton's method) ---
constexpr size_t isqrt(size_t n) {
  if (n == 0)
    return 0;
  size_t x = n;
  for (int i = 0; i < 32; ++i) {
    if (x == 0)
      return 0;
    x = (x + n / x) / 2;
  }
  while (x * x > n)
    --x;
  return x;
}

// --- Helper: round down to nearest architecture-aligned value ---
// Produces values from the sequence: ..., 256, 384, 512, 768, 1024, 1536, ...
// These are 2^k and 3*2^(k-1), which align with cache boundaries.
constexpr size_t round_to_arch_aligned(size_t n) {
  if (n == 0)
    return 0;
  size_t p = 1;
  while (p * 2 <= n)
    p *= 2;
  size_t mid = 3 * (p / 2); // = 3 * 2^(k-1) where p = 2^k
  if (mid <= n && mid > p)
    return mid;
  return p;
}

// ============================================================================
// Qn exact threshold
// ============================================================================
// The brute-force Qn path allocates n*(n-1)/2 doubles for pairwise differences.
// It wins when this working set fits in L2 cache. Use 50% of the smallest L2
// as budget (leaving room for sorted_x array and other working data).
//
// Formula: n_max = floor(sqrt(0.5 * L2 / sizeof(double)))
//                = floor(sqrt(L2 / 8))
//
// Empirical validation (Phase 2A):
//   On Apple M3 Pro (E-core L2=4MB, P-core L2=16MB):
//   - JM beats brute-force at n >= 128 (native optimized)
//
// We use 128 as an empirical default for modern x86_64 to stay in L2.
constexpr size_t QN_EXACT_THRESHOLD = 128;

// ============================================================================
// Sn stack threshold
// ============================================================================
// Stack arrays: 2 * threshold * sizeof(T). At 2048 doubles: 32KB (fits L1d).
// Empirical: no measurable impact (Phase 2B). Keep at 2048.
constexpr size_t SN_STACK_THRESHOLD = 2048;

// ============================================================================
// Parallel thresholds (serial -> RcppParallel/TBB)
// ============================================================================
// Thread spawning overhead on TBB: ~10-50 microseconds.
// Below the threshold, serial execution avoids this overhead.
// Above the threshold, multi-core parallelism amortizes the cost.
//
// Empirical (Phase 2C, Apple M3 Pro 10-core):
//   Sn: threshold=12288 has best geom mean (0.013311 us/elem)
//   Qn: threshold=8192 avoids n=6144-8192 overhead while still
//       parallelizing at n>8192 (including n=12288+).
constexpr size_t SN_PARALLEL_THRESHOLD = 12288;
constexpr size_t QN_PARALLEL_THRESHOLD = 8192;

// ============================================================================
// Sort thresholds
// ============================================================================
// Three-tier sort dispatch: std::sort -> Boost spreadsort -> TBB parallel_sort
//
// Empirical (Phase 2D, measured end-to-end via Sn):
//   Float: threshold=6144 best (total=2077.7, geom=0.011469)
//          At n=8192: TBB = 152.0 us vs spreadsort = 104.8 us (45% overhead)
//   Integer: threshold=8192 best (total=1535.5, geom=0.008675)
//            Boost integer_sort's radix approach is efficient up to ~8K.
constexpr size_t SORT_BOOST_THRESHOLD = 256; // std::sort -> Boost spreadsort
constexpr size_t SORT_TBB_FLOAT_THRESHOLD = 6144; // spreadsort -> TBB (float)
constexpr size_t SORT_TBB_INT_THRESHOLD = 8192;   // spreadsort -> TBB (integer)

} // namespace fastqnsn

#endif // FASTQNSN_THRESHOLDS_H
