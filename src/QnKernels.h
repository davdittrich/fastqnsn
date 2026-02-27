#ifndef FASTQNSN_QN_KERNELS_H
#define FASTQNSN_QN_KERNELS_H

#include <algorithm>
#include <cmath>
#include <immintrin.h>
#include <vector>

#ifdef __ARM_NEON
#include <arm_neon.h>
#endif

namespace fastqnsn_dynamic {

// Scalar baseline
template <typename T>
void qn_brute_force_scalar(const T *sorted_x, size_t n, double *diffs) {
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    for (size_t j = 0; j < i; ++j) {
      diffs[k++] = (double)sorted_x[i] - (double)sorted_x[j];
    }
  }
}

// AVX2 optimized version for double
#if defined(__x86_64__) || defined(_M_X64)
__attribute__((target("avx2"))) inline void
qn_brute_force_avx2(const double *sorted_x, size_t n, double *diffs) {
#else
inline void qn_brute_force_avx2(const double *sorted_x, size_t n,
                                double *diffs) {
#endif
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    double xi = sorted_x[i];
    __m256d v_xi = _mm256_set1_pd(xi);

    size_t j = 0;
    // Process 4 elements at a time
    for (; j + 3 < i; j += 4) {
      __m256d v_xj = _mm256_loadu_pd(&sorted_x[j]);
      __m256d v_diff = _mm256_sub_pd(v_xi, v_xj);
      // Since x is sorted, xi >= xj, so v_diff is always >= 0. No need for abs.
      _mm256_storeu_pd(&diffs[k], v_diff);
      k += 4;
    }

    // Remainder
    for (; j < i; ++j) {
      diffs[k++] = xi - sorted_x[j];
    }
  }
}

// Neon optimized version for double
inline void qn_brute_force_neon(const double *sorted_x, size_t n,
                                double *diffs) {
#ifdef __ARM_NEON
  size_t k = 0;
  for (size_t i = 1; i < n; ++i) {
    double xi = sorted_x[i];
    float64x2_t v_xi = vdupq_n_f64(xi);

    size_t j = 0;
    // Process 2 elements at a time (Neon double vectors are 128-bit)
    for (; j + 1 < i; j += 2) {
      float64x2_t v_xj = vld1q_f64(&sorted_x[j]);
      float64x2_t v_diff = vsubq_f64(v_xi, v_xj);
      vst1q_f64(&diffs[k], v_diff);
      k += 2;
    }

    // Remainder
    for (; j < i; ++j) {
      diffs[k++] = xi - sorted_x[j];
    }
  }
#endif
}

} // namespace fastqnsn_dynamic

#endif // FASTQNSN_QN_KERNELS_H
