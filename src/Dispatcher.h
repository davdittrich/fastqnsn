#ifndef FASTQNSN_DISPATCHER_H
#define FASTQNSN_DISPATCHER_H

#include "QnKernels.h"
#include "RuntimeConfig.h"

namespace fastqnsn_dynamic {

class Dispatcher {
public:
  template <typename T>
  static void qn_brute_force(const T *sorted_x, size_t n, double *diffs) {
    auto &config = RuntimeConfig::get();

    // SIMD kernels currently only optimized for double
    if constexpr (std::is_same_v<T, double>) {
      if (config.hw.simd_level == SIMDLevel::AVX2 && n >= 64 && n <= 1024) {
        // AVX2 shows benefit in this range in our A/B tests
        qn_brute_force_avx2(sorted_x, n, diffs);
        return;
      }
      if (config.hw.simd_level == SIMDLevel::Neon && n >= 64 && n <= 1024) {
        qn_brute_force_neon(sorted_x, n, diffs);
        return;
      }
    }

    // Default to scalar for other types or if no SIMD available
    qn_brute_force_scalar(sorted_x, n, diffs);
  }

  // Future expansion for Sn, Parallel, etc.
};

} // namespace fastqnsn_dynamic

#endif // FASTQNSN_DISPATCHER_H
