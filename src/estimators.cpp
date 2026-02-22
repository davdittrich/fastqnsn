#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include "kernels.h"
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <tbb/parallel_sort.h>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// --- UTILITIES ---

inline double lowmedian_ptr(double *arr, size_t n) {
  if (n == 0)
    return 0.0;
  size_t m = (size_t)(std::floor(((double)n + 1.0) / 2.0) - 1.0);
  std::nth_element(arr, arr + m, arr + n);
  return arr[m];
}

// --- SN ESTIMATOR WORKER (LARGE N) ---

struct SnWorker : public Worker {
  const double *sorted_x;
  size_t n;
  double *results;

  SnWorker(const double *sorted_x, size_t n, double *results)
      : sorted_x(sorted_x), n(n), results(results) {}

  void operator()(size_t begin, size_t end) {
    for (size_t i = begin; i < end; ++i) {
      // S_n selection logic
      // The kth element in { |x_i - x_j| } where j != i
      // This is handled by zig_sn_select_k
      size_t target_rank = (n - 1) / 2; // 0-indexed target rank for median
      results[i] = zig_sn_select_k(sorted_x, n, i, target_rank + 1);
    }
  }
};

// [[Rcpp::export]]
double C_sn_fast(NumericVector x) {
  size_t n = x.size();
  if (n < 2)
    return NA_REAL;

  // Direct access to R vector memory
  const double *x_ptr = x.begin();

  if (n < 50) {
    std::vector<double> inner_medians(n);
    zig_sn_naive(x_ptr, n, inner_medians.data());
    return lowmedian_ptr(inner_medians.data(), n);
  }

  std::vector<double> sorted_x(n);
  std::copy(x_ptr, x_ptr + n, sorted_x.begin());

  if (n > 5000)
    tbb::parallel_sort(sorted_x.begin(), sorted_x.end());
  else
    std::sort(sorted_x.begin(), sorted_x.end());

  std::vector<double> inner_medians(n);
  SnWorker worker(sorted_x.data(), n, inner_medians.data());

  if (n > 2000)
    parallelFor(0, n, worker);
  else
    worker(0, n);

  return lowmedian_ptr(inner_medians.data(), n);
}

// --- QN ESTIMATOR HELPERS ---

struct QnCounter : public Worker {
  const double *sorted_x_ptr;
  double trial;
  size_t n;
  long long sumP;

  QnCounter(const double *sorted_x_ptr, double trial, size_t n)
      : sorted_x_ptr(sorted_x_ptr), trial(trial), n(n), sumP(0) {}

  QnCounter(const QnCounter &other, Split)
      : sorted_x_ptr(other.sorted_x_ptr), trial(other.trial), n(other.n),
        sumP(0) {}

  void operator()(size_t begin, size_t end) {
    // Since the Zig implementation is O(n) for the whole array,
    // we can't easily split it into parallel segments while maintaining O(n).
    // However, we can split the work by giving each worker a range of 'i'.
    // But for Q_n, the sliding window is most efficient when done in one pass.
    // Let's implement a thread-safe way or just use the O(n) Zig call on the
    // whole array. Given that O(n) is extremely fast, serial execution of the
    // counting might be fine unless n is very large. To keep it parallel, each
    // worker will do a sub-segment.

    uint64_t count = 0;
    size_t j = 0;
    // Adjust j to the start of the range 'begin'
    while (j < begin && sorted_x_ptr[begin] - sorted_x_ptr[j] > trial) {
      j++;
    }

    for (size_t i = begin; i < end; ++i) {
      while (j < i && sorted_x_ptr[i] - sorted_x_ptr[j] > trial) {
        j++;
      }
      count += (i - j);
    }
    sumP += count;
  }

  void join(const QnCounter &other) { sumP += other.sumP; }
};

// For Qn refinement loop (whimed refinement)
// This is now partially handled by Zig for efficiency and branchlessness.
// We'll call whimed_ptr in C++ (it's O(n)) and use Zig for the bounds updates.

extern double whimed_ptr(double *a, int *iw, size_t n);

// [[Rcpp::export]]
double C_qn_fast(NumericVector x) {
  size_t n = x.size();
  if (n < 2)
    return NA_REAL;

  const double *x_ptr = x.begin();
  std::vector<double> sorted_x(n);
  std::copy(x_ptr, x_ptr + n, sorted_x.begin());
  if (n > 5000)
    tbb::parallel_sort(sorted_x.begin(), sorted_x.end());
  else
    std::sort(sorted_x.begin(), sorted_x.end());

  size_t h = n / 2 + 1;
  uint64_t k_target = (uint64_t)h * (h - 1) / 2;

  // For Qn, we use the advanced Johnson-Mizoguchi selector implemented in Zig.
  // This requires O(n) work space for the algorithm.
  std::vector<double> work(n);
  std::vector<int32_t> iweight(n);
  std::vector<int32_t> left(n);
  std::vector<int32_t> right(n);

  return zig_qn_jm_select(sorted_x.data(), n, k_target, work.data(),
                          iweight.data(), left.data(), right.data());
}

// Simplified whimed_ptr for completeness (keeping it in C++ for now to avoid
// porting selector)
double whimed_ptr(double *a, int *iw, size_t n) {
  if (n == 0)
    return 0.0;
  if (n == 1)
    return a[0];

  long long wtotal = 0;
  for (size_t i = 0; i < n; ++i)
    wtotal += iw[i];
  long long target = wtotal / 2;

  size_t left = 0, right = n - 1;
  while (left < right) {
    double pivot = a[left + (right - left) / 2];
    size_t i_lt = left;
    size_t i_eq = left;
    for (size_t j = left; j <= right; ++j) {
      if (a[j] < pivot) {
        std::swap(a[i_eq], a[j]);
        std::swap(iw[i_eq], iw[j]);
        if (i_eq != i_lt) {
          std::swap(a[i_lt], a[i_eq]);
          std::swap(iw[i_lt], iw[i_eq]);
        }
        i_lt++;
        i_eq++;
      } else if (a[j] == pivot) {
        std::swap(a[i_eq], a[j]);
        std::swap(iw[i_eq], iw[j]);
        i_eq++;
      }
    }
    long long w_lt = 0;
    for (size_t k = left; k < i_lt; ++k)
      w_lt += iw[k];
    long long w_eq = 0;
    for (size_t k = i_lt; k < i_eq; ++k)
      w_eq += iw[k];
    if (target < w_lt) {
      right = i_lt - 1;
    } else if (target < (w_lt + w_eq)) {
      return pivot;
    } else {
      target -= (w_lt + w_eq);
      left = i_eq;
    }
  }
  return a[left];
}
