#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]
#include "constants.h"
#include "sort_utils.h"
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <type_traits>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// --- UTILITIES ---

template <typename T> inline double lowmedian_ptr(T *arr, size_t n) {
  if (n == 0)
    return 0.0;
  size_t h = (n - 1) / 2;
  std::nth_element(arr, arr + h, arr + n);
  return static_cast<double>(arr[h]);
}

// --- SN ESTIMATOR WORKER ---

template <typename T> struct SnWorker : public Worker {
  const T *sorted_x;
  size_t n;
  T *results;

  SnWorker(const T *sorted_x, size_t n, T *results)
      : sorted_x(sorted_x), n(n), results(results) {}

  void operator()(size_t begin, size_t end) {
    if (begin >= end)
      return;
    int32_t h = static_cast<int32_t>(n / 2);

    int32_t i = static_cast<int32_t>(begin);
    int32_t L = std::max(0, i - h);
    int32_t L_max_limit = std::min(i, static_cast<int32_t>(n) - 1 - h);

    while (L < L_max_limit &&
           std::max(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]) >
               std::max(sorted_x[i] - sorted_x[L + 1],
                        sorted_x[L + 1 + h] - sorted_x[i])) {
      L++;
    }

    for (; i < static_cast<int32_t>(end); ++i) {
      int32_t L_min = std::max(0, i - h);
      int32_t L_max = std::min(i, static_cast<int32_t>(n) - 1 - h);
      if (L < L_min)
        L = L_min;

      while (L < L_max && std::max(sorted_x[i] - sorted_x[L],
                                   sorted_x[L + h] - sorted_x[i]) >
                              std::max(sorted_x[i] - sorted_x[L + 1],
                                       sorted_x[L + 1 + h] - sorted_x[i])) {
        L++;
      }
      results[i] =
          std::max(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]);
    }
  }
};

template <typename T> double C_sn_impl(const T *x_ptr, size_t n) {
  if (n < 2)
    return NA_REAL;

  if (n <= 2000) {
    T sorted_x[2000];
    for (size_t i = 0; i < n; i++) {
      if constexpr (std::is_floating_point_v<T>) {
        if (!std::isfinite(x_ptr[i]))
          return NA_REAL;
      }
      sorted_x[i] = x_ptr[i];
    }
    fastqnsn::optimized_sort(sorted_x, sorted_x + n);

    T inner_medians[2000];
    int32_t h = static_cast<int32_t>(n / 2);
    int32_t L = 0;
    for (int32_t i = 0; i < static_cast<int32_t>(n); ++i) {
      int32_t L_min = std::max(0, i - h);
      int32_t L_max = std::min(i, static_cast<int32_t>(n) - 1 - h);
      if (L < L_min)
        L = L_min;
      while (L < L_max && std::max(sorted_x[i] - sorted_x[L],
                                   sorted_x[L + h] - sorted_x[i]) >
                              std::max(sorted_x[i] - sorted_x[L + 1],
                                       sorted_x[L + 1 + h] - sorted_x[i])) {
        L++;
      }
      inner_medians[i] =
          std::max(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]);
    }
    double raw = lowmedian_ptr(inner_medians, n);
    return raw * CONST_SN * get_sn_factor(n);
  }

  std::vector<T> sorted_x(n);
  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (!std::isfinite(x_ptr[i]))
        return NA_REAL;
    }
    sorted_x[i] = x_ptr[i];
  }

  fastqnsn::optimized_sort(sorted_x.begin(), sorted_x.end());

  std::vector<T> inner_medians(n);
  SnWorker<T> worker(sorted_x.data(), n, inner_medians.data());

  if (n > 10000)
    parallelFor(0, n, worker, std::max<size_t>(10000, n / 4));
  else
    worker(0, n);

  double raw = lowmedian_ptr(inner_medians.data(), n);
  return raw * CONST_SN * get_sn_factor(n);
}

// --- QN ESTIMATOR HELPERS ---

template <typename T>
inline T whimed_cpp(T *a, int32_t *iw, size_t n, int64_t target) {
  if (n == 0)
    return T(0);
  if (n == 1)
    return a[0];

  size_t l = 0, r = n - 1;
  int64_t t = target;

  while (l < r) {
    T pivot = a[l + (r - l) / 2];
    size_t i = l, j = l;
    while (j <= r) {
      if (a[j] < pivot) {
        std::swap(a[i], a[j]);
        std::swap(iw[i], iw[j]);
        i++;
      }
      j++;
    }
    int64_t wleft = 0;
    for (size_t idx = l; idx < i; ++idx)
      wleft += iw[idx];

    if (wleft > t) {
      r = (i > l) ? i - 1 : l;
    } else {
      size_t i_eq = i, j_eq = i;
      while (j_eq <= r) {
        if (a[j_eq] == pivot) {
          std::swap(a[i_eq], a[j_eq]);
          std::swap(iw[i_eq], iw[j_eq]);
          i_eq++;
        }
        j_eq++;
      }
      int64_t weq = 0;
      for (size_t idx = i; idx < i_eq; ++idx)
        weq += iw[idx];

      if (wleft + weq > t)
        return pivot;
      else {
        t -= (wleft + weq);
        l = i_eq;
      }
    }
  }
  return a[l];
}

template <typename T> struct QnCountWorker : public Worker {
  const T *x;
  size_t n;
  double trial;
  uint64_t sumP = 0;
  uint64_t sumQ = 0;

  QnCountWorker(const T *x, size_t n, double trial)
      : x(x), n(n), trial(trial) {}
  QnCountWorker(const QnCountWorker &other, Split)
      : x(other.x), n(other.n), trial(other.trial) {}

  void operator()(size_t begin, size_t end) {
    if (begin >= end)
      return;
    size_t i = begin;
    if (i == 0)
      i = 1;
    if (i >= end)
      return;

    size_t jp = std::upper_bound(x, x + i, (double)x[i] - trial) - x;
    size_t jq = std::lower_bound(x, x + i, (double)x[i] - trial) - x;

    for (; i < end; ++i) {
      double target = (double)x[i] - trial;
      while (jp < i && (double)x[jp] <= target)
        jp++;
      while (jq < jp && (double)x[jq] < target)
        jq++;
      sumP += (i - jp);
      sumQ += (i - jq);
    }
  }

  void join(const QnCountWorker &other) {
    sumP += other.sumP;
    sumQ += other.sumQ;
  }
};

template <typename T> struct QnRefineWorker : public Worker {
  const T *x;
  size_t n;
  double trial;
  bool is_sumP;
  int32_t *bounds;

  QnRefineWorker(const T *x, size_t n, double trial, bool is_sumP,
                 int32_t *bounds)
      : x(x), n(n), trial(trial), is_sumP(is_sumP), bounds(bounds) {}

  void operator()(size_t begin, size_t end) {
    if (begin >= end)
      return;
    size_t i = begin;
    if (i == 0)
      i = 1;
    if (i >= end)
      return;

    size_t j = is_sumP ? (std::upper_bound(x, x + i, (double)x[i] - trial) - x)
                       : (std::lower_bound(x, x + i, (double)x[i] - trial) - x);

    for (; i < end; ++i) {
      double target = (double)x[i] - trial;
      if (is_sumP) {
        while (j < i && (double)x[j] <= target)
          j++;
        int32_t jj_bound = static_cast<int32_t>(i - j);
        if (jj_bound < bounds[i])
          bounds[i] = jj_bound;
      } else {
        while (j < i && (double)x[j] < target)
          j++;
        int32_t jj_bound = static_cast<int32_t>(i - j + 1);
        if (jj_bound > bounds[i])
          bounds[i] = jj_bound;
      }
    }
  }
};

template <typename T> double C_qn_impl(const T *x_ptr, size_t n) {
  if (n < 2)
    return NA_REAL;

  if (n <= 45) {
    T sorted_x[45];
    for (size_t i = 0; i < n; i++) {
      if constexpr (std::is_floating_point_v<T>) {
        if (!std::isfinite(x_ptr[i]))
          return NA_REAL;
      }
      sorted_x[i] = x_ptr[i];
    }
    fastqnsn::optimized_sort(sorted_x, sorted_x + n);

    size_t num_pairs = n * (n - 1) / 2;

    double diffs[990]; // 45 * 44 / 2 = 990
    size_t k_idx = 0;
    for (size_t i = 1; i < n; ++i) {
      for (size_t j = 0; j < i; ++j) {
        diffs[k_idx++] = (double)sorted_x[i] - (double)sorted_x[j];
      }
    }
    size_t h = n / 2 + 1;
    size_t k_target = h * (h - 1) / 2;
    std::nth_element(diffs, diffs + k_target - 1, diffs + num_pairs);
    double raw = diffs[k_target - 1];
    return raw * CONST_QN * get_qn_factor(n);
  }

  std::vector<T> sorted_x(n);
  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (!std::isfinite(x_ptr[i]))
        return NA_REAL;
    }
    sorted_x[i] = x_ptr[i];
  }

  fastqnsn::optimized_sort(sorted_x.begin(), sorted_x.end());

  size_t h = n / 2 + 1;
  uint64_t k_target = (uint64_t)h * (h - 1) / 2;

  std::vector<double> work(n);
  std::vector<int32_t> iweight(n);
  std::vector<int32_t> left(n, 1);
  std::vector<int32_t> right(n);
  for (size_t i = 0; i < n; ++i)
    right[i] = static_cast<int32_t>(i);

  uint64_t nL = 0;
  uint64_t nR = (uint64_t)n * (n - 1) / 2;

  while (nR - nL > n) {
    size_t m = 0;
    for (size_t i = 1; i < n; ++i) {
      if (left[i] <= right[i]) {
        int32_t w = right[i] - left[i] + 1;
        int32_t jj = left[i] + w / 2;
        work[m] = (double)sorted_x[i] - (double)sorted_x[i - jj];
        iweight[m] = w;
        m += 1;
      }
    }

    double trial = whimed_cpp(work.data(), iweight.data(), m,
                              static_cast<int64_t>((nR - nL) / 2));

    QnCountWorker<T> countWorker(sorted_x.data(), n, trial);
    if (n > 10000)
      parallelReduce(1, n, countWorker, std::max<size_t>(10000, n / 4));
    else
      countWorker(1, n);

    if (k_target <= countWorker.sumP) {
      QnRefineWorker<T> refineWorker(sorted_x.data(), n, trial, true,
                                     right.data());
      if (n > 10000)
        parallelFor(1, n, refineWorker, std::max<size_t>(10000, n / 4));
      else
        refineWorker(1, n);
      nR = countWorker.sumP;
    } else if (k_target > countWorker.sumQ) {
      QnRefineWorker<T> refineWorker(sorted_x.data(), n, trial, false,
                                     left.data());
      if (n > 10000)
        parallelFor(1, n, refineWorker, std::max<size_t>(10000, n / 4));
      else
        refineWorker(1, n);
      nL = countWorker.sumQ;
    } else {
      return trial * CONST_QN * get_qn_factor(n);
    }
  }

  std::vector<double> final_diffs;
  final_diffs.reserve(nR - nL);
  for (size_t i = 1; i < n; ++i) {
    for (int32_t jj = left[i]; jj <= right[i]; ++jj) {
      final_diffs.push_back((double)sorted_x[i] - (double)sorted_x[i - jj]);
    }
  }
  std::nth_element(final_diffs.begin(),
                   final_diffs.begin() + (k_target - nL - 1),
                   final_diffs.end());
  double raw = final_diffs[k_target - nL - 1];
  return raw * CONST_QN * get_qn_factor(n);
}

// --- R EXPORTS ---

// [[Rcpp::export]]
double C_sn_fast(NumericVector x) { return C_sn_impl(x.begin(), x.size()); }

// [[Rcpp::export]]
double C_sn_int_fast(IntegerVector x) { return C_sn_impl(x.begin(), x.size()); }

// [[Rcpp::export]]
double C_qn_fast(NumericVector x) { return C_qn_impl(x.begin(), x.size()); }

// [[Rcpp::export]]
double C_qn_int_fast(IntegerVector x) { return C_qn_impl(x.begin(), x.size()); }
