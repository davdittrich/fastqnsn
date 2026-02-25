#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(BH)]]
#include "constants.h"
#include "sort_utils.h"
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <memory>
#include <new>
#include <type_traits>

using namespace Rcpp;
using namespace RcppParallel;

// --- FLOYD-RIVEST SELECTION ---
// Faster than std::nth_element for average cases (~30% fewer comparisons)
template <typename Iter>
void floyd_rivest_select(Iter begin, Iter nth, Iter end) {
  using T = typename std::iterator_traits<Iter>::value_type;
  while (end - begin > 5) {
    auto n = end - begin;
    auto k = nth - begin;
    double nr = (double)n;
    double z = std::log(nr);
    double s = 0.5 * std::exp(2.0 * z / 3.0);
    double sd = 0.5 * std::sqrt(z * s * (nr - s) / nr) *
                (2.0 * k - n + 1 < 0 ? -1.0 : 1.0);
    auto newl = begin + (std::max)((decltype(k))0,
                                   (decltype(k))(k - (double)k * s / nr + sd));
    auto newr = begin + (std::min)((decltype(k))(n - 1),
                                   (decltype(k))(k + (nr - k) * s / nr + sd));
    floyd_rivest_select(begin + (newl - begin), nth,
                        begin + (newr - begin) + 1);
    begin = newl;
    end = newr + 1;
  }
  // Insertion sort for small partitions
  for (auto i = begin + 1; i < end; ++i) {
    T val = std::move(*i);
    auto j = i;
    while (j > begin && *(j - 1) > val) {
      *j = std::move(*(j - 1));
      --j;
    }
    *j = std::move(val);
  }
}

// --- UTILITIES ---

template <typename T> inline double lowmedian_ptr(T *arr, size_t n) {
  if (n == 0)
    return 0.0;
  size_t h = (n - 1) / 2;
  floyd_rivest_select(arr, arr + h, arr + n);
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

    int32_t low = L;
    int32_t high = L_max_limit;
    int32_t best_L = L;
    while (low <= high) {
      int32_t mid = low + (high - low) / 2;
      if (mid == L_max_limit) {
        best_L = mid;
        break;
      }
      double v_mid = std::max(sorted_x[i] - sorted_x[mid],
                              sorted_x[mid + h] - sorted_x[i]);
      double v_next = std::max(sorted_x[i] - sorted_x[mid + 1],
                               sorted_x[mid + 1 + h] - sorted_x[i]);
      if (v_mid > v_next) {
        low = mid + 1;
        best_L = mid + 1;
      } else {
        high = mid - 1;
      }
    }
    L = best_L;

    for (; i < static_cast<int32_t>(end); ++i) {
      int32_t L_min = std::max(0, i - h);
      int32_t L_max = std::min(i, static_cast<int32_t>(n) - 1 - h);
      if (L < L_min)
        L = L_min;

      T candidate =
          std::max(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]);
      while (L < L_max) {
        T next = std::max(sorted_x[i] - sorted_x[L + 1],
                          sorted_x[L + 1 + h] - sorted_x[i]);
        if (candidate <= next)
          break;
        L++;
        candidate = next;
      }
      results[i] = candidate;
    }
  }
};

template <typename T> double C_sn_impl(const T *x_ptr, size_t n) {
  if (n < 2)
    return NA_REAL;
  if (n > 6060000000ULL)
    Rcpp::stop("fastqnsn Error: sample size n > 6.06 * 10^9 natively overflows "
               "64-bit pair boundaries. 128-bit architecture required.");

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
      T candidate =
          std::max(sorted_x[i] - sorted_x[L], sorted_x[L + h] - sorted_x[i]);
      while (L < L_max) {
        T next = std::max(sorted_x[i] - sorted_x[L + 1],
                          sorted_x[L + 1 + h] - sorted_x[i]);
        if (candidate <= next)
          break;
        L++;
        candidate = next;
      }
      inner_medians[i] = candidate;
    }
    double raw = lowmedian_ptr(inner_medians, n);
    return raw * CONST_SN * get_sn_factor(n);
  }

  // Arena allocation: sorted_x(n*T) + inner_medians(n*T)
  std::unique_ptr<T[]> sn_arena;
  try {
    sn_arena = std::make_unique<T[]>(2 * n);
  } catch (const std::bad_alloc &e) {
    Rcpp::stop(
        "fastqnsn Out of Memory: failed to allocate %zu bytes for Sn arena.",
        2 * n * sizeof(T));
  }
  T *sorted_x = sn_arena.get();
  T *inner_medians = sn_arena.get() + n;

  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (!std::isfinite(x_ptr[i]))
        return NA_REAL;
    }
    sorted_x[i] = x_ptr[i];
  }

  fastqnsn::optimized_sort(sorted_x, sorted_x + n);

  SnWorker<T> worker(sorted_x, n, inner_medians);

  if (n > 10000)
    parallelFor(0, n, worker, 2048);
  else
    worker(0, n);

  double raw = lowmedian_ptr(inner_medians, n);
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
      while (__builtin_expect(jp < i && (double)x[jp] <= target, 0))
        jp++;
      while (__builtin_expect(jq < jp && (double)x[jq] < target, 0))
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
        while (__builtin_expect(j < i && (double)x[j] <= target, 0))
          j++;
        int32_t jj_bound = static_cast<int32_t>(i - j);
        if (jj_bound < bounds[i])
          bounds[i] = jj_bound;
      } else {
        while (__builtin_expect(j < i && (double)x[j] < target, 0))
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
  if (n > 6060000000ULL)
    Rcpp::stop("fastqnsn Error: sample size n > 6.06 * 10^9 natively overflows "
               "64-bit pair boundaries. 128-bit architecture required.");

  if (n <= 3000) {
    std::unique_ptr<T[]> sorted_x(new T[n]);
    for (size_t i = 0; i < n; i++) {
      if constexpr (std::is_floating_point_v<T>) {
        if (!std::isfinite(x_ptr[i]))
          return NA_REAL;
      }
      sorted_x[i] = x_ptr[i];
    }
    std::sort(sorted_x.get(), sorted_x.get() + n);

    size_t num_pairs = n * (n - 1) / 2;
    std::unique_ptr<double[]> diffs(new double[num_pairs]);
    size_t k_idx = 0;
    for (size_t i = 1; i < n; ++i) {
      for (size_t j = 0; j < i; ++j) {
        diffs[k_idx++] = (double)sorted_x[i] - (double)sorted_x[j];
      }
    }
    size_t h = n / 2 + 1;
    size_t k_target = h * (h - 1) / 2;
    floyd_rivest_select(diffs.get(), diffs.get() + k_target - 1,
                        diffs.get() + num_pairs);
    double raw = diffs[k_target - 1];
    return raw * CONST_QN * get_qn_factor(n);
  }

  // Arena allocation: sorted_x(n*T) + work(n*double) + iweight(n*int32) +
  // left(n*int32) + right(n*int32)
  size_t arena_bytes =
      n * sizeof(T) + n * sizeof(double) + 3 * n * sizeof(int32_t);
  std::unique_ptr<char[]> arena;
  try {
    arena = std::make_unique<char[]>(arena_bytes);
  } catch (const std::bad_alloc &e) {
    Rcpp::stop(
        "fastqnsn Out of Memory: failed to allocate %zu bytes for Qn arena.",
        arena_bytes);
  }
  char *ptr = arena.get();
  T *sorted_x = reinterpret_cast<T *>(ptr);
  ptr += n * sizeof(T);
  double *work = reinterpret_cast<double *>(ptr);
  ptr += n * sizeof(double);
  int32_t *iweight = reinterpret_cast<int32_t *>(ptr);
  ptr += n * sizeof(int32_t);
  int32_t *left = reinterpret_cast<int32_t *>(ptr);
  ptr += n * sizeof(int32_t);
  int32_t *right = reinterpret_cast<int32_t *>(ptr);

  for (size_t i = 0; i < n; i++) {
    if constexpr (std::is_floating_point_v<T>) {
      if (!std::isfinite(x_ptr[i]))
        return NA_REAL;
    }
    sorted_x[i] = x_ptr[i];
  }

  fastqnsn::optimized_sort(sorted_x, sorted_x + n);

  size_t h = n / 2 + 1;
  uint64_t k_target = (uint64_t)h * (h - 1) / 2;

  std::fill_n(left, n, 1);
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
        work[m] = (float)((double)sorted_x[i] - (double)sorted_x[i - jj]);
        iweight[m] = w;
        m += 1;
      }
    }

    double trial =
        whimed_cpp(work, iweight, m, static_cast<int64_t>((nR - nL) / 2));

    QnCountWorker<T> countWorker(sorted_x, n, trial);
    if (n > 10000)
      parallelReduce(1, n, countWorker, 2048);
    else
      countWorker(1, n);

    if (k_target <= countWorker.sumP) {
      QnRefineWorker<T> refineWorker(sorted_x, n, trial, true, right);
      if (n > 10000)
        parallelFor(1, n, refineWorker, 2048);
      else
        refineWorker(1, n);
      nR = countWorker.sumP;
    } else if (k_target > countWorker.sumQ) {
      QnRefineWorker<T> refineWorker(sorted_x, n, trial, false, left);
      if (n > 10000)
        parallelFor(1, n, refineWorker, 2048);
      else
        refineWorker(1, n);
      nL = countWorker.sumQ;
    } else {
      return trial * CONST_QN * get_qn_factor(n);
    }
  }

  std::unique_ptr<double[]> final_diffs;
  try {
    final_diffs = std::make_unique<double[]>(nR - nL);
  } catch (const std::bad_alloc &e) {
    Rcpp::stop("fastqnsn Out of Memory: failed to allocate %zu bytes for final "
               "Qn diffs.",
               (size_t)(nR - nL) * sizeof(double));
  }
  size_t fd_idx = 0;
  for (size_t i = 1; i < n; ++i) {
    for (int32_t jj = left[i]; jj <= right[i]; ++jj) {
      final_diffs[fd_idx++] = (double)sorted_x[i] - (double)sorted_x[i - jj];
    }
  }
  floyd_rivest_select(final_diffs.get(),
                      final_diffs.get() + (k_target - nL - 1),
                      final_diffs.get() + (nR - nL));
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
