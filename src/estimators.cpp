#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include "constants.h"
#include <RcppParallel.h>
#include <algorithm>
#include <cmath>
#include <tbb/parallel_sort.h>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// --- SN INNER KERNEL ---

inline double sn_index_select(const double *x, size_t n, size_t i, size_t h) {
  const double val_i = x[i];
  int32_t na = (int32_t)i;
  int32_t nb = (int32_t)(n - 1 - i);

  size_t k = h;
  int32_t a_start = 1;
  int32_t a_end = na;
  int32_t b_start = 1;
  int32_t b_end = nb;

  while (true) {
    if (a_start > a_end)
      return x[i + (size_t)(b_start + (int32_t)k - 1)] - val_i;
    if (b_start > b_end)
      return val_i - x[i - (size_t)(a_start + (int32_t)k - 1)];
    if (k == 1) {
      double valA = val_i - x[i - (size_t)a_start];
      double valB = x[i + (size_t)b_start] - val_i;
      return (valA < valB) ? valA : valB;
    }

    size_t half = k / 2;
    int32_t stepA = std::min((int32_t)half, a_end - a_start + 1);
    int32_t stepB = std::min((int32_t)half, b_end - b_start + 1);

    double valA = val_i - x[i - (size_t)(a_start + stepA - 1)];
    double valB = x[i + (size_t)(b_start + stepB - 1)] - val_i;

    if (valA < valB) {
      k -= (size_t)stepA;
      a_start += stepA;
    } else {
      k -= (size_t)stepB;
      b_start += stepB;
    }
  }
}

// --- UTILITIES ---

inline double lowmedian_ptr(double *arr, size_t n) {
  if (n == 0) return 0.0;
  size_t h = (n - 1) / 2;
  std::nth_element(arr, arr + h, arr + n);
  return arr[h];
}

// --- SN ESTIMATOR WORKER ---

struct SnWorker : public Worker {
  const double *sorted_x;
  size_t n;
  double *results;

  SnWorker(const double *sorted_x, size_t n, double *results)
      : sorted_x(sorted_x), n(n), results(results) {}

  void operator()(size_t begin, size_t end) {
    size_t h = n / 2;
    for (size_t i = begin; i < end; ++i) {
      results[i] = sn_index_select(sorted_x, n, i, h);
    }
  }
};

// [[Rcpp::export]]
double C_sn_fast(NumericVector x) {
  size_t n = x.size();
  if (n < 2) return NA_REAL;

  if (Rcpp::any(Rcpp::is_na(x)).is_true()) return NA_REAL;

  const double *x_ptr = x.begin();

  if (n <= 512) {
    double sorted_x[512];
    std::copy(x_ptr, x_ptr + n, sorted_x);
    std::sort(sorted_x, sorted_x + n);

    double inner_medians[512];
    size_t h_inner = n / 2;
    for (size_t i = 0; i < n; ++i) {
      inner_medians[i] = sn_index_select(sorted_x, n, i, h_inner);
    }
    double raw = lowmedian_ptr(inner_medians, n);
    return raw * 1.19259855312321 * get_sn_factor(n);
  }

  std::vector<double> sorted_x(n);
  std::copy(x_ptr, x_ptr + n, sorted_x.begin());

  if (n > 5000)
    tbb::parallel_sort(sorted_x.begin(), sorted_x.end());
  else
    std::sort(sorted_x.begin(), sorted_x.end());

  std::vector<double> inner_medians(n);
  SnWorker worker(sorted_x.data(), n, inner_medians.data());

  if (n > 1000)
    parallelFor(0, n, worker);
  else
    worker(0, n);

  double raw = lowmedian_ptr(inner_medians.data(), n);
  return raw * 1.19259855312321 * get_sn_factor(n);
}

// --- QN ESTIMATOR HELPERS ---

inline double whimed_cpp(double *a, int32_t *iw, size_t n, int64_t target) {
  if (n == 0) return 0.0;
  if (n == 1) return a[0];

  size_t l = 0, r = n - 1;
  int64_t t = target;

  while (l < r) {
    double pivot = a[l + (r - l) / 2];
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
    for (size_t idx = l; idx < i; ++idx) wleft += iw[idx];

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
      for (size_t idx = i; idx < i_eq; ++idx) weq += iw[idx];

      if (wleft + weq > t) return pivot;
      else {
        t -= (wleft + weq);
        l = i_eq;
      }
    }
  }
  return a[l];
}

struct QnCountWorker : public Worker {
  const double *x;
  size_t n;
  double trial;
  uint64_t sumP = 0;
  uint64_t sumQ = 0;

  QnCountWorker(const double *x, size_t n, double trial)
      : x(x), n(n), trial(trial) {}
  QnCountWorker(const QnCountWorker &other, Split)
      : x(other.x), n(other.n), trial(other.trial) {}

  void operator()(size_t begin, size_t end) {
    if (begin >= end) return;
    size_t i = begin;
    if (i == 0) i = 1;
    if (i >= end) return;

    size_t jp = std::upper_bound(x, x + i, x[i] - trial) - x;
    size_t jq = std::lower_bound(x, x + i, x[i] - trial) - x;

    for (; i < end; ++i) {
      while (jp < i && (x[i] - x[jp]) >= trial) jp++;
      sumP += (i - jp);
      while (jq < i && (x[i] - x[jq]) > trial) jq++;
      sumQ += (i - jq);
    }
  }

  void join(const QnCountWorker &other) {
    sumP += other.sumP;
    sumQ += other.sumQ;
  }
};

struct QnRefineWorker : public Worker {
  const double *x;
  size_t n;
  double trial;
  bool is_sumP;
  int32_t *bounds;

  QnRefineWorker(const double *x, size_t n, double trial, bool is_sumP, int32_t *bounds)
      : x(x), n(n), trial(trial), is_sumP(is_sumP), bounds(bounds) {}

  void operator()(size_t begin, size_t end) {
    if (begin >= end) return;
    size_t i = begin;
    if (i == 0) i = 1;
    if (i >= end) return;

    size_t j = is_sumP ? (std::upper_bound(x, x + i, x[i] - trial) - x)
                       : (std::lower_bound(x, x + i, x[i] - trial) - x);

    for (; i < end; ++i) {
      while (j < i && (is_sumP ? (x[i] - x[j]) >= trial : (x[i] - x[j]) > trial)) j++;
      if (is_sumP) {
        int32_t jj_bound = (int32_t)(i - j);
        if (jj_bound < bounds[i]) bounds[i] = jj_bound;
      } else {
        int32_t jj_bound = (int32_t)(i - j + 1);
        if (jj_bound > bounds[i]) bounds[i] = jj_bound;
      }
    }
  }
};

double qn_jm_select_cpp(const double *x, size_t n, uint64_t k_target, double *work, int32_t *iweight, int32_t *left, int32_t *right) {
  for (size_t i = 0; i < n; ++i) {
    left[i] = 1;
    right[i] = (int32_t)i;
  }

  uint64_t nL = 0;
  uint64_t nR = (uint64_t)n * (n - 1) / 2;

  while (nR - nL > n) {
    size_t m = 0;
    for (size_t i = 1; i < n; ++i) {
      if (left[i] <= right[i]) {
        int32_t w = right[i] - left[i] + 1;
        int32_t jj = left[i] + w / 2;
        work[m] = x[i] - x[i - jj];
        iweight[m] = w;
        m += 1;
      }
    }

    double trial = whimed_cpp(work, iweight, m, (int64_t)((nR - nL) / 2));

    uint64_t sumP = 0;
    uint64_t sumQ = 0;
    size_t jp = 0, jq = 0;
    for (size_t i = 1; i < n; ++i) {
      while (jp < i && (x[i] - x[jp]) >= trial) jp++;
      sumP += (i - jp);
      while (jq < i && (x[i] - x[jq]) > trial) jq++;
      sumQ += (i - jq);
    }

    if (k_target <= sumP) {
      size_t j = 0;
      for (size_t i = 1; i < n; ++i) {
        while (j < i && (x[i] - x[j]) >= trial) j++;
        int32_t jj_bound = (int32_t)(i - j);
        if (jj_bound < right[i]) right[i] = jj_bound;
      }
      nR = sumP;
    } else if (k_target > sumQ) {
      size_t j = 0;
      for (size_t i = 1; i < n; ++i) {
        while (j < i && (x[i] - x[j]) > trial) j++;
        int32_t jj_bound = (int32_t)(i - j + 1);
        if (jj_bound > left[i]) left[i] = jj_bound;
      }
      nL = sumQ;
    } else {
      return trial;
    }
  }

  std::vector<double> final_diffs;
  final_diffs.reserve(nR - nL);
  for (size_t i = 1; i < n; ++i) {
    for (int32_t jj = left[i]; jj <= right[i]; ++jj) {
      final_diffs.push_back(x[i] - x[i - jj]);
    }
  }
  std::nth_element(final_diffs.begin(), final_diffs.begin() + (k_target - nL - 1), final_diffs.end());
  return final_diffs[k_target - nL - 1];
}

// [[Rcpp::export]]
double C_qn_fast(NumericVector x) {
  size_t n = x.size();
  if (n < 2) return NA_REAL;

  if (Rcpp::any(Rcpp::is_na(x)).is_true()) return NA_REAL;

  const double *x_ptr = x.begin();

  if (n <= 128) {
    double sorted_x[128];
    std::copy(x_ptr, x_ptr + n, sorted_x);
    std::sort(sorted_x, sorted_x + n);

    size_t num_pairs = n * (n - 1) / 2;
    double diffs[128 * 127 / 2];
    size_t k = 0;
    for (size_t i = 1; i < n; ++i) {
      for (size_t j = 0; j < i; ++j) {
        diffs[k++] = sorted_x[i] - sorted_x[j];
      }
    }
    size_t h = n / 2 + 1;
    size_t k_target = h * (h - 1) / 2;
    std::nth_element(diffs, diffs + k_target - 1, diffs + num_pairs);
    double raw = diffs[k_target - 1];
    return raw * 2.21914446598508 * get_qn_factor(n);
  }

  std::vector<double> sorted_x(n);
  std::copy(x_ptr, x_ptr + n, sorted_x.begin());
  if (n > 5000)
    tbb::parallel_sort(sorted_x.begin(), sorted_x.end());
  else
    std::sort(sorted_x.begin(), sorted_x.end());

  size_t h = n / 2 + 1;
  uint64_t k_target = (uint64_t)h * (h - 1) / 2;

  if (n <= 2000) {
    std::vector<double> work(n);
    std::vector<int32_t> iweight(n);
    std::vector<int32_t> left(n);
    std::vector<int32_t> right(n);
    double raw = qn_jm_select_cpp(sorted_x.data(), n, k_target, work.data(),
                                  iweight.data(), left.data(), right.data());
    return raw * 2.21914446598508 * get_qn_factor(n);
  }

  std::vector<int32_t> left(n, 1);
  std::vector<int32_t> right(n);
  for (size_t i = 0; i < n; ++i) right[i] = (int32_t)i;

  std::vector<double> work(n);
  std::vector<int32_t> iweight(n);
  uint64_t nL = 0;
  uint64_t nR = (uint64_t)n * (n - 1) / 2;

  while (nR - nL > n) {
    size_t m = 0;
    for (size_t i = 1; i < n; ++i) {
      if (left[i] <= right[i]) {
        int32_t w = right[i] - left[i] + 1;
        int32_t jj = left[i] + w / 2;
        work[m] = sorted_x[i] - sorted_x[i - jj];
        iweight[m] = w;
        m += 1;
      }
    }

    double trial = whimed_cpp(work.data(), iweight.data(), m, (int64_t)((nR - nL) / 2));

    QnCountWorker countWorker(sorted_x.data(), n, trial);
    parallelReduce(1, n, countWorker);

    if (k_target <= countWorker.sumP) {
      QnRefineWorker refineWorker(sorted_x.data(), n, trial, true, right.data());
      parallelFor(1, n, refineWorker);
      nR = countWorker.sumP;
    } else if (k_target > countWorker.sumQ) {
      QnRefineWorker refineWorker(sorted_x.data(), n, trial, false, left.data());
      parallelFor(1, n, refineWorker);
      nL = countWorker.sumQ;
    } else {
      return trial * 2.21914446598508 * get_qn_factor(n);
    }
  }

  std::vector<double> final_diffs;
  final_diffs.reserve(nR - nL);
  for (size_t i = 1; i < n; ++i) {
    for (int32_t jj = left[i]; jj <= right[i]; ++jj) {
      final_diffs.push_back(sorted_x[i] - sorted_x[i - jj]);
    }
  }
  std::nth_element(final_diffs.begin(), final_diffs.begin() + (k_target - nL - 1), final_diffs.end());
  double raw = final_diffs[k_target - nL - 1];
  return raw * 2.21914446598508 * get_qn_factor(n);
}
