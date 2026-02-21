#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <tbb/parallel_sort.h>
#include <algorithm>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace RcppParallel;

// --- UTILITIES ---

inline double lowmedian_ptr(double* arr, size_t n) {
    if (n == 0) return 0.0;
    size_t m = (size_t)(std::floor(((double)n + 1.0) / 2.0) - 1.0);
    std::nth_element(arr, arr + m, arr + n);
    return arr[m];
}

// --- NAIVE SCALAR IMPLEMENTATIONS (FOR SMALL N) ---

double C_sn_naive(const double* x, size_t n) {
    std::vector<double> inner_medians(n);
    std::vector<double> diffs;
    diffs.reserve(n - 1);
    for (size_t i = 0; i < n; ++i) {
        diffs.clear();
        for (size_t j = 0; j < n; ++j) {
            if (i != j) {
                double d = std::abs(x[i] - x[j]);
                if (std::isfinite(d)) diffs.push_back(d);
            }
        }
        if (diffs.empty()) {
            inner_medians[i] = 0.0;
            continue;
        }
        size_t m = (diffs.size() - 1) / 2;
        std::nth_element(diffs.begin(), diffs.begin() + m, diffs.end());
        inner_medians[i] = diffs[m];
    }
    return lowmedian_ptr(inner_medians.data(), n);
}

double C_qn_naive(const double* x, size_t n) {
    size_t count = n * (n - 1) / 2;
    std::vector<double> diffs;
    diffs.reserve(count);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            double d = std::abs(x[i] - x[j]);
            if (std::isfinite(d)) diffs.push_back(d);
        }
    }
    size_t h = n / 2 + 1;
    size_t k = h * (h - 1) / 2 - 1; 
    if (diffs.empty()) return 0.0;
    if (k >= diffs.size()) k = diffs.size() - 1;
    std::nth_element(diffs.begin(), diffs.begin() + k, diffs.end());
    return diffs[k];
}

// --- SN ESTIMATOR WORKER (LARGE N) ---

struct SnWorker : public Worker {
    const double* sorted_x;
    size_t n;
    double* results;

    SnWorker(const double* sorted_x, size_t n, double* results)
        : sorted_x(sorted_x), n(n), results(results) {}

    void operator()(size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            if (i == 0) { results[0] = sorted_x[n / 2] - sorted_x[0]; continue; }
            if (i == n - 1) { results[n - 1] = sorted_x[n - 1] - sorted_x[(n + 1) / 2 - 1]; continue; }

            size_t nA = i;
            size_t nB = n - 1 - i;
            size_t target_rank = (n - 1) / 2;
            
            // Lambdas for O(log n) selection
            auto getA = [&](size_t idx) { return sorted_x[i] - sorted_x[i - 1 - idx]; };
            auto getB = [&](size_t idx) { return sorted_x[i + 1 + idx] - sorted_x[i]; };

            if (i < (n + 1) / 2) {
                size_t k = target_rank + 1, lA = 0, rA = nA, lB = 0, rB = nB;
                while (true) {
                    if (lA == rA) { results[i] = getB(lB + k - 1); break; }
                    if (lB == rB) { results[i] = getA(lA + k - 1); break; }
                    if (k == 1) { results[i] = (std::min)(getA(lA), getB(lB)); break; }
                    size_t mid = k / 2;
                    size_t stepA = (std::min)(mid, rA - lA), stepB = (std::min)(mid, rB - lB);
                    if (getA(lA + stepA - 1) < getB(lB + stepB - 1)) { lA += stepA; k -= stepA; }
                    else { lB += stepB; k -= stepB; }
                }
            } else {
                size_t k = target_rank + 1, lA = 0, rA = nB, lB = 0, rB = nA;
                // Swap roles of A and B for symmetry
                auto getA_sym = [&](size_t idx) { return sorted_x[i + 1 + idx] - sorted_x[i]; };
                auto getB_sym = [&](size_t idx) { return sorted_x[i] - sorted_x[i - 1 - idx]; };
                while (true) {
                    if (lA == rA) { results[i] = getB_sym(lB + k - 1); break; }
                    if (lB == rB) { results[i] = getA_sym(lA + k - 1); break; }
                    if (k == 1) { results[i] = (std::min)(getA_sym(lA), getB_sym(lB)); break; }
                    size_t mid = k / 2;
                    size_t stepA = (std::min)(mid, rA - lA), stepB = (std::min)(mid, rB - lB);
                    if (getA_sym(lA + stepA - 1) < getB_sym(lB + stepB - 1)) { lA += stepA; k -= stepA; }
                    else { lB += stepB; k -= stepB; }
                }
            }
        }
    }
};

// [[Rcpp::export]]
double C_sn_fast(NumericVector x) {
    size_t n = x.size();
    if (n < 2) return NA_REAL;
    if (n < 50) return C_sn_naive(x.begin(), n);

    std::vector<double> sorted_x = as<std::vector<double>>(x);
    if (n > 10000) tbb::parallel_sort(sorted_x.begin(), sorted_x.end());
    else std::sort(sorted_x.begin(), sorted_x.end());
    
    std::vector<double> inner_medians(n);
    SnWorker worker(sorted_x.data(), n, inner_medians.data());
    if (n > 500) parallelFor(0, n, worker);
    else worker(0, n);
    
    return lowmedian_ptr(inner_medians.data(), n);
}

// --- QN ESTIMATOR HELPERS ---

double whimed_ptr(double* a, int* iw, size_t n) {
    if (n == 0) return 0.0;
    if (n == 1) return a[0];
    long long wtotal = 0;
    for (size_t i = 0; i < n; ++i) wtotal += iw[i];
    long long target = wtotal / 2;
    size_t left = 0, right = n - 1;
    while (left < right) {
        double pivot = a[left + (right - left) / 2];
        size_t i = left;
        for (size_t j = left; j <= right; ++j) {
            if (a[j] < pivot) {
                std::swap(a[i], a[j]); std::swap(iw[i], iw[j]);
                i++;
            }
        }
        size_t split = i;
        long long wleft = 0;
        for (size_t k = left; k < split; ++k) wleft += iw[k];
        if (wleft > target) { right = split - 1; }
        else { target -= wleft; left = split; }
    }
    return a[left];
}

struct QnCounter : public Worker {
    const double* sorted_x;
    double trial;
    size_t n;
    long long sumP;
    long long sumQ;

    QnCounter(const double* sorted_x, double trial, size_t n)
        : sorted_x(sorted_x), trial(trial), n(n), sumP(0), sumQ(0) {}

    QnCounter(const QnCounter& other, Split)
        : sorted_x(other.sorted_x), trial(other.trial), n(other.n), sumP(0), sumQ(0) {}

    void operator()(size_t begin, size_t end) {
        for (size_t i = begin; i < end; ++i) {
            int j_p = 0;
            while (j_p < (int)n && (sorted_x[i] - sorted_x[n - 1 - j_p]) < trial) j_p++;
            sumP += j_p;
            int j_q = n;
            while (j_q > 0 && (sorted_x[i] - sorted_x[n - j_q]) > trial) j_q--;
            sumQ += j_q;
        }
    }

    void join(const QnCounter& other) {
        sumP += other.sumP;
        sumQ += other.sumQ;
    }
};

// [[Rcpp::export]]
double C_qn_fast(NumericVector x) {
    size_t n = x.size();
    if (n < 2) return NA_REAL;
    if (n < 100) return C_qn_naive(x.begin(), n);

    std::vector<double> sorted_x = as<std::vector<double>>(x);
    if (n > 10000) tbb::parallel_sort(sorted_x.begin(), sorted_x.end());
    else std::sort(sorted_x.begin(), sorted_x.end());
    
    size_t h = n / 2 + 1;
    long long k_target = (long long)h * (h - 1) / 2;
    std::vector<int> left(n), right(n, n);
    for(size_t i = 0; i < n; ++i) left[i] = n - i + 1;
    
    long long nL = (long long)n * (n + 1) / 2;
    long long nR = (long long)n * n;
    long long knew = k_target + nL; 
    
    double qn_val = 0;
    bool found = false;

    // Optimization: Reuse vectors to avoid O(n) allocations in refinement loop
    std::vector<double> work;
    std::vector<int> weights;
    work.reserve(n);
    weights.reserve(n);

    for (int iter = 0; iter < 100; ++iter) {
        if (nR - nL <= (long long)n) break;
        work.clear();
        weights.clear();
        for (size_t i = 1; i < n; ++i) {
            if (left[i] <= right[i]) {
                int w = right[i] - left[i] + 1;
                work.push_back(sorted_x[i] - sorted_x[n - (left[i] + w/2)]);
                weights.push_back(w);
            }
        }
        if (work.empty()) break;
        double trial = whimed_ptr(work.data(), weights.data(), work.size());
        
        QnCounter counter(sorted_x.data(), trial, n);
        if (n > 1000) parallelReduce(0, n, counter);
        else counter(0, n);
        
        if (knew <= counter.sumP) {
            for (int i = 0; i < (int)n; ++i) {
                int jj = 0;
                while (jj < (int)n && (sorted_x[i] - sorted_x[n - 1 - jj]) < trial) jj++;
                right[i] = jj;
            }
            nR = counter.sumP;
        } else if (knew > counter.sumQ) {
            for (int i = 0; i < (int)n; ++i) {
                int jj = n;
                while (jj > 0 && (sorted_x[i] - sorted_x[n - jj]) > trial) jj--;
                left[i] = jj + 1;
            }
            nL = counter.sumQ;
        } else {
            qn_val = trial; found = true; break;
        }
    }
    
    if (!found) {
        std::vector<double> candidates;
        for (size_t i = 1; i < n; ++i) {
            for (int jj_idx = left[i]; jj_idx <= right[i]; ++jj_idx) {
                double val = sorted_x[i] - sorted_x[n - jj_idx];
                if (std::isfinite(val)) candidates.push_back(val);
            }
        }
        if (candidates.empty()) return 0.0;
        size_t rank = knew - nL - 1;
        if (rank >= candidates.size()) rank = candidates.size() - 1;
        std::nth_element(candidates.begin(), candidates.begin() + rank, candidates.end());
        qn_val = candidates[rank];
    }
    return qn_val;
}
