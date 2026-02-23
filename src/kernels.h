#ifndef KERNELS_H
#define KERNELS_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

uint64_t zig_qn_count(const double *x_ptr, size_t n, double trial);

void zig_sn_naive(const double *x_ptr, size_t n, double *out_ptr);

double zig_sn_select_k(const double *x_ptr, size_t n, size_t i, size_t k);
double zig_sn_deterministic(const double *x_ptr, size_t n, size_t i);

double zig_qn_jm_select(const double *x_ptr, size_t n, uint64_t k,
                        double *work_ptr, int32_t *iweight_ptr,
                        int32_t *left_ptr, int32_t *right_ptr);

double zig_whimed(double *a_ptr, int32_t *iw_ptr, size_t n, int64_t target);

#ifdef __cplusplus
}
#endif

#endif
