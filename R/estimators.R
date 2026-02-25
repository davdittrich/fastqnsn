#' @useDynLib fastqnsn, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppParallel
NULL

#' Robust Sn Scale Estimator
#'
#' Computes the robust Rousseeuw-Croux Sn scale estimator. It is a highly robust
#' alternative to the standard deviation and the Median Absolute Deviation (MAD).
#'
#' @param x A numeric vector.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#'
#' @details
#' The Sn estimator is defined as:
#' \deqn{S_n = c \cdot \text{lomed}_i \left\{ \text{himed}_j |x_i - x_j| \right\}}{Sn = c * lomed_i { himed_j |x_i - x_j| }}
#' This implementation uses a fast hybrid architecture:
#' \itemize{
#'   \item For sample sizes n <= 8192, a high-speed deterministic C++ kernel is used.
#'   \item For larger samples (n > 8192), a multi-threaded algorithm via RcppParallel is used.
#'   \item Optimized sorting: Uses Boost Spreadsort for small/medium samples and TBB parallel sort for large samples.
#' }
#' It incorporates modern consistency constants and finite-sample bias correction factors
#' as described by Akinshin (2022).
#'
#' @return A numeric scalar representing the estimated robust scale.
#'
#' @references
#' Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the Median Absolute Deviation. Journal of the American Statistical Association.
#' Akinshin, A. (2022). Finite-sample Rousseeuw-Croux scale estimators. arXiv:2209.12268.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(100)
#' # Typical usage
#' sn(x)
#'
#' # Resistant to extreme outliers
#' x_contaminated <- c(x, 1000, 2000, 3000)
#' sn(x_contaminated)
#'
#' @export
sn <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) {
    return(NA_real_)
  }
  if (n > 6060000000) {
    stop("fastqnsn Error: sample size n > 6.06 * 10^9 natively overflows 64-bit pair boundaries. 128-bit architecture required.")
  }
  if (is.integer(x)) {
    return(C_sn_int_fast(x))
  }
  return(C_sn_fast(x))
}

#' Robust Qn Scale Estimator
#'
#' Computes the robust Rousseeuw-Croux Qn scale estimator. It is based on the
#' first quartile of interpoint distances and is more efficient than MAD and Sn.
#'
#' @param x A numeric vector.
#' @param na.rm A logical value indicating whether NA values should be stripped before the computation proceeds.
#'
#' @details
#' The Qn estimator is defined as:
#' \deqn{Q_n = d \cdot \{ |x_i - x_j| ; i < j \}_{(k)}}{Qn = d * {|x_i - x_j| ; i < j}(k)}
#' This implementation uses a fast hybrid architecture:
#' \itemize{
#'   \item For sample sizes n <= 2048, an optimized deterministic C++ kernel is used.
#'   \item For larger samples (n > 2048), a multi-threaded implementation using the Johnson-Mizoguchi selection algorithm.
#'   \item Optimized sorting: Uses Boost Spreadsort for small/medium samples and TBB parallel sort for large samples.
#' }
#' It corrects the legacy typo in consistency constants (legacy 2.2219 vs modern 2.2191)
#' and applies finite-sample bias corrections from Akinshin (2022).
#'
#' @return A numeric scalar representing the estimated robust scale.
#'
#' @references
#' Rousseeuw, P. J., & Croux, C. (1993). Alternatives to the Median Absolute Deviation. Journal of the American Statistical Association.
#' Akinshin, A. (2022). Finite-sample Rousseeuw-Croux scale estimators. arXiv:2209.12268.
#' Johnson, D. B., & Mizoguchi, T. (1978). Selecting the Kth element in X + Y. SIAM Journal on Computing.
#'
#' @examples
#' set.seed(42)
#' x <- rnorm(100)
#' # Typical usage
#' qn(x)
#'
#' # Highly efficient on clean data
#' qn(x) / sd(x) # Should be close to 1
#'
#' @export
qn <- function(x, na.rm = FALSE) {
  if (na.rm) x <- x[!is.na(x)]
  n <- length(x)
  if (n < 2) {
    return(NA_real_)
  }
  if (n > 6060000000) {
    stop("fastqnsn Error: sample size n > 6.06 * 10^9 natively overflows 64-bit pair boundaries. 128-bit architecture required.")
  }
  if (is.integer(x)) {
    return(C_qn_int_fast(x))
  }
  return(C_qn_fast(x))
}
