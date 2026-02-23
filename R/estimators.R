#' @useDynLib fastqnsn, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @import RcppParallel
NULL

# Bias correction factors for Sn (n <= 100)
# Data from Akinshin (2022)
factors_sn <- c(
  NA, 0.7430, 1.8498, 0.9551, 1.3486, 0.9941, 1.1983, 1.0050, 1.1318, 1.0069,
  1.0959, 1.0063, 1.0742, 1.0051, 1.0601, 1.0038, 1.0501, 1.0028, 1.0430, 1.0022,
  1.0374, 1.0014, 1.0331, 1.0009, 1.0297, 1.0007, 1.0269, 1.0004, 1.0245, 1.0001,
  1.0226, 0.9999, 1.0209, 0.9997, 1.0195, 0.9998, 1.0183, 0.9996, 1.0172, 0.9997,
  1.0162, 0.9996, 1.0154, 0.9996, 1.0146, 0.9996, 1.0139, 0.9995, 1.0132, 0.9995,
  1.0126, 0.9995, 1.0123, 0.9995, 1.0117, 0.9995, 1.0113, 0.9996, 1.0109, 0.9996,
  1.0105, 0.9995, 1.0102, 0.9996, 1.0099, 0.9997, 1.0095, 0.9996, 1.0092, 0.9997,
  1.0090, 0.9997, 1.0088, 0.9996, 1.0085, 0.9997, 1.0084, 0.9997, 1.0081, 0.9997,
  1.0079, 0.9997, 1.0076, 0.9997, 1.0076, 0.9997, 1.0074, 0.9997, 1.0072, 0.9997,
  1.0070, 0.9997, 1.0069, 0.9997, 1.0067, 0.9998, 1.0066, 0.9997, 1.0065, 0.9998
)

# Bias correction factors for Qn (n <= 100)
factors_qn <- c(
  NA, 0.3995, 0.9939, 0.5133, 0.8441, 0.6122, 0.8589, 0.6700, 0.8736, 0.7201,
  0.8890, 0.7575, 0.9023, 0.7855, 0.9125, 0.8078, 0.9211, 0.8260, 0.9279, 0.8410,
  0.9338, 0.8537, 0.9389, 0.8644, 0.9430, 0.8737, 0.9468, 0.8819, 0.9501, 0.8890,
  0.9530, 0.8953, 0.9557, 0.9010, 0.9579, 0.9060, 0.9600, 0.9106, 0.9619, 0.9148,
  0.9636, 0.9185, 0.9652, 0.9220, 0.9667, 0.9252, 0.9680, 0.9281, 0.9692, 0.9309,
  0.9704, 0.9333, 0.9715, 0.9357, 0.9724, 0.9378, 0.9733, 0.9399, 0.9742, 0.9418,
  0.9750, 0.9435, 0.9757, 0.9453, 0.9765, 0.9469, 0.9771, 0.9484, 0.9777, 0.9498,
  0.9784, 0.9511, 0.9789, 0.9523, 0.9794, 0.9536, 0.9800, 0.9547, 0.9805, 0.9558,
  0.9809, 0.9568, 0.9814, 0.9578, 0.9818, 0.9587, 0.9822, 0.9597, 0.9826, 0.9605,
  0.9829, 0.9614, 0.9833, 0.9621, 0.9836, 0.9629, 0.9840, 0.9636, 0.9843, 0.9644
)

predict_sn_factor <- function(n) if (n %% 2 == 1) 1 + 0.707 / n - 7.181 / n^2 else 1 + 0.043 / n - 6.288 / n^2
predict_qn_factor <- function(n) if (n %% 2 == 1) 1 - 1.594 / n + 3.22 / n^2 else 1 - 3.672 / n + 11.087 / n^2

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
#'   \item For sample sizes n <= 10000, a high-speed deterministic C++ kernel is used.
#'   \item For larger samples (n > 10000), a multi-threaded algorithm via RcppParallel is used.
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
#'   \item For sample sizes n <= 3000, an optimized deterministic C++ kernel is used.
#'   \item For larger samples (n > 3000), a multi-threaded implementation using the Johnson-Mizoguchi selection algorithm.
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
  if (is.integer(x)) {
    return(C_qn_int_fast(x))
  }
  return(C_qn_fast(x))
}
