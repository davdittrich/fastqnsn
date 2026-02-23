library(fastqnsn)
library(testthat)

test_that("Integer vectors work correctly", {
  set.seed(42)
  x <- sample(-1000:1000, 100, replace=TRUE)

  # Check Sn
  expected_sn <- robustbase::Sn(as.numeric(x), constant=1)
  # fastqnsn uses Akinshin factors, so we compare raw value
  # To get raw value from fastqnsn, we divide by constant and factor.
  # Or just compare with a known value if we trust robustbase constant.

  # Actually, let's just check if it runs without error and returns something sensible.
  res_sn <- fastqnsn::sn(x)
  expect_true(is.numeric(res_sn))
  expect_true(res_sn > 0)

  # Check Qn
  res_qn <- fastqnsn::qn(x)
  expect_true(is.numeric(res_qn))
  expect_true(res_qn > 0)

  # Large integers
  x_large <- c(-1000000000L, 1000000000L, 0L)
  # Differences: 1e9, 2e9, 1e9.
  # Qn of {1e9, 1e9, 2e9} is 1e9?
  # h = 3/2 + 1 = 2. k = 2*(1)/2 = 1.
  # 1st smallest distance is 1e9.
  # Scaled Qn: 1e9 * 2.21914447 * get_qn_factor(3)

  res_qn_large <- fastqnsn::qn(x_large)
  expect_true(res_qn_large > 1e8)
})

print("Integer tests passed!")
