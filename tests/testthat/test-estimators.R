library(testthat)
library(fastqnsn)

test_that("sn handles basic numeric vectors correctly", {
  set.seed(42)
  x <- rnorm(100)
  res <- sn(x)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(res > 0)
})

test_that("qn handles basic numeric vectors correctly", {
  set.seed(42)
  x <- rnorm(100)
  res <- qn(x)
  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(res > 0)
})

test_that("na.rm works as expected", {
  x <- c(1, 2, 3, 4, 5, NA)
  expect_true(is.na(sn(x)))
  expect_true(is.na(qn(x)))
  
  expect_false(is.na(sn(x, na.rm = TRUE)))
  expect_false(is.na(qn(x, na.rm = TRUE)))
})

test_that("edge cases for small n are handled", {
  expect_true(is.na(sn(numeric(0))))
  expect_true(is.na(qn(numeric(0))))
  
  expect_true(is.na(sn(1)))
  expect_true(is.na(qn(1)))
  
  # n=2 is the minimum
  expect_false(is.na(sn(c(1, 2))))
  expect_false(is.na(qn(c(1, 2))))
})

test_that("outlier resistance is maintained", {
  set.seed(1)
  x <- rnorm(50)
  s1 <- sn(x)
  q1 <- qn(x)
  
  # Massive outlier
  x_out <- c(x, 1e10)
  expect_lt(abs(sn(x_out) - s1), 0.5) # Should be very close
  expect_lt(abs(qn(x_out) - q1), 0.5)
})

test_that("large scale data consistency", {
  # Triggering O(n log n) and parallel logic
  x <- rnorm(1500)
  expect_silent(sn(x))
  expect_silent(qn(x))
})
