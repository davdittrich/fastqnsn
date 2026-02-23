library(testthat)
library(fastqnsn)
library(robustbase)

test_that("Integer vectors work correctly", {
  set.seed(42)
  x_int <- as.integer(sample.int(1000, 100))

  expect_equal(sn(x_int), sn(as.numeric(x_int)))
  expect_equal(qn(x_int), qn(as.numeric(x_int)))
})

test_that("Small integer vectors match robustbase with corrections", {
  x <- c(1L, 5L, 2L, 8L, 3L)
  expect_true(is.finite(sn(x)))
  expect_true(is.finite(qn(x)))
})
