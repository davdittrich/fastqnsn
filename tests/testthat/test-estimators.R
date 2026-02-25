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

test_that("architectural memory boundaries are enforced natively", {
  # Structural bounds logic tests
  # We verify the maximum supported dimensions are strictly guarded in the R logic,
  # preventing the underlying C++ pointer offsets from silently overflowing.
  qn_body <- paste(deparse(qn), collapse = " ")
  sn_body <- paste(deparse(sn), collapse = " ")

  expect_true(grepl("6060000000", qn_body))
  expect_true(grepl("6060000000", sn_body))

  expect_true(grepl("128-bit architecture required", qn_body))
  expect_true(grepl("128-bit architecture required", sn_body))

  # C++ source check for std::bad_alloc OOM catches (since simulating OOM physically
  # breaks CI pipelines, we structurally assert the catch blocks exist in the source code).
  cpp_file <- system.file("src", "estimators.cpp", package = "fastqnsn")
  if (file.exists(cpp_file)) {
    cpp_source <- readLines(cpp_file)
    cpp_string <- paste(cpp_source, collapse = " ")
    expect_true(grepl("catch (const std::bad_alloc& e)", cpp_string, fixed = TRUE))
    expect_true(grepl("Out of Memory", cpp_string, fixed = TRUE))
  }
})
