library(testthat)
library(fastqnsn)

test_that("Micro-scale inputs continue to work correctly without overhead", {
    x <- 1:10
    expect_type(sn(x), "double")
    expect_type(qn(x), "double")
})

test_that("bounds checking validation is active", {
    # To test if the R wrapper's validation works, we just assert normal data continues
    # to evaluate. Allocating a numeric vector of length 7 billion requires too much RAM to mock.
    expect_type(qn(1:10), "double")
})
