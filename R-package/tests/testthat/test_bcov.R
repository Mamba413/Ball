library(Ball)
context("bd and bd.test function")
skip_on_cran()

test_that("Error if computation result for ball covariance is wrong!", {
  target_value <- 0.034214
  names(target_value) <- "bcov"
  expect_equal(bcov(1:10, 1:10), target_value)
  expect_equal(bcov.test(1:10, 1:10, R = 0), target_value)
})
