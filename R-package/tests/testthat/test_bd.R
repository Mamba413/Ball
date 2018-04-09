library(Ball)
context("bd and bd.test function")
skip_on_cran()

test_that("Error if computation result for ball divergence is wrong!", {
  target_value <- 2.4032
  names(target_value) <- "bd"
  expect_equal(bd(1:15, size = c(5, 5, 5)), target_value)
  expect_equal(bd.test(1:15, size = c(5, 5, 5), R = 0), target_value)
})


# test abnormal is not necessary. 
# If there exist any abnormal values, .C function will return detailed error message
test_that("Error if input data contain abnormal values", {
  x1 <- rnorm(20)
  x2 <- rnorm(20)
  x1[1] <- NA
  expect_error(bd.test(x1, x2))
  x1[1] <- Inf
  expect_error(bd.test(x1, x2))
})


test_that("Multi-thread support is valid!", {
  cat("Multi-thread computation via permutation.\n")
  x <- c(rnorm(100, mean = 0.1), rnorm(100, mean = 0))
  x1 <- matrix(rnorm(100 * 2), ncol = 2)
  x2 <- matrix(rnorm(100 * 2), ncol = 2)
  #
  res <- bd.test(x, size = c(100, 100), R = 1999, num.threads = 2)
  expect_length(res$permuted_stat, 1999)
  res <- bd.test(x1, x2, R = 1999, num.threads = 2)
  expect_length(res$permuted_stat, 1999)
  cat("Multi-thread computation via statistics.\n")
  x <- c(rnorm(500, mean = 0.1), rnorm(500, mean = 0))
  x1 <- matrix(rnorm(500 * 2), ncol = 2)
  x2 <- matrix(rnorm(500 * 2), ncol = 2)
  #
  res <- bd.test(x, size = c(500, 500), R = 199, num.threads = 2)
  expect_length(res$permuted_stat, 199)
  res <- bd.test(x1, x2, R = 199, num.threads = 2)
  expect_length(res$permuted_stat, 199)
})