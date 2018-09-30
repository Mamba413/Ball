library(Ball)
context("bd and bd.test function")
skip_on_cran()

test_that("Error if computation result for ball covariance is wrong!", {
  target_value <- 0.034214
  names(target_value) <- "bcov"
  expect_equal(bcov(1:10, 1:10), target_value)
  expect_equal(bcov.test(1:10, 1:10, R = 0), target_value)
})


test_that("Multi-thread support is valid!", {
  cat("Multi-thread computation via permutation for univariate test of independence problem.\n")
  Y <- rnorm(400)
  X <- rnorm(400)
  fit1 <- bcov.test(Y, X, R = 300, num.threads = 1)
  fit2 <- bcov.test(Y, X, R = 300, num.threads = 2)
  fit3 <- bcov.test(Y, X, R = 300, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
  
  cat("Multi-thread computation via permutation for multivariate test of independence problem.\n")
  Y <- matrix(rnorm(200*10), ncol = 10)
  X <- matrix(rnorm(200*10), ncol = 10)
  fit1 <- bcov.test(Y, X, R = 300, num.threads = 1)
  fit2 <- bcov.test(Y, X, R = 300, num.threads = 2)
  fit3 <- bcov.test(Y, X, R = 300, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})