library(testthat)
library(Ball)
context("bd and bd.test function")
skip_on_cran()

test_that("Error if computation result for ball divergence is wrong!", {
  target_value <- 2.4032
  names(target_value) <- "kbd.sum"
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


test_that("bd, bd.test function return unmatched ball divergence statistic", {
  dat <- lapply(rep(50, 3), function(i) {
    rnorm(i)
  })
  res1 <- bd(dat, kbd.type = "sum")
  res2 <- bd(dat, kbd.type = "max")
  res3 <- bd(dat, kbd.type = "maxsum")
  expect_equal(names(res1), "kbd.sum")
  expect_equal(names(res2), "kbd.max")
  expect_equal(names(res3), "kbd.maxsum")
  
  res1 <- bd.test(dat, kbd.type = "sum")
  res2 <- bd.test(dat, kbd.type = "max")
  res3 <- bd.test(dat, kbd.type = "maxsum")
  expect_equal(names(res1[["statistic"]]), "kbd.sum")
  expect_equal(names(res2[["statistic"]]), "kbd.max")
  expect_equal(names(res3[["statistic"]]), "kbd.maxsum")
  
  expect_equal(names(res1[["p.value"]]), "kbd.sum.pvalue")
  expect_equal(names(res2[["p.value"]]), "kbd.max.pvalue")
  expect_equal(names(res3[["p.value"]]), "kbd.maxsum.pvalue")
})


# test_that("Multi-thread support is valid!", {
#   cat("Multi-thread computation via permutation.\n")
#   set.seed(4)
#   num <- 150
#   x <- c(rnorm(num, mean = 0.1), rnorm(num, mean = 0))
#   x1 <- matrix(rnorm(num * 2), ncol = 2)
#   x2 <- matrix(rnorm(num * 2), ncol = 2)
#   #
#   cat("Univariate case: \n")
#   t1 <- system.time(res1 <- bd.test(x, size = c(num, num), R = 1999, num.threads = 1))
#   t2 <- system.time(res2 <- bd.test(x, size = c(num, num), R = 1999, num.threads = 2))
#   expect_equal(res1$statistic, res2$statistic)
#   expect_equal(res1$p.value < 0.05, res2$p.value < 0.05)
#   expect_true(t1[3]/t2[3] > 1.1)
#   #
#   cat("Multivariate case: \n")
#   t1 <- system.time(res1 <- bd.test(x1, x2, R = 1999, num.threads = 1))
#   t2 <- system.time(res2 <- bd.test(x1, x2, R = 1999, num.threads = 2))
#   expect_equal(res1$statistic, res2$statistic)
#   expect_equal(res1$p.value < 0.05, res2$p.value < 0.05)
#   expect_true(t1[3]/t2[3] > 1.1)
# })