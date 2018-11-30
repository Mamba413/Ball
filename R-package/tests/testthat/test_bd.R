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


test_that("Multi-thread support is valid!", {
  cat("Multi-thread computation via permutation for univariate K-sample problem.\n")
  n1 <- 200
  n2 <- 200
  n3 <- 200
  x <- rnorm(n1)
  y <- rnorm(n2)
  z <- rnorm(n3)
  fit1 <- bd.test(list(x, y, z), R = 399, num.threads = 1)
  fit2 <- bd.test(list(x, y, z), R = 399, num.threads = 2)
  fit3 <- bd.test(list(x, y, z), R = 399, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
  
  cat("Multi-thread computation via permutation for multivariate K-sample problem.\n")
  Y <- matrix(rnorm(100*10), ncol = 10)
  X <- matrix(rnorm(100*10), ncol = 10)
  Z <- matrix(rnorm(100*10), ncol = 10)
  
  fit1 <- bd.test(list(X, Y, Z), R = 399, size = c(100, 100, 100), num.threads = 1)
  fit2 <- bd.test(list(X, Y, Z), R = 399, size = c(100, 100, 100), num.threads = 2)
  fit3 <- bd.test(list(X, Y, Z), R = 399, size = c(100, 100, 100), num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})

test_that("Multi-thread computation via permutation for univariate two-sample problem", {
  n1 <- 200
  n2 <- 200
  x <- rnorm(n1)
  y <- rnorm(n2)
  fit1 <- bd.test(list(x, y), R = 399, num.threads = 1)
  fit2 <- bd.test(list(x, y), R = 399, num.threads = 2)
  fit3 <- bd.test(list(x, y), R = 399, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})

test_that("Multi-thread computation via permutation for multivariate two-sample problem", {
  Y <- matrix(rnorm(100*10), ncol = 10)
  X <- matrix(rnorm(100*10), ncol = 10)
  fit1 <- bd.test(list(X, Y), R = 399, size = c(100, 100, 100), num.threads = 1)
  fit2 <- bd.test(list(X, Y), R = 399, size = c(100, 100, 100), num.threads = 2)
  fit3 <- bd.test(list(X, Y), R = 399, size = c(100, 100, 100), num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})

test_that("output of formula interface is incorrect", {
  dat <- data.frame("x" = rnorm(100), "y" = as.factor(c(0, 1)))
  res1 <- bd.test(x ~ y, data = dat)
  expect_equal(strsplit(res1[["data.name"]], "\n")[[1]][1], "x by y")
})

test_that("Compare the outputs of formula interface and default interface (two-sample)", {
  dat <- data.frame("x" = rnorm(100), "y" = as.factor(c(0, 1)))
  res1 <- bd.test(x ~ y, data = dat)
  res2 <- bd.test(dat[["x"]][seq(1, 100, 2)], dat[["x"]][seq(2, 100, 2)])
  res1[["data.name"]] <- ""
  res2[["data.name"]] <- ""
  expect_equal(res1, res2)
})

test_that("Compare the outputs of formula interface and default interface (K-sample)", {
  res1 <- bd.test(Sepal.Width ~ Species, data = iris)
  iris_list <- split(iris[["Sepal.Width"]], iris[["Species"]])
  res2 <- bd.test(iris_list)
  res1[["data.name"]] <- ""
  res2[["data.name"]] <- ""
  expect_equal(res1, res2)
})
