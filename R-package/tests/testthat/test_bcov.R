library(Ball)
context("bd and bd.test function")
skip_on_cran()

test_that("Error if computation result for ball covariance is wrong!", {
  target_value <- 0.034214
  names(target_value) <- "bcov"
  expect_equal(bcov(1:10, 1:10), target_value)
  expect_equal(bcov.test(1:10, 1:10, R = 0), target_value)
})


test_that("Multi-thread computation via permutation for univariate test of independence problem", {
  Y <- rnorm(400)
  X <- rnorm(400)
  fit1 <- bcov.test(Y, X, R = 300, num.threads = 1)
  fit2 <- bcov.test(Y, X, R = 300, num.threads = 2)
  fit3 <- bcov.test(Y, X, R = 300, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})


test_that("Multi-thread computation via permutation for multivariate test of independence problem.", {
  Y <- matrix(rnorm(200*10), ncol = 10)
  X <- matrix(rnorm(200*10), ncol = 10)
  fit1 <- bcov.test(Y, X, R = 300, num.threads = 1)
  fit2 <- bcov.test(Y, X, R = 300, num.threads = 2)
  fit3 <- bcov.test(Y, X, R = 300, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})

test_that("Multi-thread computation via permutation for mutual independence test", {
  x <- lapply(rep(60, 3), rnorm)
  fit1 <- bcov.test(x, R = 299, num.threads = 1)
  fit2 <- bcov.test(x, R = 299, num.threads = 2)
  fit3 <- bcov.test(x, R = 299, num.threads = 4)
  expect_equal(fit1[["complete.info"]][["statistic"]], fit2[["complete.info"]][["statistic"]])
  expect_equal(fit1[["complete.info"]][["statistic"]], fit3[["complete.info"]][["statistic"]])
})

test_that("validity of formula interface", {
  res1 <- bcov.test(~ CONT + INTG, data = USJudgeRatings)
  expect_equal(strsplit(res1[["data.name"]], "\n")[[1]][1], "CONT and INTG")
  res1 <- bcov.test(~ CONT + INTG + DMNR, data = USJudgeRatings)
  expect_equal(strsplit(res1[["data.name"]], "\n")[[1]][1], "CONT and INTG and DMNR")
})

test_that("formula interface is identical to default interface (independence)", {
  res1 <- bcov.test(~ CONT + INTG, data = USJudgeRatings)
  res2 <- bcov.test(USJudgeRatings[["CONT"]], USJudgeRatings[["INTG"]])
  res1[["data.name"]] <- ""
  res2[["data.name"]] <- ""
  expect_equal(res1, res2)
})

test_that("formula interface is identical to default interface (mutual independence)", {
  res1 <- bcov.test(~ CONT + INTG + DMNR, data = USJudgeRatings)
  res2 <- bcov.test(list(USJudgeRatings[["CONT"]], 
                         USJudgeRatings[["INTG"]], 
                         USJudgeRatings[["DMNR"]]))
  res1[["data.name"]] <- ""
  res2[["data.name"]] <- ""
  expect_equal(res1, res2)
})

