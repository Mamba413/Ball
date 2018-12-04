library(testthat)
library(Ball)
context("bd.gwas.test function")
skip_on_cran()

test_that("The statistic value of bd.gwas.test is different to bd.test", {
  set.seed(1)
  num <- 50
  snp_num <- 100
  p <- 10
  x <- matrix(rnorm(num * p), nrow = num)
  snp <- sapply(1:snp_num, function(i) { 
    set.seed(i)
    sample(c(0, 1), size = num, replace = TRUE)
  })
  # gwas mode:
  res <- bd.gwas.test(x = x, snp = snp, seed = 4)
  # original mode:
  statistic_vec <- c()
  for (i in 1:snp_num) {
    statistic_vec[i] <- bd.test(x[order(snp[, i]), ], size = as.vector(table(snp[, i])), R = 0)
  }
  expect_equal(res[["statistic"]], expected = statistic_vec)
})

test_that("p-value is unreasonable", {
  set.seed(1)
  num <- 50
  snp_num <- 100
  p <- 10
  x1 <- matrix(rnorm(num * p), nrow = num/2)
  x2 <- matrix(rnorm(num * p, mean = 2), nrow = num/2)
  x <- rbind(x1, x2)
  snp <- sapply(1:snp_num, function(i) { 
    set.seed(i)
    sample(c(0, 1), size = num, replace = TRUE)
  })
  snp[, 1] <- rep(c(0, 1), each = num/2)
  res <- bd.gwas.test(x = x, snp = snp, seed = 2)
  # detect first SNP
  expect_true(res[["p.value"]][1] < 0.0005)
  # check with average p-value 
  mean_pvalue <- mean(res[["p.value"]][-1])
  expect_true(mean_pvalue > 0.47 && mean_pvalue < 0.53)
  # check with ks test
  suppressWarnings(unif_test_pvalue <- stats::ks.test(res[["p.value"]][-1], "punif", 0, 1)[["p.value"]])
  expect_true(unif_test_pvalue > 0.05)
})

test_that("Result is irreproducible", {
  set.seed(1)
  num <- 50
  snp_num <- 100
  p <- 10
  x <- matrix(rnorm(num * p), nrow = num)
  snp <- sapply(1:snp_num, function(i) { 
    set.seed(i)
    sample(c(0, 1), size = num, replace = TRUE)
  })
  res1 <- bd.gwas.test(x = x, snp = snp, seed = 4)
  res2 <- bd.gwas.test(x = x, snp = snp, seed = 4)
  expect_equal(res1, res2)
})

test_that("Result with different num.thread is irreproducible", {
  set.seed(1)
  num <- 50
  snp_num <- 100
  p <- 10
  x <- matrix(rnorm(num * p), nrow = num)
  snp <- sapply(1:snp_num, function(i) { 
    set.seed(i)
    sample(c(0, 1), size = num, replace = TRUE)
  })
  res1 <- bd.gwas.test(x = x, snp = snp, seed = 4)
  res2 <- bd.gwas.test(x = x, snp = snp, seed = 4, num.threads = 2)
  expect_equal(res1, res2)
})
