library(testthat)
library(Ball)
context("bd.gwas.test function")
skip_on_cran()

test_that("Error if computation result for ball divergence is wrong! (two-sample)", {
  set.seed(1)
  num <- 100
  snp_num <- 5
  x <- rnorm(num)
  snp <- sapply(1:snp_num, function(i) {
    sample(0:1, size = num, replace = TRUE)
  })
  res <- bd.gwas.test(x = x, snp = snp, num.threads = 1, num.permutations = 0, verbose = FALSE)
  bd_gwas_stats <- res[["statistic"]]
  for(i in 1:snp_num) {
    label <- snp[, i]
    new_x <- x[order(label)]
    size <- as.vector(table(label))
    bd_value <- prod(size) * bd(x = new_x, size = size) / sum(size)
    expect_equal(bd_gwas_stats[i], as.double(bd_value))
  }
})

test_that("Error if computation result for ball divergence is wrong! (K-sample)", {
  set.seed(1)
  num <- 100
  snp_num <- 5
  x <- as.matrix(rnorm(num))
  snp <- sapply(1:snp_num, function(i) {
    sample(0:2, size = num, replace = TRUE)
  })
  res <- bd.gwas.test(x = x, snp = snp, num.threads = 1, num.permutations = 0, verbose = FALSE)
  bd_gwas_stats <- res[["statistic"]]
  for(i in 1:snp_num) {
    label <- snp[, i]
    ulabel <- sort(unique(label))
    bd_value <- 0
    for (ulabel1 in ulabel) {
      for (ulabel2 in setdiff(ulabel, ulabel1)) {
        x1 <- x[label == ulabel1, , drop = FALSE]
        x2 <- x[label == ulabel2, , drop = FALSE]
        size <- c(nrow(x1), nrow(x2))
        bd_value <- bd_value + (prod(size) * bd(x1, x2) / sum(size))
      }
    }
    bd_value <- bd_value / 2
    expect_equal(bd_gwas_stats[i], as.double(bd_value))
  }
})

test_that("Error if computation result for ball divergence is wrong when multi-thread computation!", {
  set.seed(1)
  num <- 100
  snp_num <- 5
  x <- rnorm(num)
  snp <- sapply(1:snp_num, function(i) {
    sample(0:2, size = num, replace = TRUE)
  })
  res <- bd.gwas.test(x = x, snp = snp, num.threads = 1, num.permutations = 29999, verbose = FALSE)
  res1 <- bd.gwas.test(x = x, snp = snp, num.threads = 2, num.permutations = 29999, verbose = FALSE)
  expect_equal(res[["statistic"]], res1[["statistic"]])
  expect_equal(res[["permuted_statistic"]], res1[["permuted_statistic"]])
  expect_equal(res[["p.value"]], res1[["p.value"]])
})
