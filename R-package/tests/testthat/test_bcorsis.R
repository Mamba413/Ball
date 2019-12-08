library(Ball)
require(mvtnorm)
context("bcorsis function")
skip_on_cran()

test_that("Ball Correlation based Interaction Pursuit is unreasonable!", {
  set.seed(1)
  p <- 3000
  n <- 150
  x <- matrix(rnorm(n = n * p), nrow = n, ncol = p)
  y <- 4 * x[, 1] * x[, 2]
  
  res <- bcorsis(x = x, y = y, method = "interaction")
  expect_true(all(c(1, 2) %in% head(res[["ix"]], 5)))
})

test_that("Ball Correlation based survival screening is unreasonable!", {
  data("genlung")
  predictor <- genlung[["covariate"]]
  surv_status <- genlung[["survival"]]
  result <- bcorsis(x = predictor,
                    y = surv_status,
                    d = "small", method = "survival")
  top_gene <- colnames(predictor)[result[["ix"]]]
  expect_equal("hsa.miR.564", top_gene[1])
})

test_that("Ball Correlation based SIS is unreasonable!", {
  set.seed(1)
  n <- 150
  p <- 3000
  x <- matrix(rnorm(n * p), nrow = n)
  error <- rnorm(n)
  y <- 3*x[, 1] + 5*(x[, 3])^2 + error
  res <- bcorsis(y = y, x = x)
  expect_true(all(c(1, 3) %in% head(res[["ix"]])))
})

test_that("Ball Correlation based iterative SIS is unreasonable!", {
  set.seed(1)
  n <- 150
  p <- 1500
  sigma_mat <- matrix(0.5, nrow = p, ncol = p)
  diag(sigma_mat) <- 1
  x <- rmvnorm(n = n, sigma = sigma_mat)
  error <- rnorm(n)
  rm(sigma_mat); gc(reset = TRUE)
  y <- 6*(x[, 1])^2 + 5*(x[, 2])^2 + 5*x[, 8] - 8*x[, 16] + error
  index <- c(1, 2, 8, 16)
  res1 <- bcorsis(y = y, x = x, method = "lm", d = 15)
  expect_true(all(index %in% head(res1[["ix"]])))
  res2 <- bcorsis(y = y, x = x, method = "gam", d = 15)
  expect_true(all(index %in% head(res2[["ix"]])))
})

test_that("Ball Correlation based iterative SIS with probability weight is unreasonable!", {
  set.seed(1)
  n <- 150
  p <- 3000
  x <- matrix(rnorm(n * p), nrow = n)
  error <- rnorm(n)
  y <- 3*x[, 1] + 5*(x[, 3])^2 + error
  res <- bcorsis(y = y, x = x, weight = "prob")
  expect_true(all(c(1, 3) %in% head(res[["ix"]])))
  expect_equal("probability", res[["weight"]])
})

test_that("Ball Correlation based SIS with Chi-square weight is unreasonable!", {
  set.seed(1)
  n <- 150
  p <- 3000
  x <- matrix(rnorm(n * p), nrow = n)
  error <- rnorm(n)
  y <- 3*x[, 1] + 5*(x[, 3])^2 + error
  res <- bcorsis(y = y, x = x, weight = "chisq")
  expect_true(all(c(1, 3) %in% head(res[["ix"]])))
  expect_equal("chisquare", res[["weight"]])
})


test_that("Ball Correlation based SIS with categorical variables is unreasonable!", {
  set.seed(1)
  n <- 150
  p <- 3000
  x <- sapply(1:p, function(i) {
    sample(0:2, size = n, replace = TRUE)
  })
  eps <- rnorm(n)
  y <- 6 * x[, 1] - 7 * x[, 2] + 5 * x[, 3] + eps
  res <- bcorsis(x = x, y = y, category = TRUE)
  expect_true(all(c(1, 2, 3) %in% head(res[["ix"]])))
  
  x <- cbind(matrix(rnorm(n * 2), ncol = 2), x)
  res <- bcorsis(x = x, y = y, category = c(-1, -2))
  expect_true(all(c(3, 4, 5) %in% head(res[["ix"]])))
  
  x <- cbind(x[, 3:5], matrix(rnorm(n * p), ncol = p))
  res <- bcorsis(x = x, y = y, category = 1:3)
  expect_true(all(c(1, 2, 3) %in% head(res[["ix"]], n = 10)))
})
