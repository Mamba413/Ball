library(Ball)
library(testthat)
require(mvtnorm)
context("bcorsis function")
skip_on_cran()

# test_that("Ball Correlation works when missing values exist!", {
#   set.seed(1)
#   n <- 150
#   p <- 200
#   x <- matrix(rnorm(n * p), nrow = n)
#   eps <- rnorm(n)
#   y <- 3 * x[, 1] + 5 * (x[, 3])^2 + eps
#   
#   x[1:3, 1] <- NA
#   x[1:3, 3] <- NA
#   res1 <- bcorsis(y = dist(y), x = x, distance = TRUE)
#   test_res1_1 <- res1[["complete.info"]][["statistic"]][c(1, 3), ]
#   test_res1_2 <- res1[["complete.info"]][["statistic"]][2, , drop = TRUE]
#   
#   res2 <- bcorsis(y = dist(y[-(1:3)]), x = x[-(1:3), ], distance = TRUE)
#   test_res2_1 <- res2[["complete.info"]][["statistic"]][c(1, 3), ]
#   test_res2_2 <- res2[["complete.info"]][["statistic"]][2, , drop = TRUE]
#   
#   expect_true(all(test_res1_1 == test_res2_1))
#   expect_true(all(test_res1_2 != test_res2_2))
  
  ### Be cautious for the distance matrix:
  # load("tests/testthat/test.RData")
  # noise_dist <- matrix(rnorm(n = prod(dim(y_dist))), nrow = nrow(y_dist))
  # noise_dist <- (noise_dist + t(noise_dist)) / 2
  # diag(noise_dist) <- 0
  # min_y_dist <- min(y_dist[upper.tri(y_dist)])
  # y_dist <- y_dist + 0.01 * min_y_dist * noise_dist
  # res <- bcorsis(x = geno.use[, 1:2], y = y_dist, distance = TRUE)
# })

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

test_that("(Univariate) Incoherence between bcov and bcor!", {
  set.seed(1)
  num <- 10
  x <- rnorm(num)
  y <- rnorm(num)
  
  bcor_value <- bcor(x, y)
  bcov_based_value <- bcov(x, y) / sqrt(bcov(x, x) * bcov(y, y))
  
  names(bcor_value) <- NULL
  names(bcov_based_value) <- NULL
  expect_equal(bcor_value, bcov_based_value)
})


test_that("(Multivariate) Incoherence between bcov and bcor!", {
  set.seed(1)
  num <- 10
  x <- as.matrix(dist(rnorm(num)))
  y <- as.matrix(dist(rnorm(num)))
  
  bcor_value <- bcor(x, y, distance = TRUE)
  bcov_based_value <- bcov(x, y, distance = TRUE) / sqrt(bcov(x, x, distance = TRUE) * bcov(y, y, distance = TRUE))
  
  names(bcor_value) <- NULL
  names(bcov_based_value) <- NULL
  expect_equal(bcor_value, bcov_based_value)
})


test_that("(Multivariate) Incoherence between bcov and bcor!", {
  set.seed(1)
  num <- 10
  x <- matrix(rnorm(num * 2), nrow = num)
  y <- matrix(rnorm(num * 2), nrow = num)
  
  bcor_value <- bcor(x, y)
  bcov_based_value <- bcov(x, y) / sqrt(bcov(x, x) * bcov(y, y))
  
  names(bcor_value) <- NULL
  names(bcov_based_value) <- NULL
  expect_equal(bcor_value, bcov_based_value)
})


