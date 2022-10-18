#' A New KBCOV Statistic
#'
#' @param X
#' @param K
#' @param n
#' @param alpha
#' @param method
#'
#' @return
#' @export
#'
#' @examples
#' library(Ball)
#' set.seed(1)
#' K=2
#' X <- list()
#' Y <- list()
#' n=1000
#' for(i in 1:K){
#'   Y[[i]] = rnorm(n)
#'   X[[i]] = as.vector(dist(Y[[i]]))
#' }
#' # when sample size is large the two statistics 
#' # have a similar statistic and test results.
#' KBCov(X, n)
#' bcov.test(Y[[1]], Y[[2]])
#' 
KBCov <- function(X,
                  n,
                  alpha = 0.05,
                  method = "limit") {
  K <- length(X)
  Kernel <- list()
  nthread = 1
  
  for (i in 1:K) {
    kernel_X = rep(0, n * (n + 1) / 2)
    Kernel[[i]] <- kbcov_margin_kernel_wrap_c(X[[i]], num = n)
  }
  
  estimation1 <- rep(0, K)
  estimation2 <- matrix(0, n, K)
  for (j in 1:K) {
    estimation1[j] <-
      1 / (n * (n - 1)) * sum(Kernel[[j]] - diag(diag(Kernel[[j]])))
    estimation2[, j] <- 1 / n * rowSums(Kernel[[j]])
  }
  
  estimation1_prod <- prod(estimation1)
  estimation2_prod <- apply(estimation2, 1, prod)
  term1 <- matrix(1, n, n)
  term2 <- (K - 1) ^ 2 * estimation1_prod
  term3 <- (K - 1) * estimation2_prod
  term5 <- matrix(0, n, n)
  term6 <- matrix(0, n, n)
  term8 <- matrix(0, n, n)
  term9 <- matrix(0, n, 1)
  for (j in 1:K) {
    term1 <- term1 * Kernel[[j]]
    term5 <- term5 + Kernel[[j]] * estimation1_prod / estimation1[j]
    term6 <-
      term6 + Kernel[[j]] * matrix(rep(estimation2_prod / estimation2[, j], n), n, n)
    term9 <-
      term9 + estimation2[, j, drop = FALSE] * estimation1_prod / estimation1[j]
    i <- j + 1
    while (i <= K) {
      term8 <-
        term8 + estimation2[, j, drop = FALSE] %*% t(estimation2[, i, drop = FALSE]) * estimation1_prod /
        (estimation1[j] *  estimation1[i]) + estimation2[, i, drop = FALSE] %*% t(estimation2[, j, drop = FALSE]) * estimation1_prod /
        (estimation1[j] * estimation1[i])
      i <- i + 1
    }
  }
  term3 <- matrix(rep(term3, n), n, n)
  term4 <- t(term3)
  term6 <- -term6
  term7 <- t(term6)
  term8 <- term8
  term9 <- matrix(rep((1 - K) * term9, n), n, n)
  term10 <- t(term9)
  H2 <-
    1 / (K * (2 * K - 1)) * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10)
  eigenvalues <-
    eigen(H2, only.values = TRUE, symmetric = TRUE)$values / n
  
  eigenvalues <- eigenvalues[eigenvalues >= 0]
  
  A2 <- 1
  B2 <- 1
  AB <- 1
  for (j in 1:K) {
    A2 <- A2 * Kernel[[j]]
    B2 <- B2 * sum(Kernel[[j]])
    AB <- AB * colSums(Kernel[[j]])
  }
  A2 <- sum(A2)
  
  AB <- sum(AB)
  kbcov = 1 / n ^ 2 * A2 + 1 / n ^ (2 * K) * B2 - 2 * AB / n ^ (K + 1)
  
  p_value = 1 - hbe(K * (2 * K - 1) * eigenvalues, n * kbcov)
  list(
    "statistic" = kbcov, 
    "p.value" = p_value
  )
}


#' A New KBCOV Statistic
#'
#' @param X
#' @param n
#' @param num.permutations 199
#' @param seed 1
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(1)
#' K=2
#' X <- list()
#' Y <- list()
#' n=1000
#' for(i in 1:K){
#'   Y[[i]] = rnorm(n)
#'   X[[i]] = as.vector(dist(Y[[i]]))
#' }
#' # when sample size is large the two statistics 
#' # have a similar statistic and test results.
#' KBCov_permuted(X, n)
#' KBCov(X, n)
#' bcov.test(Y[[1]], Y[[2]])
#' 
KBCov_permuted <- function(X,
                           n,
                           num.permutations = 299, 
                           seed = 1) 
{
  K <- length(X)
  Kernel <- list()
  nthread = 1
  
  for (i in 1:K) {
    kernel_X = rep(0, n * (n + 1) / 2)
    Kernel[[i]] <- kbcov_margin_kernel_wrap_c(X[[i]], num = n)
  }
  ma_value <- compute_MA_by_gram_matrix(Kernel)
  
  permuted_ma = numeric(num.permutations)
  for (r in 1:num.permutations) {
    set.seed(1 + r)
    for (i in 1:K) {
      index <- sample(1:n, size = n, replace = FALSE)
      Kernel[[i]] <- Kernel[[i]][index, index]
    }
    permuted_ma[r] <- compute_MA_by_gram_matrix(Kernel)
  }
  
  p_value = mean(c(ma_value <= permuted_ma, TRUE))
  list(
    "statistic" = ma_value, 
    "p.value" = p_value
  )
}

compute_MA_by_gram_matrix <- function(Kernel) {
  A2 <- 1
  B2 <- 1
  AB <- 1
  K <- length(Kernel)
  n <- nrow(Kernel[[1]])
  for (j in 1:K) {
    A2 <- A2 * Kernel[[j]]
    B2 <- B2 * sum(Kernel[[j]])
    AB <- AB * colSums(Kernel[[j]])
  }
  A2 <- sum(A2)
  AB <- sum(AB)
  ma_value <- 1 / n ^ 2 * A2 + 1 / n ^ (2 * K) * B2 - 2 * AB / n ^ (K + 1)
  ma_value
}


