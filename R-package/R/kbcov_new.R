#' Title
#'
#' @return
#'
#' @examples
#' K=2
#' n=200
#' x <- list(rnorm(n), rnorm(n))
#' bcov(x[[1]], x[[2]])
#' dst_x <- list()
#' for(i in 1:K){ dst_x[[i]] = as.vector(dist(x[[i]])) }
#' kbcov_R(dst_x)
kbcov_R <- function(X) {
  Kernel <- list()
  nthread = 1
  K <- length(X)
  for(i in 1:K){
    kernel_X = rep(0,n*(n+1)/2)
    A <- .C("bdd_matrix_bias", kernel_X = as.double(kernel_X), 
            X = as.double(X[[i]]), n = as.integer(n), 
            nthread = as.integer(nthread), 
            weight_type = as.integer(1))$kernel_X
    Kernel[[i]] = matrix(0,n,n)
    Kernel[[i]][upper.tri(Kernel[[i]],diag=TRUE)] = A
    Kernel[[i]][lower.tri(Kernel[[i]],diag=FALSE)] =
      Kernel[[i]][upper.tri(Kernel[[i]],diag=FALSE)]
    # Kernel[[i]] <- Kernel[[i]] + t(Kernel[[i]])
    # diag(Kernel[[i]]) <- diag(Kernel[[i]]) / 2
  }
  
  A2 <- 1
  B2 <- 1
  AB <- 1
  for (j in 1:K) {
    A2 <- A2*Kernel[[j]]
    B2 <- B2*sum(Kernel[[j]])
    AB <- AB*colSums(Kernel[[j]])
  }
  A2 <- sum(A2)
  AB <- sum(AB)
  kbcov = 1/n^2*A2+1/n^(2*K)*B2-2*AB/n^(K+1)
  kbcov
}