# library(Ball)
# library(momentchi2)
# 
# K=3
# X <- list()
# Y <- list()
# n=800
# for(i in 1:K){
#   set.seed(i)
#   Y[[i]] = rnorm(n)
#   X[[i]] = as.vector(dist(Y[[i]]))
# }
# KBCov_old(X,K,n)
# bcov.test(Y[[1]],Y[[2]])
# 
# 
# KBCov_old <- function(X, K, n, method = "limit") {
#   Kernel <- list()
#   Y <- list()
#   
#   for(i in 1:K){
#     Kernel[[i]] <- kbcov_margin_kernel_wrap_c(X[[i]], num = n)
#     Y[[i]] <- matrix(0,n,n)
#     Y[[i]][lower.tri(Y[[i]],diag=FALSE)] = X[[i]]
#     Y[[i]] = t(Y[[i]])
#     Y[[i]][lower.tri(Y[[i]],diag=FALSE)] = X[[i]]
#   }
#   kbcov = bcov.test(Y, num.permutations = 0,distance = T)
# 
#   estimation1 <- rep(0,K)
#   estimation2 <- matrix(0,n,K)
#   for (j in 1:K) {
#     estimation1[j] <- 1/(n*n)*sum(Kernel[[j]])
#     estimation2[,j] <- 1/n*rowSums(Kernel[[j]])
#   }
#   
#   estimation1_prod <- prod(estimation1)
#   estimation2_prod <- apply(estimation2, 1, prod)
#   term1 <- matrix(1, n, n)
#   term2 <- (K - 1)^2 * estimation1_prod
#   term3 <- (K - 1) * estimation2_prod
#   term5 <- matrix(0, n, n)
#   term6 <- matrix(0, n, n)
#   term8 <- matrix(0, n, n)
#   term9 <- matrix(0, n, 1)
#   for (j in 1:K) {
#     term1 <- term1 * Kernel[[j]]
#     term5 <- term5 + Kernel[[j]] * estimation1_prod/estimation1[j]
#     term6 <- term6 + Kernel[[j]] * matrix(rep(estimation2_prod/estimation2[,j], n), n, n)
#     term9 <- term9 + estimation2[, j, drop = FALSE] * estimation1_prod/estimation1[j]
#      
#     while (i <= K) {
#       term8 <- term8 + estimation2[, j, drop = FALSE] %*% t(estimation2[ , i, drop = FALSE]) * estimation1_prod/(estimation1[j] *  estimation1[i]) + estimation2[, i, drop = FALSE] %*% t(estimation2[, j, drop = FALSE]) * estimation1_prod/(estimation1[j] * estimation1[i])
#       i <- i + 1
#     }
#   }
#   term3 <- matrix(rep(term3, n), n, n)
#   term4 <- t(term3)
#   term6 <- -term6
#   term7 <- t(term6)
#   term8 <- term8
#   term9 <- matrix(rep((1 - K) * term9, n), n, n)
#   term10 <- t(term9)
#   H2 <- 1/(K * (2 * K - 1)) * (term1 + term2 + term3 + term4 + term5 + term6 + term7 + term8 + term9 + term10)
#   eigenvalues <- eigen(H2, only.values = TRUE, symmetric = TRUE)$values/n
#   
# 
#  
#   p_value = 1 - hbe(K*(2*K-1)*eigenvalues[eigenvalues>0], n*kbcov)
#   return(c(kbcov,p_value))
# }
