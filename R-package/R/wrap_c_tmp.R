#' @title compute joint kernel matrix
#' @noRd
#' @examples 
#' library(Ball)
#' set.seed(1)
#' num <- 100
#' x <- rnorm(num)
#' dst_x <- as.vector(dist(x))
#' res1 <- kbcov_cross_kernel_wrap_c(dst_x, k = 1, num = num)
#' res2 <- kbcov_margin_kernel_wrap_c(dst_x, num = num)
#' all(res1 == res2)
#' y <- rnorm(num)
#' dst_y <- as.vector(dist(y))
#' res <- kbcov_cross_kernel_wrap_c(c(dst_x, dst_y), k = 2, num = num)
#' dim(res)
#' term1 <- kbcov_joint_kernel_wrap_c(c(dst_x, dst_y), k = 2, num = num)
#' term2_1 <- kbcov_margin_kernel_wrap_c(dst_x, num = num)
#' term2_2 <- kbcov_margin_kernel_wrap_c(dst_y, num = num)
#' term3 <- kbcov_cross_kernel_wrap_c(c(dst_x, dst_y), k = 2, num = num)
#' bcov_v2 <- mean(term1) + mean(term2_1 * term2_2) - 2 * mean(term3)
#' bcov_v1 <- bcov(x, y)
#' bcov_v1 == bcov_v2
kbcov_cross_kernel_wrap_c <- function(x, 
                                      k, 
                                      num, 
                                      distance = TRUE, 
                                      num.threads = 1, 
                                      weight = "constant") 
{
  K <- as.integer(k)
  N <- as.integer(num)
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  weight_type <- as.integer(which(WEIGHT_TYPE == weight))
  
  xy <- as.double(x)
  kern_mat <- double(num^(k + 1))
  res <- .C("cross_kernel_matrix_bias_crude", kern_mat, xy, K, N, 
            num.threads, weight_type)
  rm(kern_mat); gc(reset = TRUE, verbose = FALSE)
  kern_mat <- array(data = 0, dim = rep(num, k + 1))
  kern_mat[1:(num^(k+1))] <- res[[1]]
  kern_mat
}

#' @title compute joint kernel matrix
#' @noRd
#' @examples 
#' rm(list = ls()); gc(reset = TRUE)
#' library(Ball)
#' set.seed(1)
#' num <- 100
#' x <- rnorm(num)
#' dst_x <- as.vector(dist(x))
#' res1 <- kbcov_joint_kernel_wrap_c(dst_x, k = 1, num = num)
#' res2 <- kbcov_margin_kernel_wrap_c(dst_x, num = num)
#' all(res1 == res2)
#' y <- rnorm(num)
#' dst_y <- as.vector(dist(y))
#' res <- kbcov_joint_kernel_wrap_c(c(dst_x, dst_y), k = 2, num = num)
#' dim(res)
#' dst_z <- as.vector(dist(rnorm(num)))
#' res <- kbcov_joint_kernel_wrap_c(c(dst_x, dst_y, dst_z), k = 3, num = num)
#' mean(res)
kbcov_joint_kernel_wrap_c <- function(x, 
                                      k, 
                                      num, 
                                      distance = TRUE, 
                                      num.threads = 1, 
                                      weight = "constant") 
{
  K <- as.integer(k)
  N <- as.integer(num)
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  weight_type <- as.integer(which(WEIGHT_TYPE == weight))
  
  xy <- as.double(x)
  kern_mat <- double((num + 1) * num / 2)
  res <- .C("joint_kernel_matrix_bias", kern_mat, xy, K, N, 
            num.threads, weight_type)
  rm(kern_mat); gc(reset = TRUE, verbose = FALSE)
  kern_mat <- matrix(0, nrow = num, ncol = num)
  kern_mat[lower.tri(kern_mat, diag = TRUE)] <- res[[1]]
  kern_mat <- kern_mat + t(kern_mat)
  diag(kern_mat) <- diag(kern_mat) / 2
  kern_mat
}

#' @title compute marginal kernel matrix
#' @noRd
kbcov_margin_kernel_wrap_c <- function(x, num, 
                                       distance = TRUE, 
                                       num.threads = 1, 
                                       weight = "constant") 
{
  N <- as.integer(num)
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  weight_type <- as.integer(which(WEIGHT_TYPE == weight))
  
  xy <- as.double(x)
  kern_mat <- double((num + 1) * num / 2)
  res <- .C("bdd_matrix_bias", kern_mat, xy, N, num.threads, weight_type)
  rm(kern_mat); gc(reset = TRUE, verbose = FALSE)
  kern_mat <- matrix(0, nrow = num, ncol = num)
  kern_mat[lower.tri(kern_mat, diag = TRUE)] <- res[[1]]
  kern_mat <- kern_mat + t(kern_mat)
  diag(kern_mat) <- diag(kern_mat) / 2
  kern_mat
}