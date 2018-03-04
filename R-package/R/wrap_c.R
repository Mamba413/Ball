#' @title compute Ball Divergence statistic
#' @description wrapper C function which compute Ball Divergence
#' @inheritParams bd.test
#' @param xy numeric vector
#' 
#' @return Ball Divergence statistic
#' @useDynLib Ball, .registration = TRUE
#' @noRd
bd_value_wrap_c <- function(xy, size, weight, dst, num.threads) {
  xy <- as.double(xy)
  bd <- as.double(numeric(1))
  weight <- as.integer(weight)
  dst <- as.integer(dst)
  num.threads <- as.integer(num.threads)
  #
  K <- as.integer(length(size))
  if(K == 2) {
    n1 <- as.integer(size[1])
    n2 <- as.integer(size[2])
    res <- .C("bd_stat", bd, xy, n1, n2, weight, dst, num.threads)
    N <- n1 + n2
  } else {
    size <- as.integer(size)
    N <- as.integer(sum(size))
    res <- .C("kbd_stat", bd, xy, size, N, K, weight, dst, num.threads)
  }
  bd <- res[[1]]
  names(bd) <- ifelse(weight, "wbd", "bd")
  list("statistic" = bd, "permuted_stat" = NULL, 
       "info" = list("N" = N, "K" = K, "size" = size, "weight" = as.logical(weight)))
}


#' @title compute Ball Divergence statistic and Ball Divergence statistic after permutation
#' @description wrapper C function which compute Ball Divergence and Ball Divergence statistic after permutation
#' @inheritParams bd.test
#' @param xy numeric vector
#' 
#' @return Ball Divergence statistic
#' @useDynLib Ball, .registration = TRUE
#' @noRd
bd_test_wrap_c <- function(xy, size, R, weight, dst, num.threads) {
  xy <- as.double(xy)
  dst <- as.integer(dst)
  R <- as.integer(R)
  weight <- as.integer(weight)
  num.threads <- as.integer(num.threads)
  #
  bd <- as.double(numeric(1))
  permuted_bd <- as.double(rep(-1, R))
  K <- as.integer(length(size))
  if(K == 2) {
    p <- as.integer(1)
    n1 <- as.integer(size[1])
    n2 <- as.integer(size[2])
    res <- .C("bd_test", bd, permuted_bd, xy, n1, n2, p, dst, R, weight, num.threads)
    
    N <- n1 + n2
  } else {
    size <- as.integer(size)
    N <- as.integer(sum(size))
    res <- .C("kbd_test", bd, permuted_bd, xy, size, N, K, dst, R, weight, num.threads)
  }
  bd <- res[[1]]
  names(bd) <- ifelse(weight, "wbd", "bd")
  permuted_bd <- res[[2]]
  permuted_bd <- permuted_bd[permuted_bd != (-1)]
  list('statistic' = bd, 'permuted_stat' = permuted_bd, 
       'info' = list('N' = N, 'K' = K, 'size' = size, 
                     'weight' = as.logical(weight), 'R' = R))
}



#' compute Ball Covariance statistic
#' @inheritParams bcov.test
#' @param x numeric vector.
#' @param y numeric vector.
#' @param n sample size. it must be integer value.
#' @param type if type == 1, function return bcov, otherwise, bcor insteaded.
#'
#' @return A list contain: Ball Covariance statistic, sample size, weight
#' @useDynLib Ball, .registration = TRUE
#' @noRd
#' 
bcov_value_wrap_c <- function(x, y, n, weight, dst, type) {
  bcov <- as.double(numeric(1))
  dst <- as.integer(dst)
  x <- as.double(x)
  y <- as.double(y)
  n <- as.integer(n)
  weight <- as.integer(weight)
  type <- ifelse(type == "bcov", 1, 2)
  type <- as.integer(type)
  res <- .C("bcov_stat", bcov, x, y, n, weight, dst, type)
  bcov <- res[[1]]
  names(bcov) <- ifelse(weight, "wbcov", "bcov")
  if(type == 2) {
    names(bcov) <- ifelse(weight, "wbcor", "bcor")
  }
  list('statistic' = bcov, "permuted_stat" = NULL,
       "info" = list("N" = res[[4]], "weight" = as.logical(res[[5]])))
}


#' compute Ball Covariance statistic and Ball Covariance statistic after permutation
#' @inheritParams bcov.test
#' @param x numeric vector.
#' @param y numeric vector.
#' @param n sample size. it must be integer value.
#'
#' @return A list contain: Ball Covariance statistic, Ball Covariance statistic after permutation, 
#' sample size, replication times, weight
#' @useDynLib Ball, .registration = TRUE
#' @noRd
#' 
bcov_test_wrap_c <- function(x, y, n, R, weight, dst, type) {
  dst <- as.integer(dst)
  x <- as.double(x)
  y <- as.double(y)
  n <- as.integer(n)
  weight <- as.integer(weight)
  R <- as.integer(R)
  type <- ifelse(type == "bcov", 1, 2)
  type <- as.integer(type)
  #
  bcov <- as.double(numeric(1))
  permuted_bcov <- as.double(numeric(R))
  res <- .C("bcov_test", bcov, permuted_bcov, x, y, n, R, weight, dst, type)
  bcov <- res[[1]]
  names(bcov) <- ifelse(weight, "wbcov", "bcov")
  list("statistic" = bcov, "permuted_stat" = res[[2]], 
       "info" = list("N" = res[[5]], "R" = res[[6]], "weight" = as.logical(res[[7]])))
}



#' compute Ball Covariance statistic for each variable
#' @inheritParams bcov.test
#' @param x numeric vector.
#' @param y numeric vector.
#' @param n sample size. it must be integer value.
#'
#' @return A list contain: Ball Covariance statistic, sample size, weight
#' @useDynLib Ball, .registration = TRUE
#' @noRd
#' 
apply_bcor_wrap <- function(x, y, n, R, weight, dst) {
  dst <- as.integer(dst)
  n <- as.integer(n)
  weight <- as.integer(weight)
  R <- as.integer(R)
  y <- as.double(y)
  type <- as.integer(2)
  #
  bcov_stat_list <- apply(x, 2, function(z) {
    # data prepare:
    if(dst) {
      z <- as.double(as.vector(as.matrix(dist(z, diag = TRUE))))
    } else {
      z <- as.double(z)
    }
    # calculate statistic or pvalue:
    bcov <- as.double(numeric(1))
    if(R == 0) {
      res <- .C("bcov_stat", bcov, z, y, n, weight, dst, type)
    } else {
      permuted_bcov <- as.double(numeric(R))
      res <- .C("bcov_test", bcov, permuted_bcov, z, y, n, R, weight, dst, type)
      res[[1]] <- calculatePvalue(res[[1]], res[[2]])
    }
    # return result:
    res[[1]]
  })
  if(R > 0) {
    bcov_stat_list <- 1 - bcov_stat_list
  }
  bcov_stat_list
}



#' @title Ball Correlation in survival.
#' @param x ordered covariate
#' @param t ordered survival event time
#' @param delta ordered survival event status
#' @param Sc Survfit object
#' @param n Sample size
#' @useDynLib Ball, .registration = TRUE
#' @noRd
#' 
#'
bcor_surv <- function(x, t, delta, Sc, n){
  # R function name should not the same as C(C++) function names
  RCT <- numeric(1)
  RCT <- .C("SRCT", as.double(t(x)), as.double(t(t)), as.double(t(delta)),
            as.double(t(Sc)), as.integer(n), RC = as.double(RCT))
  RCT[["RC"]]
}