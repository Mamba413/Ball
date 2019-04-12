#' @title compute Ball Divergence statistic
#' @description wrapper C function which compute Ball Divergence
#' @inheritParams bd.test
#' @param xy numeric vector
#' 
#' @return Ball Divergence statistic
#' @useDynLib Ball, .registration = TRUE
#' @noRd
bd_value_wrap_c <- function(xy, size, weight, distance, num.threads) {
  xy <- as.double(xy)
  bd <- as.double(numeric(1))
  weight <- as.integer(weight)
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  K <- as.integer(length(size))
  size <- as.integer(size)
  N <- as.integer(sum(size))
  res <- .C("bd_stat", bd, xy, size, N, K, weight, distance, num.threads)
  #
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
bd_test_wrap_c <- function(xy, size, num.permutations, weight, distance, num.threads) {
  xy <- as.double(xy)
  distance <- as.integer(distance)
  r <- as.integer(num.permutations)
  num.threads <- as.integer(num.threads)
  #
  K <- as.integer(length(size))
  stat_num <- ifelse(K == 2, 2, 6)
  bd <- as.double(numeric(stat_num))
  p_value <- as.double(numeric(stat_num))
  size <- as.integer(size)
  N <- as.integer(sum(size))
  res <- .C("bd_test", bd, p_value, xy, size, N, K, distance, r, num.threads)
  #
  stat_name_bd <- c("bd", "wbd")
  stat_name_kbd <- c("kbd.sum", "wkbd.sum", "kbd.max", "wkbd.max", "kbd.maxsum", "wkbd.maxsum")
  bd <- res[[1]]
  p_value <- res[[2]]
  if (K == 2) {
    names(bd) <- stat_name_bd
    bd <- bd[1]
    p_value <- p_value[1]
  } else {
    names(bd) <- stat_name_kbd
    bd <- bd[c(1, 3, 5)]
    p_value <- p_value[c(1, 3, 5)]
  }
  names(p_value) <- paste0(names(bd), ".pvalue")
  list('statistic' = bd, 'p.value' = p_value, 
       'info' = list('N' = N, 'K' = K, 'size' = size, 
                     'weight' = as.logical(weight), 'num.permutations' = num.permutations))
}


#' #' compute Ball Covariance statistic
#' #' @inheritParams bcov.test
#' #' @param x numeric vector.
#' #' @param y numeric vector.
#' #' @param n sample size. it must be integer value.
#' #' @param type if type == 1, function return bcov, otherwise, bcor insteaded.
#' #'
#' #' @return A list contain: Ball Covariance statistic, sample size, weight
#' #' @useDynLib Ball, .registration = TRUE
#' #' @noRd
#' #' 
#' bcov_value_wrap_c <- function(x, y, n, weight, distance, type, num.threads) {
#'   bcov <- as.double(numeric(1))
#'   distance <- as.integer(distance)
#'   x <- as.double(x)
#'   y <- as.double(y)
#'   n <- as.integer(n)
#'   weight_cp <- weight
#'   weight <- as.integer(weight)
#'   num.threads <- as.integer(num.threads)
#'   type <- ifelse(type == "bcov", 1, 2)
#'   type <- as.integer(type)
#'   if (type == 2) {
#'     res <- .C("bcov_stat", bcov, x, y, n, weight, distance, type, num.threads)
#'     bcov <- res[[1]]
#'   } else {
#'     bcov <- as.double(numeric(3))
#'     p_value <- as.double(numeric(3))
#'     r <- as.integer(numeric(1))
#'     weight <- as.integer(numeric(1))
#'     weight_cp <- examine_weight_arguments(weight_cp, "bcov.test")
#'     res <- .C("bcov_test", bcov, p_value, x, y, n, r, distance, num.threads)
#'     bcov <- res[[1]]
#'   }
#'   if (type == 1) {
#'     if (weight_cp == "none") {
#'       bcov <- bcov[1]
#'       names(bcov) <- "bcov"
#'     } else if(weight_cp == "prob") {
#'       bcov <- bcov[2]
#'       names(bcov) <- "bcov.prob"
#'     } else {
#'       bcov <- bcov[3]
#'       names(bcov) <- "bcov.chisq"
#'     } 
#'   } else {
#'     names(bcov) <- ifelse(weight, "wbcor", "bcor")
#'   }
#'   list('statistic' = bcov, "permuted_stat" = NULL,
#'        "info" = list("N" = res[[4]], "weight" = as.logical(weight)))
#' }


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
bcov_test_wrap_c <- function(x, y, n, num.permutations, distance, num.threads) {
  x <- as.double(x)
  y <- as.double(y)
  n <- as.integer(n)
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  r <- as.integer(num.permutations)
  #
  bcov <- as.double(numeric(3))
  p_value <- as.double(numeric(3))
  res <- .C("bcov_test", bcov, p_value, x, y, n, r, distance, num.threads)
  bcov <- res[[1]]
  p_value <- res[[2]]
  names(bcov) <- c("bcov", "bcov.prob", "bcov.chisq")
  names(p_value) <- paste0(names(bcov), ".pvalue")
  list('statistic' = bcov, 'p.value' = p_value,
       'info' = list("N" = res[[5]], "num.permutations" = res[[6]]))
}


#' compute K Ball Covariance statistic and Ball Covariance statistic after permutation
#' @inheritParams bcov.test
#' @param x numeric vector.
#' @param n sample size. it must be integer value.
#'
#' @return A list contain: Ball Covariance statistic, Ball Covariance statistic after permutation, 
#' sample size, replication times, weight
#' @useDynLib Ball, .registration = TRUE
#' @noRd
#' 
kbcov_test_wrap_c <- function(x, K, n, num.permutations, distance, num.threads) {
  x <- as.double(x)
  K <- as.integer(K)
  n <- as.integer(n)
  r <- as.integer(num.permutations)
  distance <- as.integer(distance) 
  num.threads <- as.integer(num.threads)
  #
  kbcov <- as.double(numeric(3))
  p_value <- as.double(numeric(3))
  res <- .C("kbcov_test", kbcov, p_value, x, K, n, r, distance, num.threads)
  bcov <- res[[1]]
  p_value <- res[[2]]
  names(bcov) <- c("bcov", "bcov.prob", "bcov.chisq")
  names(p_value) <- paste0(names(bcov), ".pvalue")
  list('statistic' = bcov, 'p.value' = p_value,
       'info' = list("N" = res[[5]], "num.permutations" = res[[6]]))
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
apply_bcor_wrap <- function(x, y, n, p, distance, weight, method, num.threads) {
  bcor_stat <- as.double(numeric(3*p))
  y <- as.double(y)
  x <- as.vector(x)
  x_number <- rep(1, p)
  f_number <- as.integer(p)
  size_num <- 2
  num <- as.integer(n)
  nthread <- integer(1)
  p <- as.integer(1)
  k <- as.integer(1)
  dst_y <- as.integer(distance)
  dst_x <- as.integer(0)
  nth <- as.integer(num.threads)
  #
  res <- .C("bcor_test", bcor_stat, y, x, x_number, f_number, size_num, num, p, k, dst_y, dst_x, nth)[[1]]
  bcor_stat <- matrix(res, ncol = 3, byrow = TRUE)
  colnames(bcor_stat) <- c("bcor", "bcor.prob", "bcor.chisq")
  screening_bcor_stat <- select_ball_stat(bcor_stat, weight = weight, fun_name = "bcorsis")
  if (method %in% c("interaction", "standard"))
  {
    return(list(bcor_stat, screening_bcor_stat))
  } else {
    return(screening_bcor_stat)
  }
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
bcor_surv <- function(x, time_value, delta, Sc, n){
  # R function name should not the same as C(C++) function names
  RCT <- numeric(1)
  RCT <- .C("SRCT_new", as.double(t(x)), as.integer(t(time_value)), as.integer(t(delta)),
            as.double(t(Sc)), as.integer(n), RC = as.double(RCT))
  RCT[["RC"]]
}