#' @title compute Ball Divergence statistic
#' @description wrapper C function which compute Ball Divergence via limit distribution
#' @inheritParams bd.test
#' @noRd
bcov_limit_wrap_c <- function(x, y, num, distance, num.threads, weight) {
  if (!is.null(y)) {
    x <- list(x, y)
  }
  N <- as.integer(num)
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  weight_type <- as.integer(which(WEIGHT_TYPE == weight))
  
  bdd_xy_eigen <- matrix(data = 1, ncol = num, nrow = num)
  for (i in 1:length(x)) {
    xy <- as.double(x[[i]])
    bdd_xy <- double((num + 1) * num / 2)
    res <- .C("bdd_matrix_bias", bdd_xy, xy, N, num.threads, weight_type)
    rm(bdd_xy); gc(reset = TRUE, verbose = FALSE)
    bdd_xy <- matrix(0, nrow = num, ncol = num)
    bdd_xy[lower.tri(bdd_xy, diag = TRUE)] <- res[[1]]
    bdd_xy <- bdd_xy + t(bdd_xy)
    diag(bdd_xy) <- diag(bdd_xy) / 2
    # bdd_xy_eigen <- bdd_xy_eigen * bdd_xy
    bdd_xy_eigen <- bdd_xy_eigen * center_bdd_matrix(bdd_xy)
  }
  # bdd_xy_eigen <- center_bdd_matrix(bdd_xy_eigen)
  
  eigenvalue <- eigen(bdd_xy_eigen, only.values = TRUE, symmetric = TRUE)$values
  eigenvalue <- eigenvalue[eigenvalue > 0] / num
  # eigenvalue <- eigenvalue / 120
  eigenvalue
}

#' @title compute Ball Divergence statistic
#' @description wrapper C function which compute Ball Divergence via limit distribution
#' @inheritParams bd.test
#' @noRd
bd_limit_wrap_c <- function(xy, size, distance, num.threads) {
  xy <- as.double(xy)
  n1 <- as.integer(size[1])
  n2 <- as.integer(size[2])
  distance <- as.integer(distance)
  num.threads <- as.integer(num.threads)
  
  num <- as.integer(sum(size))
  bdd_xy <- double((num + 1) * num / 2)
  res <- .C("bdd_matrix_bias_two_group", bdd_xy, xy, n1, n2, num.threads)
  bdd_xy <- matrix(0, nrow = num, ncol = num)
  bdd_xy[lower.tri(bdd_xy, diag = TRUE)] <- res[[1]]
  bdd_xy <- bdd_xy + t(bdd_xy)
  diag(bdd_xy) <- diag(bdd_xy) / 2
  
  bdd_xy <- center_bdd_matrix(bdd_xy)
  eigenvalue <- eigen(bdd_xy, only.values = TRUE, symmetric = TRUE)$values
  eigenvalue <- eigenvalue[eigenvalue > 0] / num
  2 * eigenvalue
}

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
  stat_num <- ifelse(K == 2, 3, 6)
  bd <- as.double(numeric(stat_num))
  p_value <- as.double(numeric(stat_num))
  size <- as.integer(size)
  N <- as.integer(sum(size))
  res <- .C("bd_test", bd, p_value, xy, size, N, K, distance, r, num.threads)

  bd <- res[[1]]
  p_value <- res[[2]]
  if (K == 2) {
    names(bd) <- BD_WEIGHT_STATS
  } else {
    names(bd) <- KBD_WEIGHT_STATS
  }
  names(p_value) <- paste0(names(bd), ".pvalue")
  
  if (weight == BD_WEIGHT_TYPE[1]) {
    if (K == 2) {
      index <- 1
    } else {
      index <- c(1, 3, 5)
    }
  } else if (weight == BD_WEIGHT_TYPE[2]) {
    if (K == 2) {
      index <- 2
    } else {
      index <- c(2, 4, 6)
    }
  } else if (weight == BD_WEIGHT_TYPE[3]) {
    if (K == 2) {
      index <- 3
    } else {
      index <- c(2, 4, 6)
    }
  }
  return_bd <- bd[index]
  return_p_value <- p_value[index]
  
  list('statistic' = return_bd, 'p.value' = return_p_value, 
       'info' = list('statistic' = bd, "p.value" = p_value, 'N' = N, 'K' = K, 'size' = size, 
                     'weight' = weight, 'num.permutations' = num.permutations))
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
  bcov <- as.double(numeric(length(WEIGHT_TYPE)))
  p_value <- as.double(numeric(length(WEIGHT_TYPE)))
  res <- .C("bcov_test", bcov, p_value, x, y, n, r, distance, num.threads)
  bcov <- res[[1]]
  p_value <- res[[2]]
  names(bcov) <- BCOV_WEIGHT_STATS
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
  kbcov <- as.double(numeric(length(BCOV_WEIGHT_STATS)))
  p_value <- as.double(numeric(length(BCOV_WEIGHT_STATS)))
  res <- .C("kbcov_test", kbcov, p_value, x, K, n, r, distance, num.threads)
  bcov <- res[[1]]
  p_value <- res[[2]]
  names(bcov) <- BCOV_WEIGHT_STATS
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
apply_bcor_wrap <- function(x, y, n, p, distance, weight, method, num.threads, category) {
  p_all <- ncol(x)
  if (length(category) != 0) {
    x_category <- x[, category, drop = FALSE]
    x <- x[, -category, drop = FALSE]
  } else {
    x_category <- matrix(0, nrow = 0, ncol = 0)
  }
  p_continuous <- ncol(x)
  p_category <- ncol(x_category)
  
  if (distance) {
    p <- as.integer(0)
  } else {
    p <- as.integer(p)
  }
  y <- as.double(y)
  num <- as.integer(n)
  dst_x <- as.integer(0)
  nth <- as.integer(num.threads)
  dst_y <- as.integer(distance)
  if (p_continuous != 0) {
    bcor_stat1 <- as.double(numeric(3 * p_continuous))
    x <- as.double(as.vector(x))
    missing_flag <- as.integer(as.vector(!is.na(x)))
    x[missing_flag == 0] <- -9999.99
    x_number <- as.integer(rep(1, p_continuous))
    f_number <- as.integer(p_continuous)
    k <- as.integer(1)
    res <- .C("bcor_test", bcor_stat1, y, x, x_number, 
              f_number, num, p, k, dst_y, dst_x, nth, missing_flag)[[1]]
    bcor_stat1 <- matrix(res, ncol = 3, byrow = TRUE)
  }
  
  if (p_category != 0) {
    bcor_stat2 <- as.double(numeric(3 * p_category))
    x <- as.double(as.vector(x_category))
    x_number <- as.integer(rep(1, p_category))
    f_number <- as.integer(p_category)
    k <- as.integer(2)
    # not support for categorical data now: 
    missing_flag <- as.integer(rep(1, length(x)))
    #
    res <- .C("bcor_test", bcor_stat2, y, x, x_number, f_number, num, p, k, dst_y, dst_x, nth, missing_flag)[[1]]
    bcor_stat2 <- matrix(res, ncol = 3, byrow = TRUE)
    
    if (p_continuous == 0) {
      bcor_stat <- bcor_stat2
    } else {
      bcor_stat <- matrix(nrow = p_all, ncol = 3)
      bcor_stat[-category, ] <- bcor_stat1
      bcor_stat[category, ] <- bcor_stat2 
    }
  } else {
    bcor_stat <- bcor_stat1
  }
  
  colnames(bcor_stat) <- BCOR_WEIGHT_STATS
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