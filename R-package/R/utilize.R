.onUnload <- function (libpath) {
  library.dynam.unload("Ball", libpath)
}

WEIGHT_TYPE <- c("constant", "probability", "chisquare", "rbf")
BD_WEIGHT_TYPE <- c("constant", "variance", "rbf")

BCOR_WEIGHT_STATS <- c("bcor.constant", "bcor.probability", "bcor.chisquare")
BCOV_WEIGHT_STATS <- c("bcov.constant", "bcov.probability", 
                       "bcov.chisquare", "bcov.rbf")

BD_WEIGHT_STATS <- c("bd.constant", "bd.variance", "bd.rbf")
KBD_WEIGHT_STATS <- c("kbd.sum.constant", "kbd.sum.variance", 
                      "kbd.max.constant", "kbd.max.variance", 
                      "kbd.maxsum.constant", "kbd.maxsum.variance")

center_bdd_matrix <- function(bdd) {
  bdd <- sweep(bdd, 2, colMeans(bdd)) - rowMeans(bdd) + mean(bdd)
  bdd
}

#' Hall-Buckley-Eagleson method
#'
#' Computes the cdf of a positively-weighted sum of chi-squared random variables with the Hall-Buckley-Eagleson (HBE) method.
#' @keywords distribution
#' @references
#' \itemize{
#'   \item P. Hall. Chi squared approximations to the distribution of a sum of independent random variables. \emph{The Annals of Probability}, 11(4):1028-1036, 1983.
#'   \item M. J. Buckley and G. K. Eagleson. An approximation to the distribution of quadratic forms in normal random variables. \emph{Australian Journal of Statistics}, 30(1):150-159, 1988.
#' }
#' @examples
#' hbe(c(1.5, 1.5, 0.5, 0.5), 10.203)            # should give value close to 0.95
#' @noRd
hbe <- function(coeff, x){
  # compute cumulants and nu
  K_1 <- sum(coeff)
  K_2 <- 2 * sum(coeff^2)
  K_3 <- 8 * sum(coeff^3)
  nu <- 8 * (K_2^3) / (K_3^2)
  
  # gamma parameters for chi-square
  gamma_k <- nu/2
  gamma_theta <- 2
  
  # need to transform the actual x value to x_chisqnu ~ chi^2(nu)
  # This transformation is used to match the first three moments
  # First x is normalised and then scaled to be x_chisqnu
  x_chisqnu_vec <- sqrt(2 * nu / K_2) * (x - K_1) + nu
  
  # now this is a chi_sq(nu) variable
  p_chisqnu_vec <- pgamma(x_chisqnu_vec, shape = gamma_k, scale = gamma_theta)
  p_chisqnu_vec
}

#' calculate Pvalue
#'
#' @param statValue Statistic Value
#' @param NullDistribution Ball statistic distribution when null hypothesis is true
#' @noRd
#' @return p-value
#'
calculatePvalue <- function(statValue, NullDistribution) {
  surpass_number <- sum(statValue < NullDistribution)
  if (surpass_number == 0) {
    p.value <- (surpass_number + 1) / (length(NullDistribution) + 1)
  } else {
    p.value <- surpass_number / length(NullDistribution)
  }
  p.value
}


#' Estimate memory consumpation
#'
#' @param n Sample size
#' @param funs Functions
#' @noRd
#' @return the available flag
#' 
memoryAvailable <- function(n, funs) {
  sysname <- Sys.info()[1]
  # memory check only available for window platform
  if(sysname == "Windows") {
    MemoryAvailable <- TRUE
    sizeLevel <- (n/1000)^2
    if( funs == 'UBI.test' ) {
      queryMemory <- 0.05 * sizeLevel
    } else if( funs == 'BI.test' ) {
      queryMemory <- 0.16 * sizeLevel
    } else if( funs == 'UBD.test' ) {
      queryMemory <- 0.02 * sizeLevel
    } else if( funs == 'BD.test' ) {
      queryMemory <- 0.08 * sizeLevel
    }
    sysMemory <- memory.limit()/1024
    if(queryMemory > sysMemory) {
      # MemoryAvailable <- FALSE
      stop("Sample size too large and system memory is not available!")
    } 
  } else {
    if(n > 8000) {
      warning("You may suffer from memory insufficient!")
    }
  }
}

#' Examine x, y arguments in bcov.test, bcov
#' @inheritParams bcov.test
#' @noRd
#' 
examine_x_y_bcor <- function(x, y) {
  dim_x <- dim(x)
  dim_y <- dim(y)
  if(is.null(dim_x) | is.null(dim_y)) {
    stop("x or y is NULL!")
  }
  n <- dim_x[1]
  if(n != dim_y[1]) {
    stop("x and y have different sample sizes!")
  }
  if(any(apply(y, 2, anyNA))) {
    stop("Missing data in y!")
  }
  if((dim_x[2] == 1) & (dim_y[2] == 1)) {
    p <- 1
  } else {
    p <- -1
  }
  c(n, p)
}


#' Examine x, y arguments in bcov.test, bcov
#' @inheritParams bcov.test
#' @noRd
#' 
examine_x_y <- function(x, y) {
  dim_x <- dim(x)
  dim_y <- dim(y)
  if(is.null(dim_x) | is.null(dim_y)) {
    stop("x or y is NULL!")
  }
  n <- dim_x[1]
  if(n != dim_y[1]) {
    stop("x and y have different sample sizes!")
  }
  if(any(apply(y, 2, anyNA))) {
    stop("Missing data in y!")
  }
  if(any(apply(x, 2, anyNA))) {
    stop("Missing data in x!")
  }
  if((dim_x[2] == 1) & (dim_y[2] == 1)) {
    p <- 1
  } else {
    p <- -1
  }
  c(n, p)
}

#' Examine x, y arguments in bcov.test, bcov
#' @inheritParams bcov.test
#' @noRd
#' 
examine_x_y_bcov <- function(x, y) {
  if (anyNA(x) || anyNA(y)) {
    stop("Missing value exist!")
  }
  if (length(x) != length(y)) {
    stop("x and y have different sample sizes!")
  }
}


#' Examine seed arguments in bcov.test, bd.test
#' @param seed A integer number
#' @return A integer number
#' @noRd
#' 
examine_seed_arguments <- function(seed) {
  if(is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  seed
}


#' Examine weight arguments in bcov.test, bd.test
#' @param weight A bool or character value
#' @param fun "bd.test" or "bcov.test"
#' @return A integer number
#' @noRd
#'
examine_weight_arguments <- function(weight) {
  if (is.logical(weight) || is.character(weight)) {
    if (is.logical(weight)) {
      weight <- ifelse(weight, "probability", "constant")
    } else {
      weight <- match.arg(arg = weight, choices = WEIGHT_TYPE)
    }
    return(weight)
  } else {
    stop("The weight arguments is invalid!")
  }
}

#' Title
#'
#' @param category 
#' @param p 
#'
#' @noRd
#' 
examine_category <- function(category, p) {
  if (is.logical(category)) {
    if (category) {
      category_index <- 1:p
    } else {
      category_index <- c()    
    }
  } else {
    stopifnot(all(category < 0) || all(category > 0))
    if (any(category > 0)) {
      category_index <- category
    } else {
      category_index <- setdiff(1:p, category)
    }
  }
  category_index
}



select_ball_stat <- function(ball_stat, weight, type = "bcov", fun_name = "bcov") {
  if (fun_name == "bcorsis") 
  {
    if (weight == "constant") {
      ball_stat <- ball_stat[, 1]
    } else if (weight == "probability") {
      ball_stat <- ball_stat[, 2]
    } else if (weight == "chisquare") {
      ball_stat <- ball_stat[, 3]
    }
  } else {
    if (weight == "constant")
    {
      ball_stat <- ball_stat[1]
      names(ball_stat) <- ifelse(type == "bcov", BCOV_WEIGHT_STATS[1], BCOR_WEIGHT_STATS[1])
    } else if (weight == "probability") {
      ball_stat <- ball_stat[2]
      names(ball_stat) <- ifelse(type == "bcov", BCOV_WEIGHT_STATS[2], BCOR_WEIGHT_STATS[2])
    } else if (weight == "chisquare") {
      ball_stat <- ball_stat[3]
      names(ball_stat) <- ifelse(type == "bcov", BCOV_WEIGHT_STATS[3], BCOR_WEIGHT_STATS[3])
    }
  }
  return(ball_stat)
}


#' Examine size arguments in bcov.test, bd.test
#' @param size A integer vector
#' @noRd
#' 
examine_size_arguments <- function(size) {
  # self examine:
  if(is.null(size)) {
    stop("size arguments is needed")
  }
  size <- as.integer(size)
  if(any(is.na(size)) | any(size <= 0) | (length(size)==1)) {
    stop("size arguments is invalid!")
  }
}


#' Examine R arguments in bcov.test, bd.test
#' @param R A integer number
#' @noRd
#' 
examine_R_arguments <- function(R) {
  if(is.null(R) | (R < 0)) {
    stop("R arguments is invalid!")
  }
}


#' Examine num.threads arguments in bcov.test, bd.test
#' @param R A integer number
#' @noRd
#' 
examine_threads_arguments <- function(num.threads) {
  if(is.null(num.threads) | (num.threads < 1)) {
    num.threads <<- 0;
    # stop("num.threads arguments is invalid!")
  }
}


#' Examine type arguments in bcov.test, bd.test
#' @param type "bcor" or "bcov"
#' @noRd
#' 
examine_type_arguments <- function(type) {
  if(all(!(type %in% c("bcov", "bcor")))) {
    type <- "bcov"
  }
  type
}


#' Examine dimension equality of input arguments x and y in bd.test
#'
#' @param x numeric matrix
#' @param y numeric matrix
#' @return dimension
#' @noRd
#' 
examine_dimension <- function(x, y) {
  dim_x <- dim(x)
  dim_y <- dim(y)
  p1 <- dim_x[2] 
  p2 <- dim_y[2]
  if(p1 != p2) {
    stop("x and y with different dimension!")
  }
  p1
}



#' get vectorized distance matrix
#' @inheritParams bd.test
#' @return vectorized distance matrix
#' @noRd
get_vectorized_distance_matrix <- function(x, y) {
  n1 <- dim(x)[1] 
  n2 <- dim(y)[1]
  n <- n1 + n2
  xy <- rbind(x, y)
  dxy <- as.vector(dist(xy))
  list(dxy, n1, n2)
}


#' @inheritParams bd.test
#' @return matrix
#' @noRd
#' 
get_matrixed_x <- function(x, y) {
  if(is.null(x)) {
    x <- y
  }
  as.matrix(x)
}



#' return candidate set size in bcorsis function
#' @inheritParams bcorsis
#' @param n sample size
#' @return size of candidate set
#' @noRd
#'
examine_candiate_size <- function(n, candidate, p) {
  if(p > n) {
    if(is.numeric(candidate)) {
      if(candidate <= 0) {
        stop("candidate argument is invalid!")
      }
      final_d <- as.integer(candidate)
    } else {
      if(candidate == "small"){
        final_d <- floor(n/log(n))
      } else if(candidate == "large") {
        final_d <- n - 1
      } else {
        stop("candidate argument is invalid!")
      }
    }
  } else {
    final_d <- p
    message("the number of covariate not larger than sample sizes, and SIS procedure is not essential")
  }
  final_d
}


get_screened_vars <- function(ids, rcory_result, final_d) {
  max_ids <- order(rcory_result, decreasing=T)
  chooseids <- max_ids[1:final_d]
  ids[chooseids]
}


preprocess_bcorsis_y <- function(y, y_p) {
  if(y_p != 1) {
    y_copy <- as.vector(dist(y))
    dst <- TRUE
  } else {
    y_copy <- as.vector(y)
    dst <- FALSE
  }
  list(y_copy, dst)
}


#' Examine method arguments in bcorsis
#' @inheritParams bcorsis
#' @noRd
#' 
examine_method_arguments <- function(method) {
  method <- head(unlist(strsplit(method, "-")), n = 1)
  if(!(method %in% c("standard", "pvalue", "interaction", "survival", "gam", "lm"))) {
    stop("method argument is invalid!")
  }
  method
}


examine_dst_method <- function(dst, method) {
  if(method %in% c("survival", "lm", "gam")) {
    if(dst) {
      messages <- " methods is not available when distance = TRUE"
      messages <- paste0(method, messages)
      stop(messages)
    }
  }
}