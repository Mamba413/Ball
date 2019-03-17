#' calculate Pvalue
#'
#' @param statValue Statistic Value
#' @param NullDistribution Ball statistic distribution when null hypothesis is true
#' @noRd
#' @return p-value
#'
calculatePvalue <- function(statValue, NullDistribution) {
  (sum(statValue < NullDistribution) + 1) / (length(NullDistribution) + 1)
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
  WEIGHT_METHODS <- c(TRUE, FALSE, "none", "chisq", "prob")
  method <- pmatch(weight, WEIGHT_METHODS)
  if (is.na(method)) 
    stop("invalid weight method")
  if (weight == TRUE) {
    weight <- "prob"
  } else if (weight == FALSE) {
    weight <- "none"
  }
  return(weight)
}


select_ball_stat <- function(ball_stat, weight, type = "bcov", fun_name = "bcov") {
  if (fun_name == "bcorsis") 
  {
    if (weight == "none") {
      ball_stat <- ball_stat[, 1]
    } else if (weight == "bcov") {
      ball_stat <- ball_stat[, 2]
    } else {
      ball_stat <- ball_stat[, 3]
    }
  } else {
    if (weight == "none")
    {
      ball_stat <- ball_stat[1]
      names(ball_stat) <- ifelse(type == "bcov", "bcov", "bcor")
    } else if (weight == "prob") {
      ball_stat <- ball_stat[2]
      names(ball_stat) <- ifelse(type == "bcov", "bcov.prob", "bcor.prob")
    } else {
      ball_stat <- ball_stat[3]
      names(ball_stat) <- ifelse(type == "bcov", "bcov.chisq", "bcor.chisq")
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
    stop("num.threads arguments is invalid!")
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
    y_copy <- as.vector(as.matrix(dist(y, diag = TRUE)))
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
      messages <- " methods is not available when dst = TRUE"
      messages <- paste0(method, messages)
      stop(messages)
    }
  }
}