#' @title Sub-sampling Ball Divergence based Equality of Distributions Test
#' 
#' @description Performs the nonparametric two-sample or \eqn{K}-sample Ball Divergence test for
#' equality of multivariate distributions
#' 
#' @aliases sbd.test
#' 
#' @author Xueqin Wang, Jin Zhu
#' 
#' @inheritParams bd.test
#' 
# @param distance.fun A R function for computing the distance. 
# When it is not supplied, distance is computed according to the default 
# \code{dist} in stats package.
#'   
#' 
#' 
#' @return If \code{num.permutations > 0}, \code{sbd.test} returns a \code{htest} class object containing the following components:
#' \item{\code{statistic}}{Ball Divergence statistic.}            
#' \item{\code{p.value}}{the \eqn{p}-value for the test.}
#' \item{\code{replicates}}{permutation replications of the test statistic.}
#' \item{\code{size}}{sample sizes.}
#' \item{\code{complete.info}}{a \code{list} mainly containing two vectors, the first vector is the Ball Divergence statistics 
#' with different aggregation strategy and weight, the second vector is the \eqn{p}-values of tests.}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' If \code{num.permutations = 0}, \code{sbd.test} returns a statistic value.
#' 
#' @rdname sbd.test
#' 
#' @details 
#' \code{sbd.test} is nonparametric test for the two-sample or \eqn{K}-sample problem. 
#' It can detect distribution difference between \eqn{K(K \geq 2)} sample even though sample size are imbalanced.
#' This test can cope well multivariate dataset or complex dataset. 
#' 
#' If only \code{x} is given, the statistic is 
#' computed from the original pooled samples, stacked in 
#' matrix where each row is a multivariate observation, or from the distance matrix 
#' when \code{distance = TRUE}. The first \code{sizes[1]} rows of \code{x} are the first sample, the next 
#' \code{sizes[2]} rows of \code{x} are the second sample, etc.
#' If \code{x} is a \code{list}, its elements are taken as the samples to be compared, 
#' and hence, this \code{list} must contain at least two numeric data vectors, matrices or data.frames.
#' 
#' \code{sbd.test} utilizes the Ball Divergence statistics (see \code{\link{bd}}) to measure dispersion and 
#' derives a \eqn{p}-value via replicating the random permutation \code{num.permutations} times. 
#' The function simply returns the test statistic 
#' when \code{num.permutations = 0}. 
#' 
#' The time complexity of \code{sbd.test} is around \eqn{O(R \times n^2)},
#' where \eqn{R} = \code{num.permutations} and \eqn{n} is sample size.
#' 
#' @note Actually, \code{sbd.test} simultaneously computing \code{"sum"}, \code{"summax"}, and \code{"max"} Ball Divergence statistics 
#' when \eqn{K \geq 3}.
#' Users can get other Ball Divergence statistics and their corresponding \eqn{p}-values 
#' in the \code{complete.info} element of output. We give a quick example below to illustrate. 
#' 
#' @seealso
#' \code{\link{bd}}
#' 
#' @references Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. Ball Divergence: Nonparametric two sample test. Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. https://projecteuclid.org/euclid.aos/1525313077
#' @references Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2021). Ball: An R Package for Detecting Distribution Difference and Association in Metric Spaces, Journal of Statistical Software, Vol.97(6), doi: 10.18637/jss.v097.i06.
#' 
#' @export
#' @examples
#' ################# Quick Start #################
#' library(Ball)
#' num1 <- 10
#' num2 <- 15
#' set.seed(1)
#' x <- rnorm(num1)
#' y <- rnorm(num2, mean = 1)
#' sdist <- sub.dist(x, y, x.sub = 1:num1, y.sub = 1:num2)
#' ## The result of subsampling programming is reasonable because it matches to the result of bd.test
#' sbd.test(sdist)
#' bd.test(dist(c(x, y)), size = c(num1, num2), distance = TRUE)
#' 
#' ################# Imbalance example #################
#' set.seed(1)
#' x <- rnorm(100)
#' y <- rnorm(120, mean = 1)
#' sdist <- sub.dist(x, y, x.sub = 1:30, y.sub = 1:30)
#' sbd.test(sdist, num.permutations = 0)
#' 
#' ################# Large scale example #################
#' num <- 2000
#' set.seed(1)
#' x <- rnorm(num)
#' y <- rnorm(num, mean = 1)
#' sdist <- sub.dist(x, y)
#' system.time(sbd_res <- sbd.test(sdist))
#' system.time(bd_res <- bd.test(dist(c(x, y)), size = c(num, num), distance = TRUE, num.threads = 1))
#' sbd_res
#' bd_res
sbd.test <- function(x, ...) UseMethod("sbd.test")

euclidean_dist <- function(x, y) {
  if (class(x)[1] == "matrix") {
    n1 <- dim(x)[1]
    n2 <- dim(y)[1]
    dist_mat <- matrix(0, nrow = n1, ncol = n2)
    for (i in 1:n1) {
      for (j in 1:n2) {
        dist_mat[i, j] <- norm(x[i, , drop = FALSE] - y[j, , drop = FALSE], 
                               type = "F")
      }
    }    
    # dist_mat <- sqrt(dist_mat)
  } else {
    n1 <- length(x)[1]
    n2 <- length(y)[1]
    dist_mat <- matrix(0, nrow = n1, ncol = n2)
    for (i in 1:n1) {
      for (j in 1:n2) {
        dist_mat[i, j] <- abs(x[i] - y[j])
      }
    }  
  }
  
  dist_mat
}

nlogn <- function(n) {
  n / log(n)
}

#' Title
#'
#' @param x 
#' @param y 
#' @param dist.fun 
#' @param x.sub 
#' @param y.sub 
#' @param sub.scale 
#'
#' @return
#' @export
#'
#' @examples
#' x <- rnorm(100)
#' y <- rnorm(100)
#' sub.dist(x, y)
sub.dist <- function(x, y, dist.fun = NULL, x.sub = NULL, y.sub = NULL, 
                     sub.scale = c("log2", "log", "n/log", "sqrt"), 
                     seed = 1) {
  set.seed(seed)
  x <- as.matrix(x)
  y <- as.matrix(y)
  
  ## TODO: checking
  x_num <- dim(x)[1]
  y_num <- dim(y)[1]
  if (is.null(x.sub) || is.null(y.sub)) {
    sub_scale <- match.arg(sub.scale)
    sub_fun <- switch (sub_scale,
      "log2" = log2,
      "log" = log,
      "nlogn" = nlogn, 
      "sqrt" = sqrt
    )
    
    x_sub_num <- round(sub_fun(x_num))
    x_sub <- sample(1:x_num, size = x_sub_num, replace = FALSE)
    x_sub <- sort(x_sub)
    
    y_sub_num <- round(sub_fun(y_num))
    y_sub <- sample(1:y_num, size = y_sub_num, replace = FALSE)
    y_sub <- sort(y_sub)
  } else {
    x_sub <- x.sub
    y_sub <- y.sub
    x_sub_num <- length(x_sub)
    y_sub_num <- length(y_sub)
  }
  
  if (is.null(dist.fun)) {
    xx <- euclidean_dist(x[x_sub, , drop = FALSE], x)
    xy <- euclidean_dist(x[x_sub, , drop = FALSE], y)
    yx <- euclidean_dist(y[y_sub, , drop = FALSE], x)
    yy <- euclidean_dist(y[y_sub, , drop = FALSE], y)
    dx <- rbind(cbind(xx, xy), cbind(yx, yy))
  } else {
    ## TODO:
  }
  
  sub_dist_info <- list("dist" = dx, 
                        "x.num" = x_num, "y.num" = y_num, 
                        "x.sub" = x_sub, "y.sub" = y_sub,
                        "x.sub.num" = x_sub_num, "y.sub.num" = y_sub_num
                        )
  class(sub_dist_info) <- "sub.dist"
  
  sub_dist_info
}

permute_sub_dist <- function(object) {
  x_num <- object[["x.num"]]
  y_num <- object[["y.num"]]
  num <- x_num + object[["y.num"]]
  permuted_index <- sample(1:num, size = num, replace = FALSE)
  
  tmp <- rep(-1, num)
  sub_index <- c(object[["x.sub"]], x_num + object[["y.sub"]])
  tmp[sub_index] <- sub_index
  tmp <- tmp[permuted_index]
  tmp <- tmp[tmp != -1]
  permute_center_index <- rank(tmp)
  
  object_new <- object
  permuted_dist <- object_new[["dist"]]
  permuted_dist <- permuted_dist[permute_center_index, permuted_index]
  object_new[["dist"]] <- permuted_dist
  
  object_new
}

#' @rdname sbd.test
#' @export
#' @method sbd.test default
sbd.test.default <- function(object, num.permutations = 99, 
                             seed = 1, num.threads = 0, 
                             kbd.type = c("sum", "maxsum", "max"), 
                             weight = c("constant", "variance"), ...) {
  weight <- match.arg(weight)
  kbd.type <- match.arg(kbd.type)
  # if (length(data_name) > 1) {
  #   data_name <- ""
  # }
  
  examine_R_arguments(num.permutations)
  
  ## examine num.thread arguments:
  examine_threads_arguments(num.threads)
  
  ## compute sbd
  x_num <- object[["x.num"]]
  x_sub_num <- object[["x.sub.num"]]
  y_num <- object[["y.num"]]
  y_sub_num <- object[["y.sub.num"]]
  
  dx <- object[["dist"]]
  # xy <- as.double(as.vector(t(dx)))
  # sbd <- as.double(numeric(1))
  # sbd_value <- .C("sbd_C", sbd, xy, 
  #                 as.integer(x_sub_num), as.integer(x_num), 
  #                 as.integer(y_sub_num), as.integer(y_num))[[1]]  
  sbd_value <- sbd_cpp(dx, 
                       as.integer(x_sub_num), as.integer(x_num), 
                       as.integer(y_sub_num), as.integer(y_num))
  # print("sbd:", sbd_value)
  
  ## main:
  if(num.permutations == 0) {
    return(sbd_value);
  } else {
    # permutation method:
    ## examine seed arguments:
    # set.seed(examine_seed_arguments(seed))
    ## hypothesis test:
    permuted_sbd_value <- numeric(num.permutations)
    for (r in 1:num.permutations) {
      set.seed(r)
      object_new <- permute_sub_dist(object)
      dx <- object_new[["dist"]]
      # xy <- as.double(as.vector(t(dx)))
      # sbd <- as.double(numeric(1))
      # permuted_sbd_value[r] <- .C("sbd_C", sbd, xy, 
      #                             as.integer(x_sub_num), as.integer(x_num), 
      #                             as.integer(y_sub_num), as.integer(y_num))[[1]]  
      permuted_sbd_value[r] <- sbd_cpp(dx, 
                                       as.integer(x_sub_num), as.integer(x_num), 
                                       as.integer(y_sub_num), as.integer(y_num))
      # cat(r, ": ", permuted_sbd_value[r], "\n")
    }
    set.seed(NULL)
  }
  # output information:
  # if (result[["info"]][["K"]] == 2) {
  #   stat <- result[["statistic"]][1]
  #   pvalue <- result[["p.value"]][1]
  # } else if (kbd.type == "sum") {
  #   stat <- result[["statistic"]][1]
  #   pvalue <- result[["p.value"]][1]
  #   stat_message <- "sum"
  # } else if (kbd.type == "max") {
  #   stat <- result[["statistic"]][2]
  #   pvalue <- result[["p.value"]][2]
  #   stat_message <- "max"
  # } else {
  #   stat <- result[["statistic"]][3]
  #   pvalue <- result[["p.value"]][3]
  #   stat_message <- "maxsum"
  # }
  # data_name <- paste(data_name, sprintf("\nnumber of observations = %s,", result[["info"]][["N"]]))
  # data_name <- paste(data_name, "group sizes:", paste0(result[["info"]][["size"]], collapse = " "))
  # if (method == "limit") {
  #   null_method <- "Limit Distribution"
  #   data_name <- paste0(data_name, "\nreplicates = ", 0)
  # } else {
  #   null_method <- "Permutation"
  #   data_name <- paste0(data_name, "\nreplicates = ", num.permutations)
  # }
  # data_name <- paste0(data_name, ", weight: ", weight)
  # if (result[["info"]][["K"]] == 3) {
  #   data_name <- paste0(data_name, ", kbd.type: ", stat_message)
  # }
  # data_name <- paste0(data_name, ", Weighted Ball Divergence = ", result[["info"]][["weight"]])
  alternative_message <- "distributions of samples are distinct"
  
  # return:
  # e <- list(
  #   statistic = stat,
  #   p.value = pvalue,
  #   replicates = num.permutations,
  #   size = result[["info"]][["size"]],
  #   complete.info = result[["info"]],
  #   alternative = alternative_message,
  #   method = sprintf("%s-sample Ball Divergence Test (%s)", result[["info"]][["K"]], null_method),
  #   data.name = data_name
  # )
  e <- list(
    statistic = sbd_value,
    permuted_statistic = permuted_sbd_value, 
    p.value = mean(c(permuted_sbd_value >= sbd_value, TRUE)),
    replicates = num.permutations, 
    alternative = alternative_message
  )
  class(e) <- "htest"
  return(e)
}


#' @rdname sbd.test
#'
#' @param formula a formula of the form \code{response ~ group} where \code{response} gives the data values and \code{group} a vector or factor of the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see \code{model.frame}) containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @method sbd.test formula
#' @export
#' @examples
#' 
#' ################  Formula interface  ################
#' ## Two-sample test
#' sbd.test(extra ~ group, data = sleep)
#' ## K-sample test
#' sbd.test(Sepal.Width ~ Species, data = iris)
#' sbd.test(Sepal.Width ~ Species, data = iris, kbd.type = "max")
sbd.test.formula <- function(formula, data, subset, na.action, ...) {
  if(missing(formula)
     || (length(formula) != 3L)
     || (length(attr(terms(formula[-2L]), "term.labels")) != 1L))
    stop("'formula' missing or incorrect")
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  ## need stats:: for non-standard evaluation
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, parent.frame())
  DNAME <- paste(names(mf), collapse = " by ")
  names(mf) <- NULL
  response <- attr(attr(mf, "terms"), "response")
  g <- factor(mf[[-response]])
  if(nlevels(g) < 2L)
    stop("grouping factor must contain at least two levels")
  DATA <- list()
  DATA[["x"]] <- split(mf[[response]], g)
  y <- do.call("sbd.test", c(DATA, list(...)))
  remind_info <- strsplit(y$data.name, split = "number of observations")[[1]][2]
  DNAME <- paste0(DNAME, "\nnumber of observations")
  y$data.name <- paste0(DNAME, remind_info)
  y
}

