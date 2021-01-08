#' @title Ball Divergence based Equality of Distributions Test
#' 
#' @description Performs the nonparametric two-sample or \eqn{K}-sample Ball Divergence test for
#' equality of multivariate distributions
#' 
#' @aliases bd.test
#' 
#' @author Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang, Jin Zhu
#' 
#' @param x a numeric vector, matrix, data.frame, or a list containing at least two numeric vectors, matrices, or data.frames.
#' @param y a numeric vector, matrix, data.frame.
#' @param num.permutations the number of permutation replications. When \code{num.permutations = 0}, the function just returns
#' the Ball Divergence statistic. Default: \code{num.permutations = 99}.
#' @param distance if \code{distance = TRUE}, the elements of \code{x} will be considered as a distance matrix. Default: \code{distance = FALSE}.
#' @param size a vector recording sample size of each group.
#' @param seed the random seed. Default \code{seed = 1}.
#' @param num.threads number of threads. If \code{num.threads = 0}, then all of available cores will be used. Default \code{num.threads = 0}.
#' @param kbd.type a character string specifying the \eqn{K}-sample Ball Divergence test statistic, 
#' must be one of \code{"sum"}, \code{"summax"}, or \code{"max"}. Any unambiguous substring can be given. 
#' Default \code{kbd.type = "sum"}.
#' @param method if \code{method = "permutation"}, a permutation procedure is carried out to compute the \eqn{p}-value;
#' if \code{ method = "limit"}, an approximate null distribution is used when \code{weight = "constant"}.
#' Any unambiguous substring can be given. Default \code{method = "permutation"}.
#' @param weight a character string specifying the weight form of Ball Divergence statistic.
#' It must be one of \code{"constant"} or \code{"variance"}. 
#' Any unambiguous substring can be given. Default: \code{weight = "constant"}.
#' @param ... further arguments to be passed to or from methods.
#' 
## @param weight not available now
#' 
#' @return If \code{num.permutations > 0}, \code{bd.test} returns a \code{htest} class object containing the following components:
#' \item{\code{statistic}}{Ball Divergence statistic.}            
#' \item{\code{p.value}}{the \eqn{p}-value for the test.}
#' \item{\code{replicates}}{permutation replications of the test statistic.}
#' \item{\code{size}}{sample sizes.}
#' \item{\code{complete.info}}{a \code{list} mainly containing two vectors, the first vector is the Ball Divergence statistics 
#' with different aggregation strategy and weight, the second vector is the \eqn{p}-values of tests.}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' If \code{num.permutations = 0}, \code{bd.test} returns a statistic value.
#' 
#' @rdname bd.test
#' 
#' @details 
#' \code{bd.test} is nonparametric test for the two-sample or \eqn{K}-sample problem. 
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
#' \code{bd.test} utilizes the Ball Divergence statistics (see \code{\link{bd}}) to measure dispersion and 
#' derives a \eqn{p}-value via replicating the random permutation \code{num.permutations} times. 
#' The function simply returns the test statistic 
#' when \code{num.permutations = 0}. 
#' 
#' The time complexity of \code{bd.test} is around \eqn{O(R \times n^2)},
#' where \eqn{R} = \code{num.permutations} and \eqn{n} is sample size.
#' 
#' @note Actually, \code{bd.test} simultaneously computing \code{"sum"}, \code{"summax"}, and \code{"max"} Ball Divergence statistics 
#' when \eqn{K \geq 3}.
#' Users can get other Ball Divergence statistics and their corresponding \eqn{p}-values 
#' in the \code{complete.info} element of output. We give a quick example below to illustrate. 
#' 
#' @seealso
#' \code{\link{bd}}
#' 
#' @references Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. Ball Divergence: Nonparametric two sample test. Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. https://projecteuclid.org/euclid.aos/1525313077
#' @references Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2018). Ball: An R package for detecting distribution difference and association in metric spaces. arXiv preprint arXiv:1811.03750. http://arxiv.org/abs/1811.03750
#' 
#' @export
#' @examples
#' ################# Quick Start #################
#' set.seed(1)
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 1)
#' # plot(density(x))
#' # lines(density(y), col = "red")
#' bd.test(x = x, y = y)
#' 
#' ################# Quick Start #################
#' x <- matrix(rnorm(100), nrow = 50, ncol = 2)
#' y <- matrix(rnorm(100, mean = 3), nrow = 50, ncol = 2)
#' # Hypothesis test with Standard Ball Divergence:
#' bd.test(x = x, y = y)
#'
#' ################# Simlated Non-Hilbert data #################
#' data("bdvmf")
#' \dontrun{
#' library(scatterplot3d)
#' scatterplot3d(bdvmf[["x"]], color = bdvmf[["group"]], 
#'               xlab = "X1", ylab = "X2", zlab = "X3")
#' }
#' # calculate geodesic distance between sample:
#' Dmat <- nhdist(bdvmf[["x"]], method = "geodesic")
#' # hypothesis test with BD :
#' bd.test(x = Dmat, size = c(150, 150), num.permutations = 99, distance = TRUE)
#'
#' ################# Non-Hilbert Real Data #################
#' # load data:
#' data("macaques")
#' # number of femala and male Macaca fascicularis:
#' table(macaques[["group"]])
#' # calculate Riemannian shape distance matrix:
#' Dmat <- nhdist(macaques[["x"]], method = "riemann")
#' # hypothesis test with BD:
#' bd.test(x = Dmat, num.permutations = 99, size = c(9, 9), distance = TRUE)
#' 
#' ################  K-sample Test  #################
#' n <- 150
#' bd.test(rnorm(n), size = c(40, 50, 60))
#' # alternative input method:
#' x <- lapply(c(40, 50, 60), rnorm)
#' res <- bd.test(x)
#' res
#' ## get all Ball Divergence statistics:
#' res[["complete.info"]][["statistic"]]
#' ## get all test result:
#' res[["complete.info"]][["p.value"]]
#' 
#' ################  Testing via approximate limit distribution  #################
#' \dontrun{
#' set.seed(1)
#' n <- 1000
#' x <- rnorm(n)
#' y <- rnorm(n)
#' res <- bd.test(x, y, method = "limit")
#' bd.test(x, y)
#' }
bd.test <- function(x, ...) UseMethod("bd.test")


#' @rdname bd.test
#' @export
#' @method bd.test default
bd.test.default <- function(x, y = NULL, num.permutations = 99, 
                            method = c("permutation", "limit"), distance = FALSE,
                            size = NULL, seed = 1, num.threads = 0, 
                            kbd.type = c("sum", "maxsum", "max"), 
                            weight = c("constant", "variance"), ...) {
  weight <- match.arg(weight)
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  kbd.type <- match.arg(kbd.type)
  method <- match.arg(method)
  if (length(data_name) > 1) {
    data_name <- ""
  }
  if(is.null(x) || is.null(y)) {
    
    if(is.null(x) & is.null(y)) {
      stop("x and y are all null!")
    }
    
    # modify input information:
    data_name <- gsub(x = data_name, pattern = " and NULL", replacement = "")
    
    if (class(x)[1] == "dist") {
      distance <- TRUE
    }
    if(distance) {
      examine_size_arguments(size)
      if (length(size) >= 2) {
        if (class(x)[1] == "dist") {
          if (attr(x, "Size") != sum(size)) { stop("size arguments is error!") }
          xy <- as.vector(x)
        } else {
          if (nrow(x) != sum(size)) { stop("size arguments is error!") }
          xy <- x[lower.tri(x)]
        }
      } else if (length(size) < 2) {
        stop("size arguments is error!")
      }
    } else {
      if (is.list(x)) {
        x <- lapply(x, as.matrix)
        if (length(unique(sapply(x, ncol))) != 1) {
          stop("data with different dimension!")
        }
        size <- sapply(x, nrow)
        x <- do.call("rbind", x)
        p <- ncol(x)
      } else if (is.vector(x)) {
        p <- 1
      } else {
        p <- ncol(x)
        if (p == 1) {
          x <- as.vector(x)
        }
      }
      if (p > 1) {
        xy <- dist(x)
        if (attr(xy, "Size") != sum(size)) { stop("size arguments is error!") }
        xy <- as.vector(xy)
        distance <- TRUE
      } else {
        xy <- x
        distance <- FALSE
      }
    }
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    p <- examine_dimension(x, y)
    # 
    if(p > 1 || method == "limit") {
      xy <- get_vectorized_distance_matrix(x, y)
      distance <- TRUE
      size <- c(xy[[2]], xy[[3]])
      xy <- xy[[1]]
    } else {
      xy <- c(x, y)
      distance <- FALSE
      size <- c(dim(x)[1], dim(y)[1])
    }
  }
  
  ## memory protect step:
  # memoryAvailable(n = sum(size), funs = 'BD.test')
  
  ## examine num.permutations arguments:
  if(method == "limit") {
    num.permutations <- 0
  } else {
    examine_R_arguments(num.permutations)
  }
  
  ## examine num.thread arguments:
  examine_threads_arguments(num.threads)
  
  ## main:
  if(num.permutations == 0) {
    result <- bd_test_wrap_c(xy, size, num.permutations = 0, weight, distance, num.threads)
    # approximately method:
    if(method == "limit") {
      if(result[["info"]][["K"]] == 2) {
        eigenvalue <- bd_limit_wrap_c(xy, size, distance, num.threads)
        result[["p.value"]] <- 1 - hbe(eigenvalue, prod(size) * result[["statistic"]] / sum(size))
      } else {
        return(result[["statistic"]])
      }
    } 
    # return statistic when num.permutations = 0:
    else {
      if (result[["info"]][["K"]] == 2) {
        return_stat <- result[["statistic"]][1] 
      } else {
        if (kbd.type == "sum") {
          return_stat <- result[["statistic"]][1]
        } else if (kbd.type == "max") {
          return_stat <- result[["statistic"]][2] 
        } else if (kbd.type == "maxsum") {
          return_stat <- result[["statistic"]][3] 
        }
      }
      return(return_stat)
    }
  } else {
    # permutation method:
    ## examine seed arguments:
    set.seed(examine_seed_arguments(seed))
    ## hypothesis test:
    result <- bd_test_wrap_c(xy, size, num.permutations, 
                             weight, distance, num.threads)
    # pvalue <- calculatePvalue(result[["statistic"]], result[["permuted_stat"]])
    set.seed(NULL)
  }
  # output information:
  if (result[["info"]][["K"]] == 2) {
    stat <- result[["statistic"]][1]
    pvalue <- result[["p.value"]][1]
  } else if (kbd.type == "sum") {
    stat <- result[["statistic"]][1]
    pvalue <- result[["p.value"]][1]
    stat_message <- "sum"
  } else if (kbd.type == "max") {
    stat <- result[["statistic"]][2]
    pvalue <- result[["p.value"]][2]
    stat_message <- "max"
  } else {
    stat <- result[["statistic"]][3]
    pvalue <- result[["p.value"]][3]
    stat_message <- "maxsum"
  }
  data_name <- paste(data_name, sprintf("\nnumber of observations = %s,", result[["info"]][["N"]]))
  data_name <- paste(data_name, "group sizes:", paste0(result[["info"]][["size"]], collapse = " "))
  if (method == "limit") {
    null_method <- "Limit Distribution"
    data_name <- paste0(data_name, "\nreplicates = ", 0)
  } else {
    null_method <- "Permutation"
    data_name <- paste0(data_name, "\nreplicates = ", num.permutations)
  }
  data_name <- paste0(data_name, ", weight: ", weight)
  if (result[["info"]][["K"]] == 3) {
    data_name <- paste0(data_name, ", kbd.type: ", stat_message)
  }
  # data_name <- paste0(data_name, ", Weighted Ball Divergence = ", result[["info"]][["weight"]])
  alternative_message <- "distributions of samples are distinct"
  
  # return:
  e <- list(
    statistic = stat,
    p.value = pvalue,
    replicates = num.permutations,
    size = result[["info"]][["size"]],
    complete.info = result[["info"]],
    alternative = alternative_message,
    method = sprintf("%s-sample Ball Divergence Test (%s)", result[["info"]][["K"]], null_method),
    data.name = data_name
  )
  class(e) <- "htest"
  return(e)
}


#' @rdname bd.test
#'
#' @param formula a formula of the form \code{response ~ group} where \code{response} gives the data values and \code{group} a vector or factor of the corresponding groups.
#' @param data an optional matrix or data frame (or similar: see \code{model.frame}) containing the variables in the formula \code{formula}. By default the variables are taken from \code{environment(formula)}.
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @method bd.test formula
#' @export
#' @examples
#' 
#' ################  Formula interface  ################
#' ## Two-sample test
#' bd.test(extra ~ group, data = sleep)
#' ## K-sample test
#' bd.test(Sepal.Width ~ Species, data = iris)
#' bd.test(Sepal.Width ~ Species, data = iris, kbd.type = "max")
bd.test.formula <- function(formula, data, subset, na.action, ...) {
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
  y <- do.call("bd.test", c(DATA, list(...)))
  remind_info <- strsplit(y$data.name, split = "number of observations")[[1]][2]
  DNAME <- paste0(DNAME, "\nnumber of observations")
  y$data.name <- paste0(DNAME, remind_info)
  y
}


#' @title Ball Divergence statistic
#' @description Compute Ball Divergence statistic, which is a generic dispersion measure in Banach spaces.
#' @author Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang
#' @inheritParams bd.test
#' @rdname bd
#' @return 
#' \item{\code{bd }}{ Ball Divergence statistic}
#' 
#' @details 
#' Given the samples not containing missing values, \code{bd} returns Ball Divergence statistics.
#' If we set \code{distance = TRUE}, arguments \code{x}, \code{y} can be a \code{dist} object or a
#' symmetric numeric matrix recording distance between samples; 
#' otherwise, these arguments are treated as data.
#' 
#' Ball divergence statistic measure the distribution difference of two datasets in Banach spaces. 
#' The Ball divergence statistic is proven to be zero if and only if two datasets are identical.
#' 
#' The definition of the Ball Divergence statistics is as follows.
#' Given two independent samples \eqn{ \{x_{1}, \ldots, x_{n}\} } with the associated probability measure \eqn{\mu} and 
#' \eqn{ \{y_{1}, \ldots, y_{m}\} } with \eqn{\nu}, where the observations in each sample are \emph{i.i.d}.
#' Let \eqn{\delta(x,y,z)=I(z\in \bar{B}(x, \rho(x,y)))}, 
#' where \eqn{\delta(x,y,z)} indicates whether \eqn{z} is located in the closed ball \eqn{\bar{B}(x, \rho(x,y))} 
#' with center \eqn{x} and radius \eqn{\rho(x, y)}. 
#' We denote:
#' \deqn{
#' A_{ij}^{X}=\frac{1}{n}\sum_{u=1}^{n}{\delta(X_i,X_j,X_u)}, \quad A_{ij}^{Y}=\frac{1}{m}\sum_{v=1}^{m}{\delta(X_i,X_j,Y_v)},
#' }
#' \deqn{
#' C_{kl}^{X}=\frac{1}{n}\sum_{u=1}^{n}{\delta(Y_k,Y_l,X_u)}, \quad C_{kl}^{Y}=\frac{1}{m}\sum_{v=1}^{m}{\delta(Y_k,Y_l,Y_v)}.
#' }
#' \eqn{A_{ij}^X} represents the proportion of samples \eqn{ \{x_{1}, \ldots, x_{n}\} } located in the 
#' ball \eqn{\bar{B}(X_i,\rho(X_i,X_j))} and \eqn{A_{ij}^Y} represents the proportion of samples \eqn{ \{y_{1}, \ldots, y_{m}\} } 
#' located in the ball \eqn{\bar{B}(X_i,\rho(X_i,X_j))}. 
#' Meanwhile, \eqn{C_{kl}^X} and \eqn{C_{kl}^Y} represent the corresponding proportions located in the ball \eqn{\bar{B}(Y_k,\rho(Y_k,Y_l))}.
#' The Ball Divergence statistic is defined as:
#' \deqn{D_{n,m}=A_{n,m}+C_{n,m}}
#' 
#' Ball Divergence can be generalized to the \emph{K}-sample test problem. Suppose we 
#' have \eqn{K} group samples, each group include \eqn{n_{k}} samples. 
#' The definition of \eqn{K}-sample Ball Divergence statistic could be 
#' to directly sum up the two-sample Ball Divergence statistics of all sample pairs (\code{kbd.type = "sum"})
#' \deqn{\sum_{1 \leq k < l \leq K}{D_{n_{k},n_{l}}},}
#' or to find one sample with the largest difference to the others (\code{kbd.type = "maxsum"})
#' \deqn{\max_{t}{\sum_{s=1, s \neq t}^{K}{D_{n_{s}, n_{t}}},}}
#' to aggregate the \eqn{K-1} most significant different two-sample Ball Divergence statistics (\code{kbd.type = "max"})
#' \deqn{\sum_{k=1}^{K-1}{D_{(k)}},}
#' where \eqn{D_{(1)}, \ldots, D_{(K-1)}} are the largest \eqn{K-1} two-sample Ball Divergence statistics among 
#' \eqn{\{D_{n_s, n_t}| 1 \leq s < t \leq K\}}. When \eqn{K=2},
#' the three types of Ball Divergence statistics degenerate into two-sample Ball Divergence statistic.
#' 
#' See \code{\link{bd.test}} for a test of distribution equality based on the Ball Divergence.
#' 
#' @seealso
#' \code{\link{bd.test}}
#' @export
#' 
#' @references Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. Ball Divergence: Nonparametric two sample test. Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. https://projecteuclid.org/euclid.aos/1525313077
#' 
#' @examples
#' ############# Ball Divergence #############
#' x <- rnorm(50)
#' y <- rnorm(50)
#' bd(x, y)
bd <- function(x, y = NULL, distance = FALSE, size = NULL, num.threads = 1, kbd.type = c("sum", "maxsum", "max")) {
  res <- bd.test(x = x, y = y, distance = distance, size = size, num.permutations = 0, kbd.type = kbd.type)
  res
}
