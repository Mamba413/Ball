#' @title Ball Covariance Test
#' @author Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu, Jin Zhu
#' @description Ball Covariance test of independence. 
#' Ball covariance are generic dependence measures in Banach spaces.
#' 
#' @inheritParams bd.test
#' @param x a numeric vector, matrix, data.frame, or a list containing at least two numeric vectors, matrices, or data.frames.
#' @param y a numeric vector, matrix, or data.frame.
#' @param num.permutations the number of permutation replications. When \code{num.permutations = 0}, the function just returns
#' the Ball Covariance statistic. Default: \code{num.permutations = 99}.
#' @param distance if \code{distance = TRUE}, the elements of \code{x} and \code{y} are considered as distance matrices.
#' @param weight a logical or character string used to choose the weight form of Ball Covariance statistic.. 
#' If input is a character string, it must be one of \code{"constant"}, \code{"probability"}, or \code{"chisquare"}. 
#' Any unambiguous substring can be given. 
#' If input is a logical value, it is equivalent to \code{weight = "probability"} if \code{weight = TRUE} while 
#' equivalent to \code{weight = "constant"} if \code{weight = FALSE}.
#' Default: \code{weight = FALSE}.
#' @param method if \code{method = "permutation"}, a permutation procedure will be carried out;
#' if \code{method = "limit"}, the p-values based on approximate Ball Covariance distribution are given
#' (Only available for Independence test and \code{weight = "constant"}).
#' 
#' @return If \code{num.permutations > 0}, \code{bcov.test} returns a \code{htest} class object containing the following components:
#' \item{\code{statistic}}{Ball Covariance statistic.}            
#' \item{\code{p.value}}{the p-value for the test.}  
#' \item{\code{replicates}}{permutation replications of the test statistic.}
#' \item{\code{size}}{sample size.} 
#' \item{\code{complete.info}}{a \code{list} mainly containing two vectors, the first vector is the Ball Covariance statistics 
#' with different weights, the second is the \eqn{p}-values of weighted Ball Covariance tests.}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' If \code{num.permutations = 0}, \code{bcov.test} returns a statistic value.
#' 
#' @details 
#' \code{bcov.test} are non-parametric tests of independence in Banach spaces. 
#' It can detect the dependence between two random objects (variables) and 
#' the mutual dependence among at least three random objects (variables).
#' 
#' If two samples are pass to arguments \code{x} and \code{y}, the sample sizes (i.e. number of rows or length of the vector) 
#' of the two variables must agree. If a \code{\link{list}} object is passed to \code{x}, this \code{list} must contain at least
#' two numeric vectors, matrices, or data.frames, and each element of this \code{list} 
#' must with the same sample size. Moreover, data pass to \code{x} or \code{y} 
#' must not contain missing or infinite values. 
#' If \code{distance = TRUE}, \code{x} is considered as a distance matrix or a list containing distance matrices, 
#' and \code{y} is considered as a distance matrix; otherwise, these arguments are treated as data.
#' 
#' \code{bcov.test} utilizes the Ball Covariance statistics (see \code{\link{bcov}}) to measure dependence and 
#' derives a \eqn{p}-value via replicating the random permutation \code{num.permutations} times.
#' 
#' See Pan et al 2018 for theoretical properties of the test, including statistical consistency.
#' 
#' @note Actually, \code{bcov.test} simultaneously computing Ball Covariance statistics with 
#' \code{"constant"}, \code{"probability"}, and \code{"chisquare"} weights.
#' Users can get other Ball Covariance statistics with different weight and their corresponding \eqn{p}-values 
#' in the \code{complete.info} element of output. We give a quick example below to illustrate. 
#' 
#' @references Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019) Ball Covariance: A Generic Measure of Dependence in Banach Space, Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1543600
#' @references Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2018). Ball: An R package for detecting distribution difference and association in metric spaces. arXiv preprint arXiv:1811.03750. http://arxiv.org/abs/1811.03750
#' 
#' @rdname bcov.test
#' 
#' @useDynLib Ball, .registration = TRUE
#' @export
#' @seealso \code{\link{bcov}}, \code{\link{bcor}}
#' @examples
#' set.seed(1)
#' 
#' ################# Quick Start #################
#' error <- runif(50, min = -0.3, max = 0.3)
#' x <- runif(50, 0, 4*pi)
#' y <- cos(x) + error
#' # plot(x, y)
#' res <- bcov.test(x, y)
#' res
#' ## get all Ball Covariance statistics:
#' res[["complete.info"]][["statistic"]]
#' ## get all test result:
#' res[["complete.info"]][["p.value"]]
#' 
#' ################# Quick Start #################
#' x <- matrix(runif(50 * 2, -pi, pi), nrow = 50, ncol = 2)
#' error <- runif(50, min = -0.1, max = 0.1)
#' y <- sin(x[,1] + x[,2]) + error
#' bcov.test(x = x, y = y, weight = "prob")
#' 
#' ################# Ball Covariance Test for Non-Hilbert Data #################
#' # load data:
#' data("ArcticLake")
#' # Distance matrix between y:
#' Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
#' # Distance matrix between x:
#' Dx <- dist(ArcticLake[["depth"]])
#' # hypothesis test with BCov:
#' bcov.test(x = Dx, y = Dy, distance = TRUE)
#' 
#' ################  Weighted Ball Covariance Test  #################
#' data("ArcticLake")
#' Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
#' Dx <- dist(ArcticLake[["depth"]])
#' # hypothesis test with weighted BCov:
#' bcov.test(x = Dx, y = Dy, distance = TRUE, weight = "prob")
#' 
#' ################# Mutual Independence Test #################
#' x <- rnorm(50)
#' y <- (x > 0) * x + rnorm(50)
#' z <- (x <= 0) * x + rnorm(50)
#' data_list <- list(x, y, z)
#' bcov.test(data_list)
#' data_list <- lapply(data_list, function(x) {
#'   as.matrix(dist(x))
#' })
#' bcov.test(data_list, distance = TRUE)
#' bcov.test(data_list, distance = FALSE, weight = "chi")
#' 
#' ################# Mutual Independence Test for Meteorology data #################
#' data("meteorology")
#' bcov.test(meteorology)
#' 
#' ################  Testing via approximate limit distribution  #################
#' set.seed(1)
#' n <- 200
#' p <- 2
#' x <- matrix(rnorm(n * p), nrow = n)
#' y <- matrix(rnorm(n * p), nrow = n)
#' bcov.test(x, y, method = "limit")
#' bcov.test(x, y)
#' 
bcov.test <- function(x, ...) UseMethod("bcov.test")


#' @rdname bcov.test
#' @export
#' @method bcov.test default
bcov.test.default <- function(x, y = NULL, num.permutations = 99, 
                              method = c("permutation", "limit"), 
                              distance = FALSE, weight = FALSE, 
                              seed = 1, num.threads = 0, ...)
{
  method <- match.arg(method)
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  if (length(data_name) > 1) {
    data_name <- ""
  }
  type <- "bcov"
  # modify input information:
  if(class(x) == "list") {
    data_name <- gsub(x = data_name, pattern = " and NULL", replacement = "")
  }
  weight <- examine_weight_arguments(weight)
  result <- bcov_test_internal_wrap(x = x, y = y, num.permutations = num.permutations, 
                                    distance = distance, weight = weight, 
                                    seed = seed, method = method, type = type, 
                                    num.threads = num.threads)
  # return result of hypothesis test:
  if(num.permutations == 0) {
    return(result)
  } else {
    if (weight == "constant") {
      stat <- result[["statistic"]][1]
      pvalue <- result[["p.value"]][1]
      weight_name <- "constant"
    } else if (weight == "probability") {
      stat <- result[["statistic"]][2]
      pvalue <- result[["p.value"]][2]
      weight_name <- "probability"
    } else if (weight == "chisquare") {
      stat <- result[["statistic"]][3]
      pvalue <- result[["p.value"]][3]
      weight_name <- "chisquare"
    } 
    
    data_name <- paste0(data_name,"\nnumber of observations = ", result[["info"]][["N"]])
    data_name <- paste0(data_name, "\nreplicates = ", num.permutations, 
                        ", weight: ", weight_name)
    test_method <- "Ball Covariance test of %sindependence"
    test_type <- ifelse(class(x) == "list" && length(x) > 2, "mutual ", "")
    test_method <- sprintf(test_method, test_type)
    # if(type == "bcor") {
    #   test_method <- gsub(pattern = "Covariance", replacement = "Correlation", x = test_method)
    #   data_name <- gsub(pattern = "Covariance", replacement = "Correlation", x = data_name)
    # }
    alternative_message <- "random variables are dependent"
    e <- list(
      statistic = stat,
      p.value = pvalue,
      replicates = num.permutations,
      size = result[["info"]][["N"]],
      complete.info = result,
      alternative = alternative_message,
      method = test_method,
      data.name = data_name
    )
    class(e) <- "htest"
    return(e)
  }
}



#' @rdname bcov.test
#'
#' @param formula a formula of the form \code{~ u + v}, where each of \code{u} and \code{v} are numeric variables giving the data values for one sample. The samples must be of the same length.
#' @param data an optional matrix or data frame (or similar: see \code{model.frame}) containing the variables in the formula formula. By default the variables are taken from environment(formula).
#' @param subset an optional vector specifying a subset of observations to be used.
#' @param na.action a function which indicates what should happen when the data contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param ... further arguments to be passed to or from methods.
#' 
#' @export
#' @method bcov.test formula
#' @importFrom stats model.frame
#'
#' @examples
#' 
#' ################  Formula interface  ################
#' ## independence test:
#' bcov.test(~ CONT + INTG, data = USJudgeRatings)
#' ## independence test with chisquare weight:
#' bcov.test(~ CONT + INTG, data = USJudgeRatings, weight = "chi")
#' ## mutual independence test:
#' bcov.test(~ CONT + INTG + DMNR, data = USJudgeRatings)
bcov.test.formula <- function(formula, data, subset, na.action, ...) {
  if(missing(formula)
     || !inherits(formula, "formula")
     || length(formula) != 2L)
    stop("'formula' missing or invalid")
  m <- match.call(expand.dots = FALSE)
  if(is.matrix(eval(m$data, parent.frame())))
    m$data <- as.data.frame(data)
  m[[1L]] <- quote(stats::model.frame)
  m$... <- NULL
  mf <- eval(m, environment(formula))
  if(length(mf) < 2L)
    stop("invalid formula")
  DNAME <- paste(names(mf), collapse = " and ")
  dat <- list()
  dat[["x"]] <- as.list(mf)
  y <- do.call("bcov.test", c(dat, list(...)))
  remind_info <- strsplit(y$data.name, split = "number of observations")[[1]][2]
  DNAME <- paste0(DNAME, "\nnumber of observations")
  y$data.name <- paste0(DNAME, remind_info)
  y
}


#' Ball covariance test internal function
#' @inheritParams bcov.test
#' @param type 
#'
#' @noRd
bcov_test_internal <- function(x, y, num.permutations = 99, distance = FALSE, weight = FALSE, 
                               seed = 4, method = 'permute', num.threads)
{
  if (class(x) == "dist" && class(y) == "dist") {
    distance <- TRUE
  }
  if (distance) {
    if (class(x) == "dist" || class(y) == "dist") {
      if (class(x) != "dist") {
        x <- as.dist(x)
      }
      if (class(y) != "dist") {
        y <- as.dist(y)
      }
      num <- attr(x, "Size")
      x <- as.vector(x)
      y <- as.vector(y)
    } else {
      num <- nrow(x)
      x <- x[lower.tri(x)]
      y <- y[lower.tri(y)]
    }
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    num <- nrow(x)
    if (ncol(x) != 1 || ncol(y) != 1) {
      x <- as.vector(dist(x))
      y <- as.vector(dist(y))
      distance <- TRUE
    } else {
      x <- as.vector(x)
      y <- as.vector(y)
    }
  }
  examine_x_y_bcov(x, y)
  
  ## memory protect step:
  # memoryAvailable(num, funs = 'BI.test')
  ## examine test type:
  # type <- examine_type_arguments(type)
  ## examine num.permutations arguments:
  if(method == "limit") {
    num.permutations <- 0
  } else {
    examine_R_arguments(num.permutations)
  }
  #
  if(num.permutations == 0) {
    result <- bcov_test_wrap_c(x = x, y = y, n = num, num.permutations = 0, 
                               distance = distance, num.threads = num.threads)
    if(method == "limit") {
      eigenvalue <- bcov_limit_wrap_c(x, y, num, distance, num.threads)
      result[["p.value"]] <- 1 - hbe(eigenvalue, num * result[["statistic"]][1])
      return(result)
    } else {
      if (weight == WEIGHT_TYPE[1]) {
        return(result[[1]][1])
      } else if (weight == WEIGHT_TYPE[2]) {
        return(result[[1]][2])
      } else if (weight == WEIGHT_TYPE[3]) {
        return(result[[1]][3])
      }
    }
  } else {
    set.seed(seed = examine_seed_arguments(seed))
    result <- bcov_test_wrap_c(x = x, y = y, n = num, num.permutations = num.permutations, distance = distance, num.threads = num.threads)
    return(result)
  }
}


#' A internal function for carry out independence test for multiple random variables
#' @inheritParams bcov.test
#' @inherit return
#' @noRd
kbcov_test_internal <- function(x, num.permutations = 99, distance = FALSE, weight = FALSE, 
                                seed = 1, method = 'permute', num.threads)
{
  ############################################################
  #################### R Version (1.1.0) #####################
  ############################################################
  # x <- lapply(x, as.matrix)
  # size_list <- sapply(x, nrow)
  # num <- unique(size_list)
  # if(length(num) > 1) {
    # stop("sample sizes of variables are not match!")
  # }
  # if(distance) {
  #   
  # } else {
  #   x <- lapply(x, dist, diag = TRUE, upper = TRUE)
  #   x <- lapply(x, as.matrix)
  # }
  # var_num <- length(x)
  # # compute statistic:
  # stat_value <- kbcov_stat(x = x, num = num, var_num = var_num, 
  #                          weight = weight, type = type)
  # if(num.permutations == 0) {
  #   names(stat_value)
  #   return(stat_value)
  # } else {
  #   seed <- examine_seed_arguments(seed)
  #   set.seed(seed)
  #   # permutation procedure:
  #   permuted_stat <- matrix(nrow = 3, ncol = num.permutations)
  #   for (r in 1:num.permutations) {
  #     x_copy <- x
  #     for (v in 1:var_num) {
  #       index <- sample(1:num, size = num, replace = FALSE)
  #       x_copy[[v]] <- x[[v]][index, index]
  #     }
  #     permuted_stat[, r] <- kbcov_stat(x = x_copy, num = num, 
  #                                      var_num = var_num, 
  #                                      weight = weight, type = type)
  #   }
  #   permuted_stat <- t(permuted_stat)
  #   # calculate pvalue:
  #   pvalue <- sapply(1:3, function(i) {
  #     calculatePvalue(stat_value[i], permuted_stat[, i])
  #   })
  #   names(pvalue) <- paste0(names(stat_value), ".pvalue")
  #   # pvalue <- calculatePvalue(stat_value, permuted_stat)
  # }
  # return result:
  # list("statistic" = stat_value, 
  #      "p.value" = pvalue, 
  #      "info" = list("N" = num, "num.permutations" = num.permutations))

  
  ############################################################
  #################### C Version (1.2.0) #####################
  ############################################################
  var_num <- length(x)
  if ((!distance) && all(sapply(x, class) == "dist")) {
    distance <- TRUE
  }
  if(distance) {
    if (class(x[[1]]) != "dist") {
      if (nrow(x[[1]]) != ncol(x[[1]])) {
        stop("The elements of input list is not a distance matrix.")
      }
    }
  } else {
    x <- lapply(x, dist)
  }
  distance <- TRUE
  num <- ifelse(class(x[[1]]) == "dist", attr(x[[1]], "Size"), nrow(x[[1]]))
  if (class(x[[1]]) == "dist") {
    x <- lapply(x, as.vector)
  } else {
    x <- lapply(x, function(xx) { 
      xx[lower.tri(xx)] 
    })
  }
  if (length(unique(sapply(x, length))) != 1) {
    stop("sample size of variables are not match!")
  } else {
    old_x <- x
    x <- unlist(x)
  }
  #
  if(method == "limit") {
    num.permutations <- 0
  } else {
    examine_R_arguments(num.permutations)
  }
  #
  if(num.permutations == 0) {
    result <- kbcov_test_wrap_c(x = x, K = var_num, n = num, num.permutations = 0, 
                                distance = distance, num.threads = num.threads)
    if(method == "limit") {
      eigenvalue <- bcov_limit_wrap_c(old_x, NULL, num, distance, num.threads)
      result[["p.value"]] <- 1 - hbe(eigenvalue, num * result[["statistic"]])
      return(result)
    } else {
      if (weight == WEIGHT_TYPE[1]) {
        return(result[[1]][1])
      } else if (weight == WEIGHT_TYPE[2]) {
        return(result[[1]][2])
      } else {
        return(result[[1]][3])
      }
    }
  } else {
    set.seed(seed = examine_seed_arguments(seed))
    result <- kbcov_test_wrap_c(x = x, K = var_num, n = num, num.permutations = num.permutations, 
                                distance = distance, num.threads = num.threads)
    return(result)
  }
  
}


#' compute extension of BCov for independence of multiple random variables
#'
#' @param x list containing distance matrix, each element is distance matrix
#' @param num sample size
#' @param var_num random variables number
#' @param weight whether used weight
#' @param type Ball Correlation or Ball Covariance. now, only Ball Covariance is considered.
#'
#' @return Ball Covariance statistic
#' @noRd
kbcov_stat <- function(x, num, var_num, weight, type) {
  compare_list <- list()
  prop_in_ball_vec <- c()
  value_diff <- numeric(1)
  value_diff_prob <- numeric(1)
  value_diff_hhg <- numeric(1)
  hhg_ball_num <- numeric(1)
  #
  stat_value <- numeric(1)
  stat_value_prob <- numeric(1)
  stat_value_hhg <- numeric(1)
  
  # stat_value_name <- c("bcov", "bcov.prob", "bcov.hhg")
  # weight_name <- c("none", "prob", "hhg")
  # names(stat_value) <- stat_value_name[which(weight_name == weight)]
  # if(type == "bcor") {
  #   names(stat_value) <- gsub(x = names(stat_value), 
  #                             pattern = "bcov", replacement = "bcor")
  # }
  
  # compute extention of BCov:
  for (i in 1:num) {
    for (j in 1:num) {
      value_diff_hhg <- 0
      value_diff_prob <- 0
      value_diff <- 0
      # compute ball statistic for ball with sample i as center and radius is d(x_{i}, x_{j})
      # d(x_{i}, x_{j}) are stored in x[[v]][i, j], x = 1, ..., var_num
      all_in_ball_vec <- rep(1, num)
      for (v in 1:var_num) {
        compare_list[[v]] <- (x[[v]][, i] <= x[[v]][i, j])
        all_in_ball_vec <- all_in_ball_vec * compare_list[[v]]
      }
      prop_in_ball_vec <- sapply(compare_list, mean)
      value_diff <- (mean(all_in_ball_vec) - prod(prop_in_ball_vec))^2
      value_diff_prob <- value_diff / prod(prop_in_ball_vec)
      if (!any(prop_in_ball_vec %in% c(0, 1))) {
        value_diff_hhg <- value_diff / (prod(prop_in_ball_vec)*(prod(1 - prop_in_ball_vec)))
        hhg_ball_num <- hhg_ball_num + 1
      }
      # aggrate statistic value:
      stat_value <- stat_value + value_diff
      stat_value_prob <- stat_value_prob + value_diff_prob
      stat_value_hhg <- stat_value_hhg + value_diff_hhg
    }
  }
  stat_value <- stat_value / (num)^2
  stat_value_prob <- stat_value_prob / (num)^2
  stat_value_hhg <- stat_value_hhg / (hhg_ball_num)^2
  # 
  c("bcov.constant" = stat_value, "bcov.probability" = stat_value_prob, "bcov.chisquare" = stat_value_hhg)
}


#' Wrap kbcov_test_internal and bcov_test_internal
#' @inheritParams bcov.test
#' @noRd
bcov_test_internal_wrap <- function(x = x, y = y, num.permutations, distance, seed, 
                                    weight, method, type, num.threads)
{
  if(class(x) == "list") {
    if (length(x) > 2)
    {
      result <- kbcov_test_internal(x = x, num.permutations = num.permutations, distance = distance, weight = weight, 
                                    seed = seed, method = method, num.threads = num.threads)
    } else {
      y <- x[[2]]
      x <- x[[1]]
      result <- bcov_test_internal(x = x, y = y, num.permutations = num.permutations, distance = distance, 
                                   weight = weight, seed = seed, method = method, 
                                   num.threads = num.threads)
    }
  } else {
    result <- bcov_test_internal(x = x, y = y, num.permutations = num.permutations, distance = distance, 
                                 weight = weight, seed = seed, method = method, 
                                 num.threads = num.threads)
  }
  result
}


#' @title Ball Covariance and Correlation Statistics
#' @description Computes Ball Covariance and Ball Correlation statistics, 
#' which are generic dependence measures in Banach spaces.
#' @inheritParams bcov.test
#' @rdname bcov
#' 
#' @details 
#' The sample sizes of the two variables must agree,  and samples must not contain missing and infinite values. 
#' If we set \code{distance = TRUE}, arguments \code{x}, \code{y} can be a \code{dist} object or a
#' symmetric numeric matrix recording distance between samples; otherwise, these arguments are treated as data.
#' 
#' \code{bcov} and \code{bcor} compute Ball Covariance and Ball Correlation statistics.
#' 
#' Ball Covariance statistics is a generic dependence measure in Banach spaces. It enjoys the following properties: 
#' \itemize{
#' \item It is nonnegative and it is equal to zero if and only if variables are unassociated; 
#' \item It is highly robust;   
#' \item It is distribution-free and model-free;   
#' \item it is interesting that the HHG is a special case of Ball Covariance statistics. 
#' }
#' Ball correlation statistics, a normalized version of Ball Covariance statistics, generalizes Pearson correlation in two fundamental ways: 
#' \itemize{
#' \item It is well-defined for random variables in arbitrary dimension in Banach spaces
#' \item BCor is equal to zero implies random variables are unassociated.
#' }
#' 
#' The definitions of the Ball Covariance and Ball Correlation statistics between two random variables are as follows.
#' Suppose, we are given pairs of independent observations 
#' \eqn{\{(x_1, y_1),...,(x_n,y_n)\}}, where \eqn{x_i} and \eqn{y_i} can be of any dimension 
#' and the dimensionality of \eqn{x_i} and \eqn{y_i} need not be the same.
#' Then, we define sample version Ball Covariance as:
#' \deqn{\mathbf{BCov}_{\omega, n}^{2}(X, Y)=\frac{1}{n^{2}}\sum_{i,j=1}^{n}{(\Delta_{ij,n}^{X,Y}-\Delta_{ij,n}^{X}\Delta_{ij,n}^{Y})^{2}} }
#' where:
#' \deqn{ \Delta_{ij,n}^{X,Y}=\frac{1}{n}\sum_{k=1}^{n}{\delta_{ij,k}^{X} \delta_{ij,k}^{Y}}, 
#' \Delta_{ij,n}^{X}=\frac{1}{n}\sum_{k=1}^{n}{\delta_{ij,k}^{X}}, 
#' \Delta_{ij,n}^{Y}=\frac{1}{n}\sum_{k=1}^{n}{\delta_{ij,k}^{Y}} }
#' \deqn{\delta_{ij,k}^{X} = I(x_{k} \in \bar{B}(x_{i}, \rho(x_{i}, x_{j}))), 
#' \delta_{ij,k}^{Y} = I(y_{k} \in \bar{B}(y_{i}, \rho(y_{i}, y_{j})))}
#' Among them, \eqn{\bar{B}(x_{i}, \rho(x_{i}, x_{j}))} is a closed ball 
#' with center \eqn{x_{i}} and radius \eqn{\rho(x_{i}, x_{j})}.
#' Similarly, we can define \eqn{ \mathbf{BCov}_{\omega,n}^2(\mathbf{X},\mathbf{X}) } 
#' and \eqn{ \mathbf{BCov}_{\omega,n}^2(\mathbf{Y},\mathbf{Y}) }. 
#' We define Ball Correlation statistic as follows.
#' \deqn{\mathbf{BCor}_{\omega,n}^2(\mathbf{X},\mathbf{Y})=
#' \mathbf{BCov}_{\omega,n}^2(\mathbf{X},\mathbf{Y})/\sqrt{\mathbf{BCov}_{\omega,n}^2(\mathbf{X},\mathbf{X})\mathbf{BCov}_{\omega,n}^2(\mathbf{Y},\mathbf{Y})}
#' }
#' 
#' We can extend \eqn{\mathbf{BCov}_{\omega,n}} to measure the mutual independence between \eqn{K} random variables:
#' \deqn{\frac{1}{n^{2}}\sum_{i,j=1}^{n}{\left[ (\Delta_{ij,n}^{X_{1}, ..., X_{K}}-\prod_{k=1}^{K}\Delta_{ij,n}^{X_{k}})^{2}\prod_{k=1}^{K}{\hat{\omega}_{k}(X_{ki},X_{kj})} \right]}}
#' where \eqn{X_{k}(k=1,\ldots,K)} are random variables and \eqn{X_{ki}} is the \eqn{i}-th observations of \eqn{X_{k}}. 
#' 
#' See \code{\link{bcov.test}} for a test of independence based on the Ball Covariance statistic.
#' 
#' @return 
#' \item{\code{bcov }}{ Ball Covariance statistic.}
#' @seealso
#' \code{\link{bcov.test}}, \code{\link{bcorsis}}
#' @export
#' 
#' @references Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019) Ball Covariance: A Generic Measure of Dependence in Banach Space, Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1543600
#' @references Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic Sure Independence Screening Procedure, Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
#' @references Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2018). Ball: An R package for detecting distribution difference and association in metric spaces. arXiv preprint arXiv:1811.03750. http://arxiv.org/abs/1811.03750
#' 
#' @examples
#' ############# Ball Covariance #############
#' num <- 50
#' x <- rnorm(num)
#' y <- rnorm(num)
#' bcov(x, y)
#' bcov(x, y, weight = "prob")
#' bcov(x, y, weight = "chisq")
bcov <- function(x, y, distance = FALSE, weight = FALSE) {
  weight <- examine_weight_arguments(weight)
  res <- bcov_test_internal_wrap(x = x, y = y, num.permutations = 0, distance = distance, seed = 1,
                                 weight = weight, method = "permute", num.threads = 0)
  res
}

