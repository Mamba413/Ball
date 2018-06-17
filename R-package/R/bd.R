#' @title Ball Divergence based Equality of Distributions Test
#' @description Performs the nonparametric two-sample or K-sample ball divergence test for
#' equality of multivariate distributions
#' @aliases bd.test
#' @author XueQin Wang, WenLiang Pan, HePing Zhang, Yuan Tian
#' @param x a numeric vector, matrix, data.frame, \code{dist} object or list contains vector, matrix or data.frame.
#' @param y a numeric vector, matrix or data.frame.
#' @param R the number of replications, when R equals to 0, the function returns
#' the sample version of ball divergence. Default: \code{R = 99}
#' @param dst if \code{dst = TRUE}, x will be considered as a distance matrix. Default: \code{dst = FALSE}
#' @param size a vector record sample size of each group.
#' @param seed the random seed. 
#' @param num.threads Number of threads. Default \code{num.threads = 2}.
#' @param kbd.type the type of K-sample test statistics. Setting \code{kbd.type = "sum"} print the statistic 
#' and \eqn{p}-value of the summation version of \eqn{K} sample ball divergence while setting \code{kbd.type = "max"} print
#' the maximum version. Further, you can obtain the information of both summation and maximum by \code{summary} function.
#' 
## @param weight not available now
## @param method if \code{method = 'permute'}, a permutation procedure will be carried out;
## if \code{ method = 'approx'}, the p-values based on approximate Ball Divergence
## distribution are given.
#' 
#' @return bd.test returns a list with class "htest" containing the following components:
#' \item{\code{statistic}}{ball divergence statistic.}            
#' \item{\code{p.value}}{the p-value for the test.}
#' \item{\code{replicates}}{replicates of the test statistic.}
#' \item{\code{size}}{sample sizes.}
#' \item{\code{complete.info}}{a \code{list} containing multiple statistics value and their corresponding $p$ value.}
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' 
#' @rdname bd.test
#' 
#' @details 
#' \code{bd.test} are ball divergence based multivariate nonparametric tests of two-sample or 
#' K-sample problem. If only \code{x} is given, the statistic is 
#' computed from the original pooled samples, stacked in 
#' matrix where each row is a multivariate observation, or from the distance matrix 
#' when \code{dst = TRUE}. The first \code{sizes[1]} rows of \code{x} are the first sample, the next 
#' \code{sizes[2]} rows of \code{x} are the second sample, etc.
#' If x is a list, its elements are taken as the samples to be compared, 
#' and hence have to be numeric data vectors, matrix or data.frame.
#' 
#' Based on sample version ball divergence (see \code{\link{bd}}), the test is implemented by 
#' permutation with R replicates. The function simply returns the test statistic 
#' when \code{R = 0}.
#' 
#' @seealso
#' \code{\link{bd}}
#' 
#' @references Pan, Wenliang; Tian, Yuan; Wang, Xueqin; Zhang, Heping. Ball Divergence: Nonparametric two sample test. Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. https://projecteuclid.org/euclid.aos/1525313077
#' 
#' @export
#' @examples
#' ################# Quick Start #################
#' x <- rnorm(50)
#' y <- rnorm(50, mean = 1)
#' # plot(density(x))
#' # lines(density(y), col = "red")
#' # ball divergence:
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
#' bd.test(x = Dmat, size = c(150, 150), R = 99, dst = TRUE)
#'
#' ################# Non-Hilbert Real Data #################
#' # load data:
#' data("macaques")
#' # number of femala and male Macaca fascicularis:
#' table(macaques[["group"]])
#' # calculate Riemannian shape distance matrix:
#' Dmat <- nhdist(macaques[["x"]], method = "riemann")
#' # hypothesis test with BD:
#' bd.test(x = Dmat, R = 99, size = c(9, 9), dst = TRUE)
#' 
#' ################  K-sample Test  #################
#' n <- 150
#' bd.test(rnorm(n), size = c(40, 50, 60))
#' # alternative input method:
#' x <- lapply(c(40, 50, 60), rnorm)
#' bd.test(x)
#' 
bd.test <- function(x, y = NULL, R = 99, dst = FALSE,
                    size = NULL, seed = 4, num.threads = 2, 
                    kbd.type = "sum") {
  weight = FALSE
  method = 'permute'
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  if(is.null(x)|is.null(y)) {
    # modify input information:
    data_name <- gsub(x = data_name, pattern = " and NULL", replacement = "")
    
    # examine input arguments x and y:
    if(is.list(x)) {
      x <- lapply(x, as.matrix)
      size <- sapply(x, nrow)
      x <- do.call("rbind", x)
    } else {
      x <- get_matrixed_x(x, y)
    }
    # examine input arguments size:
    examine_size_arguments(x, size)
    # 
    if(dst) {
      xy <- as.vector(x)
    } else {
      p <- ncol(x)
      if(p > 1) {
        xy <- as.vector(as.matrix(dist(x, diag = TRUE)))
        dst <- TRUE
      } else {
        xy <- x
        dst <- FALSE
      }
    }
  } else {
    x <- as.matrix(x)
    y <- as.matrix(y)
    # examine dimension:
    p <- examine_dimension(x, y)
    # 
    if(p > 1) {
      xy <- get_vectorized_distance_matrix(x, y)
      dst <- TRUE
      size <- c(xy[[2]], xy[[3]])
      xy <- xy[[1]]
    } else {
      xy <- rbind(x, y)
      dst <- FALSE
      size <- c(dim(x)[1], dim(y)[1])
    }
  }
  ## memory protect step:
  memoryAvailable(n = sum(size), funs = 'BD.test')
  
  ## examine R arguments:
  if(method == "approx") {
    R <- 0
  } else {
    examine_R_arguments(R)
  }
  
  ## examine num.thread arguments:
  examine_threads_arguments(num.threads)
  
  ## main:
  if(R == 0) {
    result <- bd_test_wrap_c(xy, size, R = 0, weight, dst, num.threads)
    # approximately method:
    if(method == "approx") {
      if(result[["info"]][["K"]] == 2) {
        pvalue <- calculatePvalue(prod(size) * result[["statistic"]] / sum(size), 
                                  BDTestNullDistribution)
      } else {
        return(result[["statistic"]])
      }
    } 
    # return statistic when R = 0:
    else {
      if (result[["info"]][["K"]] == 2) {
        if (weight) {
          return_stat <- result[["statistic"]][2]
        } else {
          return_stat <- result[["statistic"]][1] 
        }
      } else {
        if (kbd.type == "sum") {
          return_stat <- result[["statistic"]][1]
        } else {
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
    result <- bd_test_wrap_c(xy, size, R, weight, dst, num.threads)
    # pvalue <- calculatePvalue(result[["statistic"]], result[["permuted_stat"]])
  }
  # output information:
  if (result[["info"]][["K"]] == 2) {
    stat <- result[["statistic"]][1]
    pvalue <- result[["p.value"]][1]
    stat_message <- ""
  } else if (kbd.type == "sum") {
    stat <- result[["statistic"]][1]
    pvalue <- result[["p.value"]][1]
    stat_message <- " (Summation Version)"
  } else {
    stat <- result[["statistic"]][3]
    pvalue <- result[["p.value"]][3]
    stat_message <- " (Maximum Version)"
  }
  data_name <- paste(data_name, sprintf("\nnumber of observations = %s,", result[["info"]][["N"]]))
  data_name <- paste(data_name, "group sizes:", paste0(result[["info"]][["size"]], collapse = " "))
  data_name <- paste0(data_name, "\nreplicates = ", R)
  # data_name <- paste0(data_name, ", Weighted Ball Divergence = ", result[["info"]][["weight"]])
  alternative_message <- "distributions of samples are distinct"
  
  # return:
  e <- list(
    statistic = stat,
    p.value = pvalue,
    replicates = R,
    size = result[["info"]][["size"]],
    complete.info = result,
    alternative = alternative_message,
    method = sprintf("Nonparametric %s-Samples Ball Divergence Test%s", result[["info"]][["K"]], stat_message),
    data.name = data_name
  )
  class(e) <- "htest"
  return(e)
}


#' @title Ball Divergence
#' @description Compute ball divergence statistic between two-sample or K-sample.
#' @author XueQin Wang, WenLiang Pan, HePing Zhang, Yuan Tian
#' @inheritParams bd.test
#' @rdname bd
#' @return 
#' \item{\code{bd }}{ sample version of ball divergence}
#' 
#' @details 
#' Given the samples not containing missing values, \code{bd} returns sample version of ball divergence.
#' If we set \code{dst = TRUE}, arguments \code{x}, \code{y} can be a \code{dist} object or a
#' symmetric numeric matrix recording distance between samples; 
#' otherwise, these arguments are treated as data.
#' 
#' Ball divergence, introduced by Pan et al(2017), is a new concept to measure the difference 
#' between two probability distributions in separable Banach space. 
#' Ball divergence of two probability measures is proven to be zero if and only if they are identical.
#' 
#' The definitions of the sample version ball divergence are as follows.
#' Given two independent samples \eqn{ \{x_{1}, ..., x_{n}\} } with the associated probability measure \eqn{\mu} and 
#' \eqn{ \{y_{1}, ..., y_{m}\} } with \eqn{\nu}, where the observations in each sample are \emph{i.i.d}.
#' 
#' Also, let \eqn{\delta(x,y,z)=I(z\in \bar{B}(x, \rho(x,y)))}, 
#' where \eqn{\delta(x,y,z)} indicates whether \eqn{z} is located in the closed ball \eqn{\bar{B}(x, \rho(x,y))} 
#' with center \eqn{x} and radius \eqn{\rho(x, y)}. 
#' We denote:
#' \deqn{
#' A_{ij}^{X}=\frac{1}{n}\sum_{u=1}^{n}{\delta(X_i,X_j,X_u)}, \quad A_{ij}^{Y}=\frac{1}{m}\sum_{v=1}^{m}{\delta(X_i,X_j,Y_v)}
#' }
#' \deqn{
#' C_{kl}^{X}=\frac{1}{n}\sum_{u=1}^{n}{\delta(Y_k,Y_l,X_u)}, \quad C_{kl}^{Y}=\frac{1}{m}\sum_{v=1}^{m}{\delta(Y_k,Y_l,Y_v)}
#' }
#' 
#' \eqn{A_{ij}^X} represents the proportion of samples \eqn{ \{x_{1}, ..., x_{n}\} } located in the 
#' ball \eqn{\bar{B}(X_i,\rho(X_i,X_j))} and \eqn{A_{ij}^Y} represents the proportion of samples \eqn{ \{y_{1}, ..., y_{m}\} } 
#' located in the ball \eqn{\bar{B}(X_i,\rho(X_i,X_j))}. 
#' Meanwhile, \eqn{C_{kl}^X} and \eqn{C_{kl}^Y} 
#' represent the corresponding proportions located in the ball \eqn{\bar{B}(Y_k,\rho(Y_k,Y_l))}.
#' 
#' we can define sample version ball divergence as:
#' \deqn{D_{n,m}=A_{n,m}+C_{n,m}}
#' 
#' BD can be generalized to the \emph{K}-sample problem, i.e. if we 
#' have \eqn{K} group samples, each group include \eqn{n^{(k)}, k=1,...,K} samples, 
#' then we can define sample version of generalized ball divergence for \emph{K}-sample problem:
#' \deqn{\sum_{1 \leq k < l \leq K}{D_{n^{(k)},n^{(l)}}}}
#' 
#' See \code{\link{bd.test}} for a test of multivariate independence based on the 
#' ball divergence.
#' 
#' @seealso
#' \code{\link{bd.test}}
#' @export
#' 
#' @references Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. (2017) Ball divergence: nonparametric two sample test, \emph{The Annals of Statistics}, to appear
#' 
#' @examples
#' ############# Ball Divergence #############
#' x <- rnorm(50)
#' y <- rnorm(50)
#' bd(x, y)
bd <- function(x, y = NULL, dst = FALSE, size = NULL, num.threads = 2, kbd.type = "sum") {
  res <- bd.test(x = x, y = y, dst = dst, size = size, R = 0)
  res
}

