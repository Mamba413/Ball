#' @title Ball Covariance Test
#' @author WenLiang Pan, XueQin Wang, HePing Zhang, HongTu Zhu, Jin Zhu
#' @description Ball covariance test of multivariate independence. 
#' Ball covariance are generic multivariate measures of dependence in banach space.
#' 
#' @inheritParams bd.test
#' @param x a numeric vector, matirx, data.frame or \code{dist} object or list contains numeric vector, matrix or data.frame.
#' @param y a numeric vector, matirx, data.frame or \code{dist} object.
#' @param R the number of replication. When \code{R = 0}, the function return ball covariance. Default: \code{R = 99}
#' @param dst if \code{dst = TRUE}, \code{x} and \code{y} will be considered as distance matrix. Default: \code{dst = FALSE}
#' @param weight when \code{weight = TRUE}, weighted ball covariance or weighted ball correlation is used instead of ball covariance
#' or ball correlation. Default: \code{weight = FALSE} 
## @param type If \code{type = 'bcor'}, ball correlation will be used instead of ball covariance.(default \code{type = 'bcov'})
## @param method if \code{method = 'permute'}, a permutation procedure will be carried out;
## if \code{method = 'approx'}, the p-values based on approximate Ball Covariance distribution are given.(Test arguments)
#' 
#' @return bcov.test returns a list with class "htest" containing the following components:
#' \item{\code{statistic}}{ball covariance or ball correlation statistic.}            
#' \item{\code{permuted_stat}}{permuted ball covariance or ball correlation statistic.} 
#' \item{\code{N}}{sample size.} 
#' \item{\code{p.value}}{the p-value for the test.}  
#' \item{\code{alternative}}{a character string describing the alternative hypothesis.}
#' \item{\code{method}}{a character string indicating what type of test was performed.}
#' \item{\code{data.name}}{description of data.}
#' 
#' @details 
#' \code{bcov.test} are non-parametric tests of multivariate independence in banach space. 
#' The test decision is obtained via permutation, with \code{R} replicates.
#' 
#' If two samples are pass to arguments \code{x} and \code{y}, the sample sizes (i.e. number of rows or length of vector) 
#' of the two variables must agree. If a \code{\link{list}} object is pass to \code{x}, 
#' each element must with same sample sizes. Moreover, data pass to \code{x} or \code{y} 
#' must not contain missing or infinite values. 
#' If we set \code{dst = TRUE}, arguments \code{x}, \code{y} can be \code{dist} object or a
#' symmetric numeric matrix recording distance between samples; 
#' otherwise these arguments are treated as data.
#' 
#' The \code{bcov.test} statistic is \code{bcov} or \code{ bcor} which are dependence measure 
#' in banach space. The \code{bcor} test statistic is based on the normalized 
#' coefficient of ball covariance. (See the manual page for \code{\link{bcov}} or \code{\link{bcor}}.)
#' 
#' For the general problem of testing independence when the distributions of \eqn{X} and 
#' \eqn{Y} are unknown, the test based on \code{bcov} can be implemented as a permutation test.
#' See (Pan et al 2017) for theoretical properties of the test, including statistical consistency.
#' 
#' 
#' @useDynLib Ball, .registration = TRUE
#' @export
#' @seealso \code{\link{bcov}}, \code{\link{bcor}}
#' @examples
#' 
#' ################# Quick Start #################
#' error <- runif(50, min = -0.3, max = 0.3)
#' x <- runif(50, 0, 4*pi)
#' y <- cos(x) + error
#' # plot(x, y)
#' bcov.test(x, y)
#' 
#' ################# Quick Start #################
#' x <- matrix(runif(50 * 2, -pi, pi), nrow = 50, ncol = 2)
#' error <- runif(50, min = -0.3, max = 0.3)
#' y <- (sin((x[,1])^2 + x[,2])) + error
#' bcov.test(x = x, y = y)
#' 
#' ################# Ball Covariance Test for Non-Hilbert Data #################
#' # load data:
#' data("ArcticLake")
#' # Distance matrix between y:
#' Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
#' # Distance matrix between x:
#' Dx <- dist(ArcticLake[["depth"]])
#' # hypothesis test with BCov:
#' bcov.test(x = Dx, y = Dy, R = 99, dst = TRUE)
#' 
#' ################  Weighted Ball Covariance Test  #################
#' data("ArcticLake")
#' Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
#' Dx <- dist(ArcticLake[["depth"]])
#' # hypothesis test with weighted BCov:
#' bcov.test(x = Dx, y = Dy, R = 99, dst = TRUE, weight = TRUE)
#' 
#' ################# Mutual Independence Test #################
#' x <- rnorm(30)
#' y <- (x > 0) * x + rnorm(30)
#' z <- (x <= 0) * x + rnorm(30)
#' data_list <- list(x, y, z)
#' bcov.test(data_list)
#' 
bcov.test <- function(x, y = NULL, R = 99, dst = FALSE, weight = FALSE, 
                      seed = 4)
{
  method = 'permute'
  data_name <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
  type <- "bcov"
  # modify input information:
  if(class(x) == "list") {
    data_name <- gsub(x = data_name, pattern = " and NULL", replacement = "")
  }
  result <- bcov_test_internal_wrap(x = x, y = y, R = R, dst = dst, weight = weight, 
                                    seed = seed, method = method, type = type)
  # return result of hypothesis test:
  data_name <- paste0(data_name,"\nnumber of observations = ", result[["N"]])
  data_name <- paste0(data_name, "\nreplicates = ", R, 
                      ", Weighted Ball Covariance = ", weight)
  test_method <- "Ball Covariance test of independence"
  if(type == "bcor") {
    test_method <- gsub(pattern = "Covariance", replacement = "Correlation", x = test_method)
    data_name <- gsub(pattern = "Covariance", replacement = "Correlation", x = data_name)
  }
  alternative_message <- "random variables are dependent"
  e <- list(
    statistic = result[["statistic"]],
    permuted_stat = result[["permuted_stat"]],
    N = result[["N"]],
    p.value = result[["p.value"]],
    alternative = alternative_message,
    method = test_method,
    data.name = data_name
  )
  class(e) <- "htest"
  return(e)
}


#' Ball covariance test internal function
#' @inheritParams bcov.test
#' @param type 
#'
#' @noRd
bcov_test_internal <- function(x, y, R = 99, dst = FALSE, weight = FALSE, 
                               seed = 4, method = 'permute', type = 'bcov')
{
  x <- as.matrix(x)
  y <- as.matrix(y)
  x_y_info <- examine_x_y(x, y)
  num <- x_y_info[1]
  p <- x_y_info[2]
  #
  if(dst == FALSE) {
    if(p != 1) {
      x <- as.vector(as.matrix(dist(x, diag = TRUE)))
      y <- as.vector(as.matrix(dist(y, diag = TRUE)))
      dst <- TRUE
    }
  } else {
    x <- as.vector(x)
    y <- as.vector(y)
  }
  ## memory protect step:
  memoryAvailable(num, funs = 'BI.test')
  ## examine test type:
  type <- examine_type_arguments(type)
  ## examine R arguments:
  if(method == "approx") {
    R <- 0
  } else {
    examine_R_arguments(R)
  }
  #
  if(R == 0) {
    result <- bcov_value_wrap_c(x = x, y = y, n = num, weight = weight, dst =  dst, type = type)
    if(method == "approx") {
      pvalue <- calculatePvalue(result[["statistic"]] * result[["info"]][["N"]], 
                                BITestNullDistribution)
    } else {
      return(result[["statistic"]])
    }
  } else {
    set.seed(seed = examine_seed_arguments(seed))
    result <- bcov_test_wrap_c(x = x, y = y, n = num, R = R, weight = weight, dst =  dst, type = type)
    pvalue <- calculatePvalue(result[["statistic"]], result[["permuted_stat"]])
  }
  list(statistic = result[["statistic"]],
       permuted_stat = result[["permuted_stat"]],
       N = num,
       p.value = pvalue)
}


#' A internal function for carry out independence test for multiple random variables
#' @inheritParams bcov.test
#' @inherit return
#' @noRd
#'
kbcov_test_internal <- function(x, R = 99, dst = FALSE, weight = FALSE, 
                                seed = 4, method = 'permute', type = 'bcov')
{
  x <- lapply(x, as.matrix)
  size_list <- sapply(x, nrow)
  num <- unique(size_list)
  if(length(num) > 1) {
    stop("sample sizes in each list must equal!")
  }
  #
  if(dst) {
    
  } else {
    x <- lapply(x, dist, diag = TRUE, upper = TRUE)
    x <- lapply(x, as.matrix)
  }
  var_num <- length(x)
  # compute statistic:
  stat_value <- kbcov_stat(x = x, num = num, var_num = var_num, 
                           weight = weight, type = type)
  if(R == 0) {
    names(stat_value)
    return(stat_value)
  } else {
    seed <- examine_seed_arguments(seed)
    set.seed(seed)
    # permutation procedure:
    permuted_stat <- numeric(R)
    for (r in 1:R) {
      x_copy <- x
      for (v in 1:var_num) {
        index <- sample(1:num, size = num, replace = FALSE)
        x_copy[[v]] <- x[[v]][index, index]
      }
      permuted_stat[r] <- kbcov_stat(x = x_copy, num = num, 
                                     var_num = var_num, 
                                     weight = weight, type = type)
    }
    # calculate pvalue:
    pvalue <- calculatePvalue(stat_value, permuted_stat)
  }
  # return result:
  list("statistic" = stat_value, 
       "permuted_stat" = permuted_stat, 
       "N" = num, 
       "p.value" = pvalue)
}


#' compute extension of BCov for independence of multiple random variables
#'
#' @param x list containing distance matrix, each element is distance matrix
#' @param num sample size
#' @param var_num random variables number
#' @param weight whether used weight
#' @param type ball correlation or ball covariance. now, only ball covariance is considered.
#'
#' @return ball covariance statistic
#' @noRd
#'
kbcov_stat <- function(x, num, var_num, weight, type) {
  compare_list <- list()
  prop_in_ball_vec <- c()
  value_diff <- numeric(1)
  #
  stat_value <- numeric(1)
  names(stat_value) <- ifelse(weight, "wbcov", "bcov")
  if(type == "bcor") {
    names(stat_value) <- gsub(x = names(stat_value), 
                              pattern = "bcov", replacement = "bcor")
  }
  # compute extention of BCov:
  for (i in 1:num) {
    for (j in 1:num) {
      # compute ball statistic for ball with sample i as center and radius is d(x_{i}, x_{j})
      # d(x_{i}, x_{j}) are stored in x[[v]][i, j], x = 1, ..., var_num
      all_in_ball_vec <- rep(1, num)
      for (v in 1:var_num) {
        compare_list[[v]] <- (x[[v]][, i] <= x[[v]][i, j])
        all_in_ball_vec <- all_in_ball_vec * compare_list[[v]]
      }
      prop_in_ball_vec <- sapply(compare_list, mean)
      value_diff <- (mean(all_in_ball_vec) - prod(prop_in_ball_vec))^2
      if(weight) {
        value_diff <- value_diff / prod(prop_in_ball_vec)
      }
      # aggrate statistic value:
      stat_value = stat_value + value_diff
    }
  }
  # 
  stat_value
}


#' Wrap kbcov_test_internal and bcov_test_internal
#' @inheritParams bcov.test
#' @noRd
#'
bcov_test_internal_wrap <- function(x = x, y = y, R, dst, seed, 
                                    weight, method, type)
{
  if(class(x) == "list") {
    result <- kbcov_test_internal(x = x, R = R, dst = dst, weight = weight, 
                                  seed = seed, 
                                  method = method, type = type)
  } else {
    result <- bcov_test_internal(x = x, y = y, R = R, dst = dst, 
                                 weight = weight, seed = seed, method = method, 
                                 type = type)
  }
  result
}


#' @title Ball Correlation and Covariance Statistics
#' @description Computes ball covariance and ball correlation statistics, 
#' which are multivariate measures of dependence in banach space.
#' @inheritParams bcov.test
#' @rdname bcov
#' 
#' @details 
#' \code{bcov} and \code{bcor} compute ball covariance and ball correlation statistics.
#' 
#' The sample sizes(number of rows or length of vector) of the two variables must agree, 
#' and samples must not contain missing values. 
#' If we set \code{dst = TRUE}, arguments \code{x}, \code{y} can be \code{dist} object or a
#' symmetric numeric matrix recording distance between samples; 
#' otherwise these arguments are treated as data.
#' 
#' Ball covariance is a generic non-parametric dependence measure in banach space, introduced by Pan et al(2017). 
#' It is noteworthy that ball covariance enjoys the following properties: 
#' 
#' (i) It is nonnegative, and holds the Cauchy-Schwartz type inequality; 
#' 
#' (ii) It is nonparametric and makes fewer restrictive data assumptions even without finite moment conditions;  
#' 
#' (iii) Its empirical version is feasible and can be used as a test statistic of independence with some desired test properties;  
#' 
#' (iv) it is interesting that the HHG dependence measure is a special case of ball covariance. 
#' 
#' Ball correlation, based on the normalized ball covariance, generalizes the idea of pearson correlation in two fundamental ways: 
#' 
#' (i) Ball correlation, \eqn{ \mathbf{BCor}_{\omega}^{2}(X, Y) }, is defined for \eqn{X} and \eqn{Y} in arbitrary dimension in banach space. 
#' 
#' (ii) Ball correlation satisfies \eqn{0 \le \mathbf{BCor}_{\omega}^{2}(X, Y) \le 1}, and \eqn{ \mathbf{BCor}_{\omega}^{2}(X, Y) } = 0 
#' only if \eqn{X} and \eqn{Y} are independent.
#' 
#' The definitions of the sample version ball covariance and ball correlation are as follows.
#' Suppose, We are given pairs of independent observations 
#' \eqn{\{(x_1, y_1),...,(x_n,y_n)\}}, where \eqn{x_i} and \eqn{y_i} can be of any dimension 
#' and the dimensionality of \eqn{x_i} and \eqn{y_i} need not be the same.
#' Then, we define sample version ball covariance as:
#' \deqn{\mathbf{BCor}_{\omega, n}^{2}(X, Y)=\frac{1}{n^{2}}\sum_{i,j=1}^{n}{(\Delta_{ij,n}^{X,Y}-\Delta_{ij,n}^{X}\Delta_{ij,n}^{Y})^{2}} }
#' where:
#' \deqn{ \Delta_{ij,n}^{X,Y}=\frac{1}{n}\sum_{k=1}^{n}{\delta_{ij,k}^{X} \delta_{ij,k}^{Y}}, 
#' \Delta_{ij,n}^{X}=\frac{1}{n}\sum_{k=1}^{n}{\delta_{ij,k}^{X}}, 
#' \Delta_{ij,n}^{Y}=\frac{1}{n}\sum_{k=1}^{n}{\delta_{ij,k}^{Y}} }
#' \deqn{\delta_{ij,k}^{X} = I(x_{k} \in \bar{B}(x_{i}, \rho(x_{i}, x_{j}))), 
#' \delta_{ij,k}^{Y} = I(y_{k} \in \bar{B}(y_{i}, \rho(y_{i}, y_{j})))}
#' 
#' Among them, \eqn{\bar{B}(x_{i}, \rho(x_{i}, x_{j}))} is a closed ball 
#' with center \eqn{x_{i}} and radius \eqn{\rho(x_{i}, x_{j})}.
#' Similarly, we can give the notations \eqn{ \mathbf{BCov}_{\omega,n}^2(\mathbf{X},\mathbf{X}) } 
#' and \eqn{ \mathbf{BCov}_{\omega,n}^2(\mathbf{Y},\mathbf{Y}) }, 
#' which are the sample version version of \eqn{ \mathbf{BCov}_{\omega}^2(\mathbf{X},\mathbf{X}) } and 
#' \eqn{ \mathbf{BCov}_{\omega}^2(\mathbf{Y},\mathbf{Y}) }. 
#' We thus define the sample version ball correlation as follows.
#' 
#' \deqn{\mathbf{BCor}_{\omega,n}^2(\mathbf{X},\mathbf{Y})=
#' \mathbf{BCov}_{\omega,n}^2(\mathbf{X},\mathbf{Y})/\sqrt{\mathbf{BCov}_{\omega,n}^2(\mathbf{X},\mathbf{X})\mathbf{BCov}_{\omega,n}^2(\mathbf{Y},\mathbf{Y})}
#' }
#' 
#' Moreover, it is natural to extend \eqn{\mathbf{BCov}_{\omega,n}} to measure mutual independence between random vectors. For example, we can redefine
#' 
#' \deqn{\frac{1}{n^2}\sum_{i,j=1}^n(\Delta_{ij,n}^{XYZ}-\Delta_{ij,n}^X\Delta_{ij,n}^Y\Delta_{ij,n}^Z)^2\hat{\omega}_1(X_i,X_j)\hat{\omega}_2(Y_i,Y_j)\hat{\omega}_2(Z_i,Z_j)}
#' 
#' to measure whether \eqn{P_{XYZ}=P_X P_Y P_Z} while the computational complexity of this
#' extension remains the same. Similarily, we extend it to \eqn{K} random variables scenario:
#' 
#' \deqn{\frac{1}{n^{2}}\sum_{i,j=1}^{n}{(\Delta_{ij,n}^{R_{1}, ..., R_{K}}-\prod_{k=1}^{K}\Delta_{ij,n}^{R_{k}})^{2}\prod_{k=1}^{K}{\hat{\omega}_{k}(R_{ki},R_{kj})}}}
#' 
#' where \eqn{R_{k}, k=1,...K} indicate random variables and \eqn{R_{ki}, i=1,...,n} denote \eqn{n} random samples of \eqn{R_{k}}. 
#' 
#' See \code{\link{bcov.test}} for a test of multivariate independence based on the 
#' ball covariance and ball correlation statistic.
#' 
#' @return 
#' \item{\code{bcov }}{ sample version of ball covariance.}
#' @seealso
#' \code{\link{bcov.test}}, \code{\link{bcorsis}}
#' @export
#' 
#' @examples
#' ############# Ball Covariance #############
#' n <- 50
#' x <- rnorm(n)
#' y <- rnorm(n)
#' bcov(x, y)
bcov <- function(x, y, dst = FALSE, weight = FALSE) {
  res <- bcov_test_internal_wrap(x = x, y = y, R = 0, dst = dst, seed = 0,
                                 weight = weight, method = "permute", type = "bcov")
  res
}


#' @inheritParams bcov.test
#' @rdname bcov
#' @return 
#' \item{\code{bcor }}{ sample version of ball correlation.}
#' @export
#' @examples
#' ############# Ball Correlation #############
#' n <- 50
#' x <- rnorm(n)
#' y <- rnorm(n)
#' bcor(x, y)
bcor <- function(x, y, dst = FALSE, weight = FALSE) {
  res <- bcov_test_internal_wrap(x = x, y = y, R = 0, dst = dst, seed = 0,
                                 weight = weight, method = "permute", type = "bcor")
  res
}