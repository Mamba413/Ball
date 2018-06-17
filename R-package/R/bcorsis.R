#' @inheritParams bcov.test
#' @rdname bcov
#' @return 
#' \item{\code{bcor }}{ sample version of ball correlation.}
#' @export
#' @examples
#' ############# Ball Correlation #############
#' num <- 50
#' x <- 1:num
#' y <- 1:num
#' bcor(x, y)
#' bcor(x, y, weight = TRUE)
#' bcor(x, y, weight = "prob")
#' bcor(x, y, weight = "chisq")
bcor <- function(x, y, dst = FALSE, weight = FALSE) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  x_y_info <- examine_x_y(x, y)
  p <- x_y_info[2]
  #
  weight <- examine_weight_arguments(weight)
  #
  if(dst == FALSE) {
    if(p != 1) {
      x <- as.double(as.vector(as.matrix(dist(x, diag = TRUE))))
      y <- as.double(as.vector(as.matrix(dist(y, diag = TRUE))))
      dst_y <- as.integer(1)
      dst_x <- as.integer(1)
    } else {
      x <- as.double(x)
      y <- as.double(y)
      dst_y <- as.integer(0)
      dst_x <- as.integer(0)
    }
  } else {
    x <- as.double(as.vector(x))
    y <- as.double(as.vector(y))
    dst_y <- as.integer(1)
    dst_x <- as.integer(1)
  }
  #
  bcor_stat <- as.double(numeric(3))
  x_number <- as.integer(1)
  f_number <- as.integer(1)
  size_number <- as.integer(1)
  num <- as.integer(x_y_info[1])
  p <- as.integer(1)
  k <- as.integer(1)
  nth <- as.integer(1)
  #
  res <- .C("bcor_test", bcor_stat, y, x, x_number, f_number, size_number, num, p, k, dst_y, dst_x, nth)
  bcor_stat <- res[[1]]
  bcor_stat <- select_ball_stat(bcor_stat, weight, type = "bcor")
  return(bcor_stat)
}


#' @title Ball Correlation Sure Independence Screening
#' @author WenLiang Pan, WeiNan Xiao, XueQin Wang, HePing Zhang, HongTu Zhu
#' @description Generic non-parametric sure independence screening procedure based on ball correlation.
#' Ball correlation is a generic multivariate measure of dependence in Banach space.
#' @inheritParams bcov.test
#' @param x a numeric matirx or data.frame included \eqn{n} rows and \eqn{p} columns. 
#' Each row is an observation vector and each column corresponding to a explanatory variable, generally \eqn{p >> n}.
#' @param d the hard cutoff rule suggests selecting \eqn{d} variables. Setting \code{d = "large"} or 
#' \code{ d = "small"} means \code{n-1} or \code{floor(n/log(n))} 
#' variables are selected. If \code{d} is a integer, 
#' \code{d} variables are selected. Default: \code{d = "small"}
#' @param weight when \code{weight = TRUE}, weighted ball correlation is used instead of ball correlation. Default: \code{ weight = FALSE}  
#' @param method method for sure independence screening procedure, include: \code{"standard"},
#' \code{"lm"}, \code{"gam"}, \code{"interaction"} and \code{"survival"}.
#' Setting \code{method = "standard"} means standard sure independence screening procedure 
#' based on ball correlation while options
#' \code{"lm"} and \code{"gam"} carry out iterative BCor-SIS procedure with ordinary 
#' linear regression and generalized additive models, respectively.
#' Options \code{"interaction"} and \code{"survival"} are designed for detecting variables 
#' with potential linear interaction or associated with censored responses. Default: \code{method = "standard"}
#' @param dst if \code{dst = TRUE}, \code{y} will be considered as a distance matrix. 
#' Arguments only available when \code{ method = "standard"} and \code{ method = "interaction"}. Default: \code{dst = FALSE}
#' @param parms parameters list only available when \code{method = "lm"} or \code{"gam"}. 
#' It contains three parameters: \code{d1}, \code{d2}, and \code{df}. \code{d1} is the
#' number of initially selected variables, \code{d2} is the number of variables collection size added in each iteration.
#' \code{df} is degree freedom of basis in generalized additive models 
#' playing a role only when \code{method = "gam"}. Default: \code{ parms = list(d1 = 5, d2 = 5, df = 3)}
#' 
#' @return 
#' \item{\code{ix }}{ the vector of indices selected by ball correlation sure independence screening procedure.} 
#' \item{\code{method }}{ the method used.} 
#' \item{\code{weight }}{ the weight used.} 
#' \item{\code{complete.info }}{ a \code{list} containing at least one \eqn{p x 3} matrix, where each row is corresponding to variable and each column is corresponding to differe ball correlation weight. If \code{method = "gam"} or \code{method = "lm"}, \code{complete.info} is empty list.} 
#' 
#' @details 
#' \code{bcorsis} implements a model-free generic screening procedure, 
#' BCor-SIS, with fewer and less restrictive assumptions. 
#' The sample sizes (number of rows or length of the vector) of the 
#' two variables \code{x} and \code{y} must agree, 
#' and samples must not contain missing values. 
#' 
#' BCor-SIS procedure for censored response is carried out when \code{method = "survival"}. At that time, 
#' the matrix or data.frame pass to argument \code{y} must have exactly two columns and the first column is 
#' event (failure) time while the second column is censored status, a dichotomous variable. 
#' 
#' If we set \code{dst = TRUE}, arguments \code{y} is considered as distance matrix, 
#' otherwise \code{y} is treated as data.
#' 
#' BCor-SIS is based on a recently developed universal dependence measure: Ball correlation (BCor). 
#' BCor efficiently measures the dependence between two random vectors, which is between 
#' 0 and 1, and 0 if and only if these two random vectors are independent under some mild conditions.
#' (See the manual page for \code{\link{bcor}}.)
#' 
#' Theory and numerical result indicate that BCor-SIS has following advantages:
#' 
#' (i) It has a strong screening consistency property without finite sub-exponential moments of the data.
#' Consequently, even when the dimensionality is an exponential order of the sample size, BCor-SIS still 
#' almost surely able to retain the efficient variables.
#' 
#' (ii) It is nonparametric and has the property of robustness.
#' 
#' (iii) It works well for complex responses and/or predictors, such as shape or survival data
#' 
#' (iv) It can extract important features even when the underlying model is complicated.
#' 
#'   
#' @seealso 
#' \code{\link{bcor}}
#' 
#' @references Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic Sure Independence Screening Procedure, Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
#' 
#' @export
#' @examples 
#' \dontrun{
#' 
#' ############### Quick Start for bcorsis function ###############
#' set.seed(1)
#' n <- 150
#' p <- 3000
#' x <- matrix(rnorm(n * p), nrow = n)
#' error <- rnorm(n)
#' y <- 3*x[, 1] + 5*(x[, 3])^2 + error
#' res <- bcorsis(y = y, x = x)
#' head(res[["ix"]])
#' 
#' ############### BCor-SIS: Censored Data Example ###############
#' data("genlung")
#' result <- bcorsis(x = genlung[["covariate"]], y = genlung[["survival"]], 
#'                   method = "survival")
#' index <- result[["ix"]]
#' top_gene <- colnames(genlung[["covariate"]])[index]
#' head(top_gene, n = 1)
#' 
#' 
#' ############### BCor-SIS: Interaction Pursuing ###############
#' set.seed(1)
#' n <- 150
#' p <- 3000
#' x <- matrix(rnorm(n * p), nrow = n)
#' error <- rnorm(n)
#' y <- 3*x[, 1]*x[, 5]*x[, 10] + error
#' res <- bcorsis(y = y, x = x, method = "interaction")
#' head(res[["ix"]])
#' 
#' ############### BCor-SIS: Iterative Method ###############
#' library(mvtnorm)
#' set.seed(1)
#' n <- 150
#' p <- 3000
#' sigma_mat <- matrix(0.5, nrow = p, ncol = p)
#' diag(sigma_mat) <- 1
#' x <- rmvnorm(n = n, sigma = sigma_mat)
#' error <- rnorm(n)
#' rm(sigma_mat); gc(reset = TRUE)
#' y <- 3*(x[, 1])^2 + 5*(x[, 2])^2 + 5*x[, 8] - 8*x[, 16] + error
#' res <- bcorsis(y = y, x = x, method = "lm", d = 15)
#' res <- bcorsis(y = y, x = x, method = "gam", d = 15)
#' res[["ix"]]
#' 
#' ############### Weighted BCor-SIS: Probability weight ###############
#' set.seed(1)
#' n <- 150
#' p <- 3000
#' x <- matrix(rnorm(n * p), nrow = n)
#' error <- rnorm(n)
#' y <- 3*x[, 1] + 5*(x[, 3])^2 + error
#' res <- bcorsis(y = y, x = x, weight = "prob")
#' head(res[["ix"]])
#' # Alternative, chisq weight:
#' res <- bcorsis(y = y, x = x, weight = "chisq")
#' head(res[["ix"]])
#' }
bcorsis <- function(x, y, d = "small", weight = FALSE, 
                    method = "standard", dst = FALSE,
                    parms = list(d1 = 5, d2 = 5, df = 3),
                    num.threads = 2)
{
  seed <- 4
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- examine_x_y(x, y)[1]
  p <- dim(x)[2]
  y_p <- dim(y)[2]
  colnames(x) <- paste0("x", 1:p)
  colnames(y) <- paste0("y", 1:y_p)
  ids <- 1:p
  
  complete_info <- list()
    
  # check weight
  weight <- examine_weight_arguments(weight)
  
  # decide candicate size
  final_d <- examine_candiate_size(n, d, p)
  
  # get arguments:
  d1 <- parms$d1
  d2 <- parms$d2
  df <- parms$df
  
  # examine method arguments
  method <- examine_method_arguments(method)
  # examine dst and method arguments
  examine_dst_method(dst = dst, method = method)
  
  if(method == "survival") {
    Xhavepickout <- bcorsis.surv(y = y, x = x, final_d = final_d, 
                                 n = n, p = p, ids = ids)
    complete_info[[1]] <- Xhavepickout[[2]]
    Xhavepickout <- Xhavepickout[[1]]
  }
  if(method %in% c("standard", "pvalue", "interaction")) {
    if(method == "pvalue") {
      # examine_R_arguments(R)
      # seed <- examine_seed_arguments(seed = seed)
      # set.seed(seed = seed)
      stop("After version 1.2.0, 'pvalue' method is no longer supported.")
    }
    # data prepare for screening:
    if(dst == FALSE) {
      if(y_p != 1) {
        y <- as.vector(as.matrix(dist(y, diag = TRUE)))
        dst <- TRUE
      }
    } else {
      y <- as.vector(y)
    }
    # BCor-SIS:
    rcory_result <- apply_bcor_wrap(x = x, y = y, n = n, p = p, 
                                    dst = dst, weight = weight, 
                                    method = method, num.threads = num.threads)
    Xhavepickout <- get_screened_vars(ids, rcory_result[[2]], final_d)
    complete_info[[1]] <- rcory_result[[1]]
    # extra method for interaction:
    Xhavepickout2 <- c()
    if(method == "interaction") {
      rcory2_result <- apply_bcor_wrap(x = (x)^2, y = y, n = n, p = p, 
                                       dst = dst, weight = weight, 
                                       method = method, num.threads = num.threads)
      Xhavepickout2 <- get_screened_vars(ids, rcory2_result[[2]], final_d)
      complete_info[[2]] <- rcory_result[[2]]
    }
    Xhavepickout <- unique(c(Xhavepickout, Xhavepickout2))
  }
  if(method %in% c("gam", "lm")) {
    # data prepare for screening:
    y_copy <- preprocess_bcorsis_y(y, y_p)
    dst <- y_copy[[2]]
    y_copy <- y_copy[[1]]
    R <- 0
    # Initial screening:
    rcory_result <- apply_bcor_wrap(x = x, y = y_copy, n = n, p = p, 
                                    dst = dst, weight = weight, 
                                    method = method, num.threads = num.threads)
    # complete_info[[1]] <- rcory_result[[1]]
    # get d1 variables as initial variables set:
    Xhavepickout <- get_screened_vars(ids, rcory_result, d1)
    Xlastpickout <- Xhavepickout
    ids <- setdiff(ids, Xhavepickout)
    # Iterative:
    if(method == 'lm'){
      while(length(Xhavepickout) < final_d)
      {
        # lm fit for x
        Xnew <- residuals(lm(x[, ids] ~ x[, Xhavepickout]))
        # lm fit for y
        y <- residuals(lm(y ~ x[, Xlastpickout]))
        
        # BCor-screening
        y_copy <- preprocess_bcorsis_y(y, y_p)[[1]]
        rcory_result <- apply_bcor_wrap(x = Xnew, y = y_copy, n = n, p = p, 
                                        dst = dst, weight = weight, 
                                        method = method, num.threads = num.threads)
        # get d2 variables for each iteration:
        Xlastpickout <- get_screened_vars(ids, rcory_result, d2)
        Xhavepickout <- c(Xhavepickout, Xlastpickout)
        ids <- setdiff(ids, Xlastpickout)
      }
    }
    if(method == 'gam'){
      while(length(Xhavepickout) < final_d)
      {
        # gam fit for x
        lastpickout_formula <- paste0(' + s(',colnames(x)[Xlastpickout], collapse = paste0(", df = ", df, ")"))
        lastpickout_formula <- paste0(lastpickout_formula, paste0(", df = ", df, ")"), collapse = "")
        lastpickout_dat <- x[, Xlastpickout]
        Xnew <- sapply(ids, function(index){
          formula_one <- paste0(colnames(x)[index], "~", lastpickout_formula)
          formula_one <- as.formula(formula_one)
          dat <- as.data.frame(cbind(x[, index], lastpickout_dat))
          colnames(dat)[1] <- colnames(x)[index]
          # colnames(dat) <- paste0("x",c(x,Xhavepickout))
          gam::gam(formula_one, data = dat)[["residuals"]]
        })
        
        # gam fit for y
        dat <- data.frame("y" = y, lastpickout_dat)
        formula_Y <- as.formula(paste("y ~ ", lastpickout_formula))
        y <- gam::gam(formula = formula_Y,data = dat)$residuals
        
        # BCor-screening
        y_copy <- preprocess_bcorsis_y(y, y_p)[[1]]
        rcory_result <- apply_bcor_wrap(x = Xnew, y = y_copy, n = n, p = p, 
                                        dst = dst, weight = weight, 
                                        method = method, num.threads = num.threads)
        # get d2 variables for each iteration:
        Xlastpickout <- get_screened_vars(ids, rcory_result, d2)
        Xhavepickout <- c(Xhavepickout, Xlastpickout)
        ids <- setdiff(ids, Xlastpickout)
      }
    }
  }
  # return:
  list("ix" = Xhavepickout, "method" = method, 
       "weight" = weight, "complete.info" = complete_info)
}


#' @title Ball Correlation Sure Independence Screening For Survival data
#' @description Utilize extension of Ball Correlation in survival to select d variables related to survival status.
#' @inheritParams bcorsis
#' @param y a numeric matirx(first column should be event time, second column should be survival status) or Surv object
#' @param standized allows the user to standardize the covariate
#' @return the ids of selected variables
#' @noRd
#'
bcorsis.surv <- function(y, x, final_d, n, p, ids, standized = TRUE){
  # prepare for screening
  time <- y[, 1]
  delta <- y[, 2]
  ord.t <- sort(time)
  ix <- order(time)
  ord.delta <- delta[ix]
  x <- x[ix,]
  if(standized) {
    x <- apply(x, 2, scale)
  }
  
  # BCor Screening(survival)
  fitc <- survival::survfit(Surv(time, 1 - delta) ~ 1)
  Sc <- fitc[["surv"]]
  if(length(unique(ord.t)) != n) {
    rep_num <- as.data.frame(table(ord.t))[, "Freq"]
    Sc <- mapply(function(x, y) {
      rep(x, y)
    }, Sc, rep_num, SIMPLIFY = FALSE)
    Sc <- unlist(Sc)
  }
  
  t_rank <- rank(ord.t, ties.method = "max") - 1
  rcory_result <- apply(x, 2, function(x){
    bcor_surv(x = x, time_value = t_rank, delta = ord.delta, Sc = Sc, n = n)
  })
  
  Xhavepickout <- get_screened_vars(ids, rcory_result, final_d)
  list(Xhavepickout, rcory_result)
}

