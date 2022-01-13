#' @inheritParams bcov.test
#' @rdname bcov
#' @return 
#' \item{\code{bcor }}{ Ball Correlation statistic.}
#' @export
#' @examples
#' ############# Ball Correlation #############
#' num <- 50
#' x <- 1:num
#' y <- 1:num
#' bcor(x, y)
#' bcor(x, y, weight = "prob")
#' bcor(x, y, weight = "chisq")
bcor <- function(x, y, distance = FALSE, weight = FALSE) {
  weight <- examine_weight_arguments(weight)
  x <- as.matrix(x)
  y <- as.matrix(y)
  x_y_info <- examine_x_y_bcor(x, y)
  p <- x_y_info[2]
  #
  if(distance == FALSE) {
    if(p != 1) {
      x <- as.double(as.vector(dist(x)))
      y <- as.double(as.vector(dist(y)))
      dst_y <- as.integer(1)
      dst_x <- as.integer(1)
    } else {
      x <- as.double(x)
      y <- as.double(y)
      dst_y <- as.integer(0)
      dst_x <- as.integer(0)
    }
  } else {
    x <- x[lower.tri(x)]
    y <- y[lower.tri(y)]
    x <- as.double(x)
    y <- as.double(y)
    dst_y <- as.integer(1)
    dst_x <- as.integer(1)
  }
  bcor_stat <- as.double(numeric(3))
  x_number <- as.integer(1)
  f_number <- as.integer(1)
  num <- as.integer(x_y_info[1])
  p <- as.integer(1)
  k <- as.integer(1)
  nth <- as.integer(1)
  missing_flag <- as.integer(rep(1, length(x)))
  res <- .C("bcor_test", bcor_stat, y, x, x_number, f_number, num, p, k, dst_y, dst_x, nth, missing_flag)
  bcor_stat <- res[[1]]
  bcor_stat <- select_ball_stat(bcor_stat, weight, type = "bcor")
  return(bcor_stat)
}


#' @title Ball Correlation based Sure Independence Screening (BCor-SIS)
#' @author Wenliang Pan, Weinan Xiao, Xueqin Wang, Hongtu Zhu, Jin Zhu
#' @description Generic non-parametric sure independence screening (SIS) procedure based on Ball Correlation.
#' Ball correlation is a generic measure of dependence in Banach spaces.
#' @inheritParams bcov.test
#' @param x a numeric matrix or data.frame included \eqn{n} rows and \eqn{p} columns. 
#' Each row is an observation vector and each column corresponding to a explanatory variable, generally \eqn{p >> n}.
#' @param d the hard cutoff rule suggests selecting \eqn{d} variables. Setting \code{d = "large"} or 
#' \code{d = "small"} means \code{n - 1} or \code{floor(n/log(n))} 
#' variables are selected. If \code{d} is a integer, \code{d} variables are selected. Default: \code{d = "small"}.
#' @param method specific method for the BCor-SIS procedure. It must be one of \code{"standard"},
#' \code{"lm"}, \code{"gam"}, \code{"interaction"}, or \code{"survival"}.
#' Setting \code{method = "standard"} means performing standard SIS procedure 
#' while the options \code{"lm"} and \code{"gam"} mean carrying out iterative SIS procedure with ordinary 
#' linear regression and generalized additive models, respectively.
#' The options \code{"interaction"} and \code{"survival"} are designed for detecting variables 
#' with potential linear interaction and associated with left censored responses, respectively. 
#' Any unambiguous substring can be given. Default: \code{method = "standard"}.
#' @param distance if \code{distance = TRUE}, \code{y} will be considered as a distance matrix. 
#' Arguments only available when \code{method = "standard"} and \code{method = "interaction"}. Default: \code{distance = FALSE}.
#' @param category a logical value or integer vector indicating columns to be selected as categorical variables.
#' If \code{category} is an integer vector, the positive/negative integers select/discard the corresponding columns;
#' If \code{category} is a logical value, \code{category = TRUE} select all columns, \code{category = FALSE} select none column.
#' Default: \code{category = FALSE}.
#' @param parms parameters list only available when \code{method = "lm"} or \code{"gam"}. 
#' It contains three parameters: \code{d1}, \code{d2}, and \code{df}. \code{d1} is the
#' number of initially selected variables, \code{d2} is the number of variables added in each iteration.
#' \code{df} is a degree freedom of basis in generalized additive models playing a role only when \code{method = "gam"}. 
#' Default: \code{parms = list(d1 = 5, d2 = 5, df = 3)}.
#' 
#' @return 
#' \item{\code{ix }}{ the indices vector corresponding to variables selected by BCor-SIS.} 
#' \item{\code{method }}{ the method used.} 
#' \item{\code{weight }}{ the weight used.} 
#' \item{\code{complete.info }}{ a \code{list} mainly containing a \eqn{p \times 3} matrix, 
#' where each row is a variable and each column is a weight Ball Correlation statistic. 
#' If \code{method = "gam"} or \code{method = "lm"}, \code{complete.info} is an empty list.} 
#' 
#' @details 
#' \code{bcorsis} performs a model-free generic sure independence screening procedure, 
#' BCor-SIS, to pick out variables from \code{x} which are potentially associated with \code{y}. 
#' BCor-SIS relies on Ball correlation, a universal dependence measure in Banach spaces.
#' Ball correlation (BCor) ranges from 0 to 1. A larger BCor implies they are likely to be associated while 
#' Bcor is equal to 0 implies they are unassociated. (See \code{\link{bcor}} for details.)
#' Consequently, BCor-SIS pick out variables with larger Bcor values with \code{y}.
#' 
#' Theory and numerical result indicate that BCor-SIS has following advantages:
#' \itemize{
#' \item BCor-SIS can retain the efficient variables even when the dimensionality (i.e., \code{ncol(x)}) is 
#' an exponential order of the sample size (i.e., \code{exp(nrow(x))});
#' \item It is distribution-free and model-free;
#' \item It is very robust;
#' \item It is works well for complex data, such as shape and survival data;
#' }
#' 
#' If \code{x} is a matrix, the sample sizes of \code{x} and \code{y} must agree.
#' If \code{x} is a \code{\link{list}} object, each element of this \code{list} must with the same sample size.
#' \code{x} and \code{y} must not contain missing or infinite values. 
#' 
#' When \code{method = "survival"}, the matrix or data.frame pass to \code{y} must have exactly two columns, where the first column is 
#' event (failure) time while the second column is a dichotomous censored status.
#' 
#' @note 
#' \code{bcorsis} simultaneously computing Ball Correlation statistics with 
#' \code{"constant"}, \code{"probability"}, and \code{"chisquare"} weights.
#' Users can get other Ball Correlation statistics with different weight in the \code{complete.info} element of output. 
#' We give a quick example below to illustrate. 
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
#' eps <- rnorm(n)
#' y <- 3 * x[, 1] + 5 * (x[, 3])^2 + eps
#' res <- bcorsis(y = y, x = x)
#' head(res[["ix"]])
#' head(res[["complete.info"]][["statistic"]])
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
#' eps <- rnorm(n)
#' y <- 3 * x[, 1] * x[, 5] * x[, 10] + eps
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
#' eps <- rnorm(n)
#' rm(sigma_mat); gc(reset = TRUE)
#' y <- 3 * (x[, 1])^2 + 5 * (x[, 2])^2 + 5 * x[, 8] - 8 * x[, 16] + eps
#' res <- bcorsis(y = y, x = x, method = "lm", d = 15)
#' res <- bcorsis(y = y, x = x, method = "gam", d = 15)
#' res[["ix"]]
#' 
#' ############### Weighted BCor-SIS: Probability weight ###############
#' set.seed(1)
#' n <- 150
#' p <- 3000
#' x <- matrix(rnorm(n * p), nrow = n)
#' eps <- rnorm(n)
#' y <- 3 * x[, 1] + 5 * (x[, 3])^2 + eps
#' res <- bcorsis(y = y, x = x, weight = "prob")
#' head(res[["ix"]])
#' # Alternative, chisq weight:
#' res <- bcorsis(y = y, x = x, weight = "chisq")
#' head(res[["ix"]])
#' 
#' ############### BCor-SIS: GWAS data ###############
#' set.seed(1)
#' n <- 150
#' p <- 3000
#' x <- sapply(1:p, function(i) {
#'   sample(0:2, size = n, replace = TRUE)
#' })
#' eps <- rnorm(n)
#' y <- 6 * x[, 1] - 7 * x[, 2] + 5 * x[, 3] + eps
#' res <- bcorsis(x = x, y = y, category = TRUE)
#' head(res[["ix"]])
#' head(res[["complete.info"]][["statistic"]])
#' 
#' x <- cbind(matrix(rnorm(n * 2), ncol = 2), x)
#' # remove the first two columns:
#' res <- bcorsis(x = x, y = y, category = c(-1, -2))
#' head(res[["ix"]])
#' 
#' x <- cbind(x[, 3:5], matrix(rnorm(n * p), ncol = p))
#' res <- bcorsis(x = x, y = y, category = 1:3)
#' head(res[["ix"]], n = 10)
#' }
bcorsis <- function(x, y, d = "small", weight = c("constant", "probability", "chisquare"), 
                    method = "standard", distance = FALSE, category = FALSE, 
                    parms = list(d1 = 5, d2 = 5, df = 3),
                    num.threads = 0)
{
  seed <- 1
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- examine_x_y_bcor(x, y)[1]
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
  
  # check category
  category_index <- examine_category(category, p)
  if (length(category_index) != 0 && (method %in% c("lm", "gam", "survival", "interaction"))) {
    stop("Handling category variables is only available when method = \"standard\".")
  }
  
  # get arguments:
  d1 <- parms$d1
  d2 <- parms$d2
  df <- parms$df
  
  # examine method arguments
  method <- examine_method_arguments(method)
  # examine dst and method arguments
  examine_dst_method(dst = distance, method = method)
  
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
    if(distance == FALSE) {
      if(y_p != 1) {
        y <- as.vector(dist(y))
        distance <- TRUE
      }
    } else {
      y <- y[lower.tri(y)]
    }
    # BCor-SIS:
    rcory_result <- apply_bcor_wrap(x = x, y = y, n = n, p = y_p, 
                                    distance = distance, weight = weight, 
                                    method = method, num.threads = num.threads, 
                                    category = category_index)
    Xhavepickout <- get_screened_vars(ids, rcory_result[[2]], final_d)
    complete_info[[1]] <- rcory_result[[1]]
    # extra method for interaction:
    Xhavepickout2 <- c()
    if(method == "interaction") {
      rcory2_result <- apply_bcor_wrap(x = (x)^2, y = y, n = n, p = y_p, 
                                       distance = distance, weight = weight, 
                                       method = method, num.threads = num.threads, 
                                       category = c())
      Xhavepickout2 <- get_screened_vars(ids, rcory2_result[[2]], final_d)
      complete_info[[2]] <- rcory_result[[2]]
    }
    Xhavepickout <- unique(c(Xhavepickout, Xhavepickout2))
  }
  if(method %in% c("gam", "lm")) {
    # data prepare for screening:
    y_copy <- preprocess_bcorsis_y(y, y_p)
    y_copy <- y_copy[[1]]
    distance <- y_copy[[2]]
    R <- 0
    # Initial screening:
    rcory_result <- apply_bcor_wrap(x = x, y = y_copy, n = n, p = y_p, 
                                    distance = distance, weight = weight, 
                                    method = method, num.threads = num.threads, 
                                    category = c())
    # complete_info[[1]] <- rcory_result[[1]]
    # get d1 variables as initial variables set:
    Xhavepickout <- get_screened_vars(ids, rcory_result, d1)
    Xlastpickout <- Xhavepickout
    ids <- setdiff(ids, Xhavepickout)
    # Iterative:
    if (method == 'lm') {
      while(length(Xhavepickout) < final_d)
      {
        # lm fit for x
        Xnew <- stats::residuals(stats::lm(x[, ids] ~ x[, Xhavepickout]))
        # lm fit for y
        y <- stats::residuals(stats::lm(y ~ x[, Xlastpickout]))
        
        # BCor-screening
        y_copy <- preprocess_bcorsis_y(y, y_p)[[1]]
        rcory_result <- apply_bcor_wrap(x = Xnew, y = y_copy, n = n, p = y_p, 
                                        distance = distance, weight = weight, 
                                        method = method, num.threads = num.threads, 
                                        category = c())
        # get d2 variables for each iteration:
        Xlastpickout <- get_screened_vars(ids, rcory_result, d2)
        Xhavepickout <- c(Xhavepickout, Xlastpickout)
        ids <- setdiff(ids, Xlastpickout)
      }
    }
    if (method == 'gam') {
      while(length(Xhavepickout) < final_d)
      {
        # gam fit for x
        lastpickout_formula <- paste0(' + gam::s(',colnames(x)[Xlastpickout], collapse = paste0(", df = ", df, ")"))
        lastpickout_formula <- paste0(lastpickout_formula, paste0(", df = ", df, ")"), collapse = "")
        lastpickout_dat <- x[, Xlastpickout]
        Xnew <- sapply(ids, function(index){
          formula_one <- paste0(colnames(x)[index], "~", lastpickout_formula)
          formula_one <- stats::as.formula(formula_one)
          dat <- as.data.frame(cbind(x[, index], lastpickout_dat))
          colnames(dat)[1] <- colnames(x)[index]
          # colnames(dat) <- paste0("x",c(x,Xhavepickout))
          suppressWarnings(residuals_value <- gam::gam(formula_one, data = dat)[["residuals"]])
          residuals_value
        })
        
        # gam fit for y
        dat <- data.frame("y" = y, lastpickout_dat)
        formula_Y <- as.formula(paste("y ~ ", lastpickout_formula))
        suppressWarnings(y <- gam::gam(formula = formula_Y, data = dat)$residuals)
        
        # BCor-screening
        y_copy <- preprocess_bcorsis_y(y, y_p)[[1]]
        rcory_result <- apply_bcor_wrap(x = Xnew, y = y_copy, n = n, p = y_p, 
                                        distance = distance, weight = weight, 
                                        method = method, num.threads = num.threads, 
                                        category = c())
        # get d2 variables for each iteration:
        Xlastpickout <- get_screened_vars(ids, rcory_result, d2)
        Xhavepickout <- c(Xhavepickout, Xlastpickout)
        ids <- setdiff(ids, Xlastpickout)
      }
    }
  }
  
  # return:
  complete_info[[2]] <- n
  complete_info[[3]] <- p
  names(complete_info) <- c("statistic", "n", "p")
  list("ix" = Xhavepickout, "method" = method, 
       "weight" = weight, "complete.info" = complete_info)
}


#' @title Ball Correlation Sure Independence Screening For Survival data
#' @description Utilize extension of Ball Correlation in survival to select d variables related to survival status.
#' @inheritParams bcorsis
#' @param y a numeric matrix (first column should be event time, second column should be survival status) or Surv object
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
  x <- x[ix, ]
  if(standized) {
    x <- apply(x, 2, scale)
  }
  
  # BCor Screening(survival)
  fitc <- survival::survfit(survival::Surv(time, 1 - delta) ~ 1)
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

