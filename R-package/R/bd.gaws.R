#' @title Fast K-sample Ball Divergence Test for GWAS Data 
#' @inheritParams bd.test
#' @param x a numeric vector, matrix, data.frame, dist object. 
#' @param snp a SNP dataset. Each column must be an integer vector.
#' @param refine a logical value. If \code{refine = TRUE}, a \eqn{p}-values refining process is applied to 
#' the SNPs which passes the pre-screening process. Default: \code{refine = TRUE} (At present, \code{refine = FALSE} is not available).
#' @param num.permutations the number of permutation replications. When \code{num.permutations = 0}, 
#' the function just returns the Ball Divergence statistic. Default: \code{num.permutations = 100 * ncol(snp)}
#' @param screening.method if \code{screening.method = "spectrum"}, the spectrum method is applied to 
#' screening the candidate SNPs, or otherwise, the permutation method is applied. Default: \code{screening.method = "permute"}.
#' @param screening.result A object return by \code{bd.gwas.test} that 
#' preserving the pre-screening result. 
#' It works only if the pre-screening is available.
#' Default: \code{screening.result = NULL}.
#' @param alpha the significance level. Default: \code{0.05 / ncol(snp)}.
#' @param verbose Show computation status and estimated runtimes. Default: \code{verbose = FALSE}.
#' 
#' @return bd.gwas.test returns a list containing the following components:
#' \item{\code{statistic}}{ball divergence statistics vector.}            
#' \item{\code{permuted.statistic}}{a data.frame containing permuted ball divergence statistic for pre-screening SNPs. 
#' If \code{refine = FALSE}, it takes value \code{NULL}.}            
#' \item{\code{eigenvalue}}{the eigenvalue of spectrum decomposition. If \code{refine = TRUE}, it takes value \code{NULL}.}
#' \item{\code{p.value}}{the p-values of ball divergence test.}
#' \item{\code{refined.snp}}{the SNPs have been refined.}
#' \item{\code{refined.p.value}}{the refined \eqn{p}-value of significant snp.}
#' \item{\code{refined.permuted.statistic}}{a data.frame containing permuted ball divergence statistics for refining \eqn{p}-values.}
#' 
#' @export
#' 
#' @references Hu, Y., Tan, H., Li, C., & Zhang, H. (2021). Identifying genetic risk variants associated with brain volumetric phenotypes via K‐sample Ball Divergence method. Genetic Epidemiology, 1–11. https://doi.org/10.1002/gepi.22423
#' 
#' @author Jin Zhu
#' 
#' @seealso 
#' \code{\link{bd}}, \code{\link{bd.test}}
#'
#' @examples
#' \donttest{
#' library(Ball)
#' set.seed(1234)
#' num <- 200
#' snp_num <- 500
#' p <- 5
#' x <- matrix(rnorm(num * p), nrow = num)
#' snp <- sapply(1:snp_num, function(i) {
#'   sample(0:2, size = num, replace = TRUE)
#' })
#' snp1 <- sapply(1:snp_num, function(i) {
#'   sample(1:2, size = num, replace = TRUE)
#' })
#' snp <- cbind(snp, snp1)
#' res <- Ball::bd.gwas.test(x = x, snp = snp)
#' mean(res[["p.value"]] < 0.05)
#' mean(res[["p.value"]] < 0.005)
#' 
#' ## only return the test statistics;
#' res <- Ball::bd.gwas.test(x = x, snp = snp, num.permutation = 0)
#' 
#' ## save pre-screening process results:
#' x <- matrix(rnorm(num * p), nrow = num)
#' snp <- sapply(1:snp_num, function(i) {
#'   sample(0:2, size = num, replace = TRUE, prob = c(1/2, 1/4, 1/4))
#' })
#' snp_screening <- Ball::bd.gwas.test(x = x, snp = snp, 
#'                                     alpha = 5*10^-4, refine = FALSE,
#'                                     num.permutations = 19999)
#' mean(res[["p.value"]] < 0.05)
#' mean(res[["p.value"]] < 0.005)
#' mean(res[["p.value"]] < 0.0005)
#' ## refine p-value according to the pre-screening process result:
#' res <- Ball::bd.gwas.test(x = x, snp = snp, alpha = 5*10^-4, 
#'                           num.permutations = 19999, 
#'                           screening.result = snp_screening[["screening.result"]])
#' }
bd.gwas.test <- function(x, snp, screening.method = c("permute", "spectrum"), 
                         refine = TRUE,
                         num.permutations, distance = FALSE, alpha, 
                         screening.result = NULL, 
                         verbose = TRUE, seed = 1, num.threads = 0, ...) 
{
  snp <- as.matrix(snp)
  num <- nrow(snp)
  snp_num <- ncol(snp)
  
  distance <- ifelse(class(x)[1] == "dist", TRUE, FALSE)
  if (distance) {
    if (class(x)[1] == "dist") {
      x <- as.vector(x)
    } else {
      x <- x[lower.tri(x)]
    }
  } else {
    x <- as.vector(stats::dist(x))
  }
  if ((0.5 * num * (num - 1)) != length(x)) {
    stop("sample size of x and snp are not match!")
  }
  x <- stats::na.fail(x)
  snp <- stats::na.fail(snp)
  
  snp_class_num <- apply(snp, 2, function(x) {
    # dplyr::n_distinct(x)
    length(unique(x))
  })
  if (any(snp_class_num == 1)) {
    stop("there are some snps only contain 1 group!")
  }
  # unique_class_num <- dplyr::n_distinct(snp_class_num)
  unique_class_num <- length(unique(snp_class_num))
  each_class_num <- as.vector(table(snp_class_num))
  
  if (length(screening.method) > 1) {
    screening.method <- "permute"
  } else {
    screening.method <- match.arg(screening.method)
  }
  
  if (screening.method == "permute") {
    if (missing(num.permutations)) {
      num.permutations <- 100 * snp_num
    }
    r <- num.permutations
  } else {
    r <- 0
  }
  
  if (r == 0) {
    verbose <- FALSE
  }
  
  eigenvalue <- NULL
  p_value <- NULL
  if (is.null(screening.result)) {
    set.seed(seed)
    statistic <- as.double(numeric(snp_num * 2))
    permuted_statistic <- as.double(numeric(r * 2 * unique_class_num))
    p_value <- as.double(numeric(snp_num * 2) + 1)
    x_index <- as.integer(numeric(num * num))
    ties <- integer(1)
    x <- as.double(x)
    snp_vec <- as.integer(snp)
    num <- as.integer(num)
    snp_num <- as.integer(snp_num)
    unique_class_num <- as.integer(unique_class_num)
    each_class_num <- as.integer(each_class_num)
    r <- as.integer(r)
    nth <- as.integer(num.threads)
    verbose_out <- as.integer(verbose)
    
    screen_res <- .C("bd_gwas_screening", statistic, permuted_statistic, p_value, x_index, ties,
                     x, snp_vec, num, snp_num, unique_class_num, each_class_num, 
                     r, nth, verbose_out)
    
    statistic <- screen_res[[1]][1:snp_num]
    permuted_statistic <- data.frame()
    x_index <- screen_res[[4]]
    ties <- screen_res[[5]]
    if (screening.method == "permute") {
      if (num.permutations > 0) {
        permuted_statistic <- data.frame(matrix(screen_res[[2]], nrow = r))
        permuted_statistic <- permuted_statistic[, seq(1, 2 * unique_class_num, by = 2), drop = FALSE]
        colnames(permuted_statistic) <- paste0("g", sort(unique(snp_class_num)))
        p_value <- screen_res[[3]][1:snp_num]
      }
    } else {
      # TODO: save the eigenvalues of spectrum method
      eigenvalue <- NULL
      p_value <- NULL
    }
    screening_result <- list(statistic, permuted_statistic, 
                             p_value, x_index, ties)
  } else {
    statistic <- screening.result[["statistic"]]
    permuted_statistic <- screening.result[["permuted_statistic"]]
    p_value <- screening.result[["p_value"]]
    x_index <- screening.result[["x_index"]]
    ties <- screening.result[["ties"]]
    screening_result <- screening.result
  }
  
  significant_snp <- NULL
  significant_snp_p_value <- NULL
  refine_permuted_statistic <- data.frame()
  refine_snp_index <- NULL
  refine_p_value_vector <- NULL
  refine_permuted_statistic_matirx <- NULL
  if (refine) {
    if (!is.null(p_value)) {
      if (missing(alpha)) {
        alpha <- 0.05 / snp_num;
      }
      refine_snp_index <- which(p_value < alpha)
      refine_snp_num <- length(refine_snp_index)
      if (refine_snp_num == 0) {
        if (verbose) {
          print("None of SNP pass the pre-screening process!")
        }
      } else {
        set.seed(seed)
        refine_permuted_statistic_matirx <- matrix(nrow = r, ncol = refine_snp_num)
        refine_p_value_vector <- numeric(refine_snp_num)
        x_index <- as.integer(x_index)
        ties <- as.integer(ties)
        x <- as.double(x)
        num <- as.integer(num)
        refine_snp_num <- as.integer(refine_snp_num)
        r <- as.integer(r)
        nth <- as.integer(num.threads)
        verbose_out <- as.integer(verbose)
        for (i in 1:length(refine_snp_index)) {
          refine_snp_statistic <- as.double(c(statistic[refine_snp_index[i]], statistic[refine_snp_index[i]]))
          refine_permuted_statistic <- as.double(numeric(r * 2))
          refine_p_value <- as.double(numeric(2))
          refine_snp_size_vec <- as.integer(table(snp[, refine_snp_index[i]]))
          refine_i_th <- as.integer(i)
          refine_k_num <- as.integer(snp_class_num[refine_snp_index[i]])
          
          refine_res <- .C("bd_gwas_refining_single", refine_snp_statistic, refine_permuted_statistic, refine_p_value, x_index, ties,
                           x, num, refine_snp_size_vec, refine_i_th, refine_k_num, 
                           refine_snp_num, r, nth, verbose_out)
          refine_permuted_statistic_matirx[, i] <- refine_res[[2]][seq(1, 2 * r, by = 2)]
          refine_p_value_vector[i] <- refine_res[[3]][1]
        }
        
        colnames(refine_permuted_statistic_matirx) <- paste0("SNP", refine_snp_index)
        p_value[refine_snp_index] <- refine_p_value_vector
        if (any(refine_p_value_vector < alpha)) {
          significant_snp <- which(p_value < alpha)
          significant_snp_p_value <- p_value[p_value < alpha]
        }
      }
    } 
  }
  
  return(list("statistic" = statistic, 
              "permuted.statistic" = permuted_statistic, 
              "eigenvalue" = eigenvalue,
              "p.value" = p_value, 
              "refined.snp" = refine_snp_index, 
              "refined.p.value" = refine_p_value_vector,
              "refined.permuted.statistic" = refine_permuted_statistic_matirx, 
              "screening.result" = screening_result))
}