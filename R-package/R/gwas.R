#' @title Fast K-sample Ball Divergence Test for GWAS Data 
#' @inheritParams bd.test
#' @param x a numeric vector, matrix, data.frame, dist object 
#' @param snp SNP data
#' 
#' @return bd.gwas.test returns a list containing the following components:
#' \item{\code{statistic}}{ball divergence statistics vector.}            
#' \item{\code{permuted.statistic}}{a data.frame containing permuted ball divergence statistic.}            
#' \item{\code{p.value}}{the p-values for the ball divergence test.}
#' 
#' @export
#' 
#' @seealso 
#' \code{\link{bd}}, \code{\link{bd.test}}
#'
#' @examples
#' \dontrun{
#' 
#' set.seed(1)
#' num <- 200
#' snp_num <- 10000
#' p <- 5
#' x <- matrix(rnorm(num * p), nrow = num)
#' snp <- sapply(1:snp_num, function(i) {
#'   sample(0:1, size = num, replace = TRUE)
#' })
#' res <- bd.gwas.test(x = x, snp = snp, num.permutations = 19999)
#' mean(res[["p.value"]] < 0.05)
#' mean(res[["p.value"]] < 0.005)
#' mean(res[["p.value"]] < 0.0005)
#' }
bd.gwas.test <- function(x, snp, method = c("permute", "spectrum"), 
                         num.permutations, distance = FALSE, 
                         seed = 1, num.threads = 0, ...) 
{
  snp <- as.matrix(snp)
  num <- nrow(snp)
  snp_num <- ncol(snp)
  
  if (length(method) > 1) {
    method <- "permute"
  } else {
    method <- match.arg(method)
  }
  
  if (method == "permute") {
    if (missing(num.permutations)) {
      r <- floor(snp_num / 0.05) * 5
    } else {
      r <- num.permutations
    }
  } else {
    r <- 0
  }

  distance <- ifelse(class(x) == "dist", TRUE, FALSE)
  if (distance) {
    if (class(x) == "dist") {
      x <- as.vector(x)
    } else {
      x <- x[lower.tri(x)]
    }
  } else {
    x <- as.vector(dist(x))
  }
  if ((0.5 * num * (num - 1)) != length(x)) {
    stop("sample size of x and snp are not match!")
  }
  
  set.seed(seed)
  statistic <- as.double(numeric(snp_num * 2))
  permuted_statistic <- as.double(numeric(r * 2))
  p_value <- as.double(numeric(snp_num * 2))
  x <- as.double(x)
  snp_vec <- as.integer(snp)
  num <- as.integer(num)
  snp_num <- as.integer(snp_num)
  r <- as.integer(r)
  nth <- as.integer(num.threads)
  
  res <- .C("bd_gwas_test", statistic, permuted_statistic, p_value, 
            x, snp_vec, num, snp_num, r, nth)
  
  if (method == "permute") {
    permuted_statistic <- data.frame(matrix(res[[2]], nrow = r, ncol = 2))
    permuted_statistic <- permuted_statistic[, 1, drop = FALSE]    
  }

  return(list("statistic" = res[[1]][1:snp_num], 
              "permuted.statistic" = permuted_statistic, 
              "p.value" = res[[3]][1:snp_num]))
}