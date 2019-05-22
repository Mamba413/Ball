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
#' @importFrom dplyr n_distinct
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
#'   sample(0:2, size = num, replace = TRUE)
#' })
#' snp1 <- sapply(1:snp_num, function(i) {
#'   sample(1:2, size = num, replace = TRUE)
#' })
#' snp <- cbind(snp, snp1)
#' res <- Ball::bd.gwas.test(x = x, snp = snp, num.permutations = 19999)
#' mean(res[["p.value"]] < 0.05)
#' mean(res[["p.value"]] < 0.005)
#' mean(res[["p.value"]] < 0.0005)
#' }
bd.gwas.test <- function(x, snp, screening.method = c("permute", "spectrum"), 
                         num.permutations, distance = FALSE, 
                         seed = 1, num.threads = 0, ...) 
{
  snp <- as.matrix(snp)
  num <- nrow(snp)
  snp_num <- ncol(snp)
  
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
  x <- na.fail(x)
  snp <- na.fail(snp)
  
  snp_class_num <- apply(snp, 2, function(x) {
    dplyr::n_distinct(x)
  })
  if (any(snp_class_num == 1)) {
    stop("there are some snps only contain 1 group!")
  }
  unique_class_num <- dplyr::n_distinct(snp_class_num)
  each_class_num <- as.vector(table(snp_class_num))
  
  if (length(screening.method) > 1) {
    screening.method <- "permute"
  } else {
    screening.method <- match.arg(screening.method)
  }
  
  if (screening.method == "permute") {
    if (missing(num.permutations)) {
      r <- floor(snp_num / 0.05) * 5
    } else {
      r <- num.permutations
    }
  } else {
    r <- 0
  }
  
  set.seed(seed)
  statistic <- as.double(numeric(snp_num * 2))
  permuted_statistic <- as.double(numeric(r * 2 * unique_class_num))
  p_value <- as.double(numeric(snp_num * 2))
  x <- as.double(x)
  snp_vec <- as.integer(snp)
  num <- as.integer(num)
  snp_num <- as.integer(snp_num)
  unique_class_num <- as.integer(unique_class_num)
  each_class_num <- as.integer(each_class_num)
  r <- as.integer(r)
  nth <- as.integer(num.threads)
  
  screen_res <- .C("bd_gwas_screening", statistic, permuted_statistic, p_value, 
                   x, snp_vec, num, snp_num, unique_class_num, each_class_num, r, nth)
  
  if (screening.method == "permute") {
    permuted_statistic <- data.frame(matrix(screen_res[[2]], nrow = r))
    permuted_statistic <- permuted_statistic[, seq(1, 2 * unique_class_num, by = 2), drop = FALSE]
    colnames(permuted_statistic) <- paste0("g", sort(unique(snp_class_num)))
  }

  return(list("statistic" = screen_res[[1]][1:snp_num], 
              "permuted.statistic" = permuted_statistic, 
              "p.value" = screen_res[[3]][1:snp_num]))
}