#' @title Fast K-sample Ball Divergence Test for GWAS Data 
#' @inheritParams bd.test
#' @param x a numeric vector, matrix, data.frame, dist object 
#' @param snp SNP data
#' @param screening.method if \code{screening.method = "spectrum"}, the spectrum method is applied to 
#' screening the candidate SNPs, or otherwise, the permutation method is applied.
#' @param alpha the significance level.
#' @param verbose Show computation status and estimated runtime. Default: \code{verbose = FALSE}.
#' 
#' @return bd.gwas.test returns a list containing the following components:
#' \item{\code{statistic}}{ball divergence statistics vector.}            
#' \item{\code{permuted.statistic}}{a data.frame containing permuted ball divergence statistic.}            
#' \item{\code{p.value}}{the p-values for the ball divergence test.}
#' 
#' @export
#' 
#' @importFrom dplyr n_distinct
#' @importFrom stats na.fail
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
#' 
#' num <- 100
#' snp_num <- 500
#' x <- rbind(matrix(rnorm(num * p), nrow = num), 
#'            matrix(rnorm(num * p, mean = 3), nrow = num), 
#'            matrix(rnorm(num * p, mean = -3), nrow = num))
#' snp <- sapply(1:snp_num, function(i) {
#'   sample(0:2, size = 3 * num, replace = TRUE)
#' })
#' snp[, 300] <- rep(0:2, each = num)
#' res <- Ball::bd.gwas.test(x = x, snp = snp)
#' 
#' set.seed(1234)
#' num <- 200
#' p <- 5
#' snp_num <- 10000
#' x <- matrix(rnorm(num * p), nrow = num)
#' snp <- sapply(1:snp_num, function(i) {
#'   sample(0:2, size = num, replace = TRUE, prob = c(1/2, 1/4, 1/4))
#' })
#' res <- Ball::bd.gwas.test(x = x, snp = snp, alpha = 5*10^-4, num.permutations = 19999, save.screening = "test.rda")
#' mean(res[["p.value"]] < 0.05)
#' mean(res[["p.value"]] < 0.005)
#' mean(res[["p.value"]] < 0.0005)
#' res <- Ball::bd.gwas.test(x = x, snp = snp, alpha = 5*10^-4, num.permutations = 19999, screening.file = "test.rda")
#' 
#' }
bd.gwas.test <- function(x, snp, screening.method = c("permute", "spectrum"), refine = TRUE,
                         num.permutations, distance = FALSE, alpha, 
                         save.screening = NULL, screening.file = NULL, 
                         verbose = TRUE, seed = 1, num.threads = 0, ...) 
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
  eigenvalue <- NULL
  p_value <- NULL
  if (is.null(screening.file)) {
    statistic <- as.double(numeric(snp_num * 2))
    permuted_statistic <- as.double(numeric(r * 2 * unique_class_num))
    p_value <- as.double(numeric(snp_num * 2) + 1)
    x <- as.double(x)
    snp_vec <- as.integer(snp)
    num <- as.integer(num)
    snp_num <- as.integer(snp_num)
    unique_class_num <- as.integer(unique_class_num)
    each_class_num <- as.integer(each_class_num)
    r <- as.integer(r)
    nth <- as.integer(num.threads)
    verbose_out <- as.integer(verbose)
    
    screen_res <- .C("bd_gwas_screening", statistic, permuted_statistic, p_value, 
                     x, snp_vec, num, snp_num, unique_class_num, each_class_num, 
                     r, nth, verbose_out)
    
    statistic <- screen_res[[1]][1:snp_num]
    permuted_statistic <- data.frame()
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
    if (!is.null(save.screening)) {
      save(statistic, permuted_statistic, p_value, file = save.screening)
    }
  } else {
    if (file.exists(screening.file)) {
      load(screening.file)
    } else {
      stop("the screening file not exists!")
    }
  }
  
  significant_snp <- NULL
  refine_permuted_statistic <- data.frame()
  refine_snp_index <- NULL
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
        refine_snp_statistic <- as.double(c(statistic[refine_snp_index], statistic[refine_snp_index]))
        refine_permuted_statistic <- as.double(numeric(r * 2 * refine_snp_num))
        refine_p_value <- as.double(numeric(2 * refine_snp_num))
        x <- as.double(x)
        num <- as.integer(num)
        r <- as.integer(r)
        nth <- as.integer(num.threads)
        verbose_out <- as.integer(verbose)
        refine_snp_num <- as.integer(refine_snp_num)
        refine_snp_size_vec <- as.integer(apply(snp[, refine_snp_index, drop = FALSE], 2, table))
        refine_snp_class_num <- as.integer(snp_class_num[refine_snp_index])
        
        refine_res <- .C("bd_gwas_refining", refine_snp_statistic, refine_permuted_statistic, refine_p_value,
                         x, num, refine_snp_num, refine_snp_size_vec, refine_snp_class_num,
                         r, nth, verbose_out)
        
        refine_permuted_statistic <- data.frame(matrix(refine_res[[2]], nrow = r))
        refine_permuted_statistic <- refine_permuted_statistic[, seq(1, 2 * refine_snp_num, by = 2), drop = FALSE]
        colnames(refine_permuted_statistic) <- paste0("SNP", sort(unique(refine_snp_index)))
        refine_p_value <- refine_res[[3]][1:refine_snp_num]
        p_value[refine_snp_index] <- refine_p_value
        if (any(refine_p_value < alpha)) {
          significant_snp <- which(p_value < alpha)
        }
      }
    } 
  }
  
  return(list("statistic" = statistic, 
              "permuted.statistic" = permuted_statistic, 
              "eigenvalue" = eigenvalue,
              "p.value" = p_value, 
              "significant.snp" = significant_snp, 
              "refined.snp" = refine_snp_index, 
              "refined.permuted.statistic" = refine_permuted_statistic))
}