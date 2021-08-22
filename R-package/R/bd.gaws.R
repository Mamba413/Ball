#' @title Fast Two Ball Divergence test for gwas data 
#' @inheritParams bd.test
#' @param x a numeric vector, matrix, data.frame, dist object 
#' @param snp SNP data
#'
#' @return bd.gwas.test returns a list containing the following components:
#' \item{\code{statistic}}{ball divergence statistic.}            
#' \item{\code{permuted_statistic}}{a data.frame containing permuted ball divergence statistic.}            
#' \item{\code{p.value}}{the p-value for the test.}
#' 
#' @export
#'
#' @examples
#' set.seed(1)
#' num <- 60
#' snp_num <- 100
#' p <- 10
#' x1 <- matrix(rnorm(num * p / 2), nrow = num/2)
#' x2 <- matrix(rnorm(num * p / 2, 3), nrow = num/2)
#' x <- rbind(x1, x2)
#' snp <- sapply(1:snp_num, function(i) { 
#'   set.seed(i)
#'   sample(c(0, 1), size = num, replace = TRUE)
#' })
#' snp[, 1] <- rep(c(0, 1), each = num/2)
#' snp[, 2] <- c(rep(0, num/2 - 1), rep(1, num/2 + 1))
#' snp[, 3] <- c(rep(0, num/2 - 2), rep(1, num/2 + 2))
#' snp[, 4] <- c(rep(0, num/2 - 4), rep(1, num/2 + 4))
#' res <- bd.gwas.test(x = x, snp = snp, num.threads = 2)
#' head(res[["p.value"]], n = 8)
bd.gwas.test <- function(x, snp, dst = FALSE, seed = 4, num.threads = 2, ...) {
  snp <- as.matrix(snp)
  num <- nrow(snp)
  snp_num <- ncol(snp)
  r <- floor(snp_num / 0.05) * 10
  group1_num <- apply(snp, 2, sum)
  if (anyNA(group1_num)) {
    stop("")
  }
  number_threshold <- floor(nrow(snp)/2)
  convert_index <- group1_num > number_threshold
  sapply(1:length(convert_index), function(i) {
    if (convert_index[i]) {
      snp[, i] <<- 1 - snp[, i]
    }
  })
  group1_num <- ifelse(convert_index, nrow(snp) - group1_num, group1_num)
  unique_group_num <- as.numeric(names(sort(table(group1_num), decreasing = TRUE)))
  
  #
  if (dst) {
    dx <- as.matrix(x)
  } else {
    dx <- as.matrix(dist(x))
  }
  
  # 
  set.seed(seed = seed)
  statistic <- as.double(numeric(snp_num))
  permuted_statistic <- as.double(numeric(r * length(unique_group_num)))
  p.value <- as.double(numeric(snp_num))
  snp_vec <- as.integer(snp)
  unique_group_num_vec <- as.integer(unique_group_num)
  unique_group_num_length <- as.integer(length(unique_group_num))
  n1_num_vec <- as.integer(group1_num)
  snp_number <- as.integer(snp_num)
  xy <- as.double(dx)
  r <- as.integer(r)
  n <- as.integer(num)
  dst <- as.integer(1)
  nthread <- as.integer(num.threads)
  
  res <- .C("bd_gwas_test", statistic, permuted_statistic, p.value, 
            snp_vec, unique_group_num_vec, unique_group_num_length, 
            n1_num_vec, snp_number, xy, r, n, dst, nthread)
  permuted_statistic <- data.frame(matrix(res[[2]], nrow = r))
  colnames(permuted_statistic) <- unique_group_num_vec
  return(list("statistic" = res[[1]], 
              "permuted_statistic" = permuted_statistic, 
              "p.value" = res[[3]]))
}

#' @title Fast Two Ball Divergence test for gwas data 
#' @noRd
#' @examples
#' set.seed(1)
#' num <- 100
#' snp_num <- 100
#' p <- 10
#' x <- matrix(rnorm(num * p), nrow = num)
#' snp <- sapply(1:snp_num, function(i) { 
#'   set.seed(i)
#'   sample(c(0, 1), size = num, replace = TRUE)
#' })
#' Sys.time()
#' res <- bd.gwas.test1(x = x, snp = snp)
#' Sys.time()
#' 
bd.gwas.test1 <- function(x, snp, dst = FALSE, seed = 4, num.threads = 1, ...) {
  snp <- as.matrix(snp)
  num <- nrow(snp)
  snp_num <- ncol(snp)
  r <- floor(snp_num / 0.05) * 10
  group1_num <- apply(snp, 2, sum)
  if (anyNA(group1_num)) {
    stop("")
  }
  number_threshold <- floor(nrow(snp)/2)
  group1_num <- ifelse(group1_num > number_threshold, nrow(snp) - group1_num, group1_num)
  unique_group_num <- as.numeric(names(sort(table(group1_num), decreasing = TRUE)))
  
  #
  if (dst) {
    dx <- as.matrix(x)
  } else {
    dx <- as.matrix(dist(x))
  }
  
  # 
  statistic <- numeric(ncol(snp))
  p.value <- numeric(ncol(snp))
  permuted_statistic <- list()
  for (i_group_num in unique_group_num) {
    snp_index <- which(group1_num == i_group_num)
    statistic[snp_index] <- apply(snp[, snp_index, drop = FALSE], 2, function(index) {
      order_index <- order(index)
      bd(x = dx[order_index, order_index], dst = TRUE, 
         size = as.vector(table(index)), num.threads = 1)
    })
    permuted_statistic[[as.character(i_group_num)]] <- sapply(1:r, function(i) {
      set.seed(i)
      random_index <- sample(1:num, size = num, replace = FALSE)
      bd(x = dx[random_index, random_index], 
         dst = TRUE, size = c(i_group_num, num - i_group_num), num.threads = 1)
    })
    p.value[snp_index] <- sapply(statistic[snp_index], calculatePvalue, 
                                 NullDistribution = permuted_statistic[[as.character(i_group_num)]])
  }
  permuted_statistic <- do.call("cbind.data.frame", permuted_statistic)
  return(list("statistic" = statistic, 
              "permuted_statistic" = permuted_statistic, 
              "p.value" = p.value))
}


