library(Ball)
library(NetworkDistance)
library(igraph)

######## Independence Test for network data ########
set.seed(1)
generate_batch_netword <- function(px)
{
  nd.wsd(lapply(px, function(y) {
    as.matrix(as_adjacency_matrix(sample_gnp(n = 10, p = y)))
  }), out.dist = FALSE)[["D"]]
}
size <- 100
x <- matrix(rt(n = size * 6, df = 3), ncol = 6)
y <- sin(x[, 1] * x[, 2] * x[, 3]) + cos(x[, 4] * x[, 5] * x[, 6])
y <- (atan(y) + pi/2)/pi
dy <- generate_batch_netword(y)
dx <- dist(x[, c(1, 2, 5, 6)])

bcov_res <- bcov.test(dx, dy, dst = TRUE)
bcov_res[["complete.info"]][["statistic"]]
bcov_res[["complete.info"]][["p.value"]]
