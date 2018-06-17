\name{Ball-package}
\alias{Ball-package}
\alias{Ball}
\docType{package}
\title{Statistical Inference and Sure Independence Screening via Ball Statistics}
\description{Hypothesis tests and sure independence screening (SIS) procedure based on ball 
statistics, including ball divergence, ball covariance, and ball correlation, 
are developed to analyze complex data. The ball divergence and ball covariance based distribution-free tests are implemented to examine equality of multivariate distributions and independence between random vectors of arbitrary dimensions. Furthermore, a generic non-parametric SIS 
procedure based on ball correlation and all of its variants are implemented to tackle the 
challenge in the context of ultra high dimensional data.}
\details{
\tabular{ll}{
Package: \tab Ball\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2017-12-28\cr
License: \tab GPL-3\cr
}

The package contains three major functions:

\code{\link{bd.test}} - The ball divergence based permutation test for equality of multivariate (or univariate) distributions.

\code{\link{bcov.test}} - The ball covariance based permutation test for independence of multivariate (or univariate) random variable.

\code{\link{bcorsis}} - The ball correlation based non-parametric sure independence screening procedure for selecting the potential significant variables in the context of ultra-high dimensional data.

See \code{vignette("Ball")} for additional information. 

}

\author{
Xueqin Wang, Wenliang Pan, Heping Zhang, Hongtu Zhu, Yuan Tian, Weinan Xiao,
Chengfeng Liu and Jin Zhu

Maintainer: Jin Zhu <zhuj37@mail2.sysu.edu.cn>
}

\references{
Pan, Wenliang; Tian, Yuan; Wang, Xueqin; Zhang, Heping. Ball Divergence: Nonparametric two sample test. Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. https://projecteuclid.org/euclid.aos/1525313077
Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic Sure Independence Screening Procedure, Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
}

\examples{

\dontrun{

################# Quick Start for bd.test function #################
x <- rnorm(50)
y <- rnorm(50, mean = 1)
# plot(density(x))
# lines(density(y), col = 'red')
# ball divergence:
bd.test(x = x, y = y)

################# Quick Start for bcov.test function #################
error <- runif(50, min = -0.3, max = 0.3)
x <- runif(50, 0, 4*pi)
y <- cos(x) + error
# plot(x, y)
bcov.test(x, y)

################# Ball Divergence Test for Non-Hilbert Data #################
# load data:
data("macaques")
# number of femala and male Macaca fascicularis:
table(macaques[["group"]])
# calculate Riemannian shape distance matrix:
Dmat <- nhdist(macaques[["x"]], method = "riemann")
# hypothesis test with BD:
bd.test(x = Dmat, R = 99, size = c(9, 9), dst = TRUE)

################# Ball Covariance Test for Non-Hilbert Data #################
# load data:
data("ArcticLake")
# Distance matrix between y:
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
# Distance matrix between x:
Dx <- dist(ArcticLake[["depth"]])
# hypothesis test with BCov:
bcov.test(x = Dx, y = Dy, R = 99, dst = TRUE)

################  K-sample Test  #################
n <- 150
bd.test(rnorm(n), size = rep(50, 3))
# alternative input method:
x <- lapply(rep(50, 3), rnorm)
bd.test(x)

################  Weighted Ball Covariance Test  #################
data("ArcticLake")
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
Dx <- dist(ArcticLake[["depth"]])
# hypothesis test with weighted BCov:
bcov.test(x = Dx, y = Dy, R = 99, dst = TRUE, weight = TRUE)

################# Mutual Independence Test #################
x <- rnorm(30)
y <- (x > 0) * x + rnorm(30)
z <- (x <= 0) * x + rnorm(30)
data_list <- list(x, y, z)
bcov.test(data_list)

################# Mutual Independence Test for Meteorology data #################
data("meteorology")
bcov.test(meteorology)

############### Quick Start for bcorsis function ###############
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
error <- rnorm(n)
y <- 3*x[, 1] + 5*(x[, 3])^2 + error
res <- bcorsis(y = y, x = x)
head(res[[1]])

############### BCor-SIS: Censored Data Example ###############
result <- bcorsis(x = genlung[["covariate"]], y = genlung[["survival"]], 
                  method = "survival")$ix
top_gene <- colnames(genlung[["covariate"]])[result]
head(top_gene, n = 1)

############### BCor-SIS: Interaction Pursuing ###############
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
error <- rnorm(n)
y <- 3*x[, 1]*x[, 5]*x[, 10] + error
res <- bcorsis(y = y, x = x, method = "interaction")
head(res[[1]])

############### BCor-SIS: Iterative Method ###############
library(mvtnorm)
set.seed(1)
n <- 150
p <- 3000
sigma_mat <- matrix(0.5, nrow = p, ncol = p)
diag(sigma_mat) <- 1
x <- rmvnorm(n = n, sigma = sigma_mat)
error <- rnorm(n)
rm(sigma_mat); gc(reset = TRUE)
y <- 3*(x[, 1])^2 + 5*(x[, 2])^2 + 5*x[, 8] - 8*x[, 16] + error
res <- bcorsis(y = y, x = x, method = "gam", d = 15)
res[[1]]
}

}



