## ---- echo=FALSE, message=FALSE, warning=FALSE---------------------------
knitr::opts_chunk$set(comment = "#", warning = FALSE, eval = TRUE, message = FALSE)
set.seed(1)
library(Ball)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("Ball")

## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github("Mamba413/Ball", build_vignettes = TRUE)

## ---- echo=FALSE---------------------------------------------------------
library(Ball)

## ------------------------------------------------------------------------
x <- rnorm(50)
y <- rnorm(50, mean = 1)
# plot(density(x), xlim = c(-5, 5))
# lines(density(y), col = 'red')

## ------------------------------------------------------------------------
bd.test(x = x, y = y)

## ------------------------------------------------------------------------
x <- matrix(rnorm(100), nrow = 50, ncol = 2)
y <- matrix(rnorm(100, mean = 3), nrow = 50, ncol = 2)

## ------------------------------------------------------------------------
bd.test(x = x, y = y)

## ------------------------------------------------------------------------
# generate random perturbation:
error <- runif(50, min = -0.3, max = 0.3)
x <- runif(50, 0, 4*pi)
y <- cos(x) + error
# plot(x, y)

## ------------------------------------------------------------------------
bcov.test(x = x, y = y)

## ------------------------------------------------------------------------
x <- matrix(runif(50 * 2, -pi, pi), nrow = 50, ncol = 2)
error <- runif(50, min = -0.1, max = 0.1)
y <- 2 * sin(x[,1] + x[,2]) + error

## ------------------------------------------------------------------------
bcov.test(x = x, y = y, weight = "prob")

## ------------------------------------------------------------------------
# load data:
data("bdvmf")

## ---- eval=FALSE, echo=FALSE---------------------------------------------
#  library(scatterplot3d)
#  scatterplot3d(bdvmf[["x"]], color = bdvmf[["group"]],
#                xlab = "X1", ylab = "X2", zlab = "X3")

## ------------------------------------------------------------------------
# calculate geodesic distance between samples:
Dmat <- nhdist(bdvmf[["x"]], method = "geodesic")
# sample sizes in each group: 150, 150
# Two-Sample Test based on BD :
bd.test(x = Dmat, size = c(150, 150), R = 99, distance = TRUE)

## ------------------------------------------------------------------------
# load data:
data("macaques")
# number of femala and male Macaca fascicularis:
# table(macaques[["group"]])  # f: 9; m: 9
# calculate Riemannian shape distance matrix:
Dmat <- nhdist(macaques[["x"]], method = "riemann")
# hypothesis test with BD:
bd.test(x = Dmat, R = 99, size = c(9, 9), distance = TRUE)

## ------------------------------------------------------------------------
data("ArcticLake")
# Distance matrix between y:
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
# Distance matrix between x:
Dx <- dist(ArcticLake[["depth"]])
# hypothesis test with BCov:
bcov.test(x = Dx, y = Dy, R = 99, distance = TRUE)

## ------------------------------------------------------------------------
n <- 150
bd.test(rnorm(n), size = rep(50, 3))

## ------------------------------------------------------------------------
data("ArcticLake")
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
Dx <- dist(ArcticLake[["depth"]])
# hypothesis test with weighted BCov:
bcov.test(x = Dx, y = Dy, num.permutations = 99, 
          distance = TRUE, weight = "constant")

## ------------------------------------------------------------------------
x <- rnorm(50)
y <- (x > 0) * x + rnorm(50)
z <- (x <= 0) * x + rnorm(50)
example1 <- list(x, y, z)

## ------------------------------------------------------------------------
h <- rnorm(50)
w <- (h)^2
x <- abs(h)
y <- h * (h < 0)
z1 <- h * (h < 0.5)
z2 <- h * (h > -0.5)
z <- cbind(z1, z2)
example2 <- list(w, x, y, z)

## ------------------------------------------------------------------------
bcov.test(x = example1, num.permutations = 199)
bcov.test(x = example2, num.permutations = 199)

## ------------------------------------------------------------------------
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
error <- rnorm(n)
y <- 3*x[, 1] + 5*(x[, 3])^2 + error

## ------------------------------------------------------------------------
res <- bcorsis(y = y, x = x)
head(res[[1]], n = 5)

## ------------------------------------------------------------------------
result <- bcorsis(x = genlung[["covariate"]], 
                  y = genlung[["survival"]], 
                  d = "small", method = "survival")
top_gene <- colnames(genlung[["covariate"]])[result[["ix"]]]
head(top_gene, n = 1)

