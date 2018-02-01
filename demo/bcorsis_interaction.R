library(Ball)

# simulation example
p <- 3000
n <- 150
x <- matrix(rnorm(n = n * p), nrow = n, ncol = p)
y <- 4*x[, 1]*x[, 2]

res <- bcorsis(x = x, y = y, method = "interaction")
head(res[["ix"]])
