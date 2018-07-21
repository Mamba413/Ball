library(Ball)

######## parallel ######## 
set.seed(1)
x <- matrix(rnorm(200 * 2), nrow = 200)
y <- matrix(rnorm(200 * 2), nrow = 200)

t1 <- Sys.time()
bcov.test(x, y, R = 1999, num.threads = 1)
Sys.time() - t1

t1 <- Sys.time()
bcov.test(x, y, R = 1999, num.threads = 2)
Sys.time() - t1

