library(Ball)

x1 <- rnorm(n = 50)
x2 <- rt(n = 50, df = 3)
x3 <- rnorm(n = 50, sd = 2)
x <- c(x1, x2, x3)

bd.test(x, size = c(50, 50, 50))
