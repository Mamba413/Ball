library(Ball)
context("bd and bd.test function")
skip_on_cran()

test_that("Error if computation result for ball divergence is wrong!", {
  target_value <- 2.4032
  names(target_value) <- "bd"
  expect_equal(bd(1:15, size = c(5, 5, 5)), target_value)
  expect_equal(bd.test(1:15, size = c(5, 5, 5), R = 0), target_value)
})


# test abnormal is not necessary. 
# If there exist any abnormal values, .C function will return detailed error message
test_that("Error if input data contain abnormal values", {
  x1 <- rnorm(20)
  x2 <- rnorm(20)
  x1[1] <- NA
  expect_error(bd.test(x1, x2))
  x1[1] <- Inf
  expect_error(bd.test(x1, x2))
})

