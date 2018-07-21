library(Ball)

data("ArcticLake")
# Distance matrix between y:
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
# Distance matrix between x:
Dx <- dist(ArcticLake[["depth"]])

# hypothesis test for BCov with probability weight:
bcov.test(x = Dx, y = Dy, R = 99, dst = TRUE, weight = TRUE)
# hypothesis test for BCov with Chi-square weight:
bcov.test(x = Dx, y = Dy, R = 99, dst = TRUE, weight = "chisq")
