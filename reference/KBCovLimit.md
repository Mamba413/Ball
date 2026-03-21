# A spectrum-based metric assocaition test (aka, Ball covariance test)

A spectrum-based metric assocaition test (aka, Ball covariance test)

## Usage

``` r
KBCovLimit(x)
```

## Arguments

- x:

  a list containing at least two numeric distance matrices.

## Value

`KBCovLimit` returns a list containing the following components:

- `statistic`:

  Ball Covariance statistic.

- `p.value`:

  the p-value for the test.

## Examples

``` r
library(Ball)
set.seed(1)
K <- 2
x <- list()
Y <- list()
n <- 200
for(i in 1:K){
  Y[[i]] <- rnorm(n)
  x[[i]] <- dist(Y[[i]])
}
kbcov_res <- KBCovLimit(x)[c("statistics", "p.value")]
#> Error in KBCovLimit(x): object 'i' not found
summary(kbcov_res)
#> Error: object 'kbcov_res' not found
bcov.test(Y[[1]],Y[[2]], method = "limit")
#> 
#>  Ball Covariance test of independence (Limit Distribution)
#> 
#> data:  Y[[1]] and Y[[2]]
#> number of observations = 200
#> replicates = 0, weight: constant
#> bcov.constant = 0.00017775, p-value = 0.174
#> alternative hypothesis: random variables are dependent
#> 

set.seed(1)
x <- list()
Y <- list()
n <- 500
K <- 3
for(i in 1:K){
  Y[[i]] <- rnorm(n)
  x[[i]] <- dist(Y[[i]])
}
kbcov_res <- KBCovLimit(x)
#> Error in KBCovLimit(x): object 'i' not found
kbcov_res[c("statistics", "p.value")]
#> Error: object 'kbcov_res' not found
summary(kbcov_res[["eigenvalues"]])
#> Error: object 'kbcov_res' not found
```
