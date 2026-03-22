# Ball Covariance Test

Ball Covariance test of independence. Ball covariance are generic
dependence measures in Banach spaces.

## Usage

``` r
bcov.test(x, ...)

# Default S3 method
bcov.test(
  x,
  y = NULL,
  num.permutations = 99,
  method = c("permutation", "limit"),
  distance = FALSE,
  weight = FALSE,
  seed = 1,
  num.threads = 0,
  ...
)

# S3 method for class 'formula'
bcov.test(formula, data, subset, na.action, ...)
```

## Arguments

- x:

  a numeric vector, matrix, data.frame, or a list containing at least
  two numeric vectors, matrices, or data.frames.

- ...:

  further arguments to be passed to or from methods.

- y:

  a numeric vector, matrix, or data.frame.

- num.permutations:

  the number of permutation replications. When `num.permutations = 0`,
  the function just returns the Ball Covariance statistic. Default:
  `num.permutations = 99`.

- method:

  if `method = "permutation"`, a permutation procedure is carried out to
  compute the \\p\\-value; if ` method = "limit"`, an approximate null
  distribution is used when `weight = "constant"`. Any unambiguous
  substring can be given. Default `method = "permutation"`.

- distance:

  if `distance = TRUE`, the elements of `x` and `y` are considered as
  distance matrices.

- weight:

  a logical or character string used to choose the weight form of Ball
  Covariance statistic.. If input is a character string, it must be one
  of `"constant"`, `"probability"`, `"chisquare"`, `"rbf"`. Any
  unambiguous substring can be given. If input is a logical value, it is
  equivalent to `weight = "probability"` if `weight = TRUE` while
  equivalent to `weight = "constant"` if `weight = FALSE`. Default:
  `weight = FALSE`.

- seed:

  the random seed. Default `seed = 1`.

- num.threads:

  number of threads. If `num.threads = 0`, then all of available cores
  will be used. Default `num.threads = 0`.

- formula:

  a formula of the form `~ u + v`, where each of `u` and `v` are numeric
  variables giving the data values for one sample. The samples must be
  of the same length.

- data:

  an optional matrix or data frame (or similar: see `model.frame`)
  containing the variables in the formula formula. By default the
  variables are taken from environment(formula).

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

## Value

If `num.permutations > 0`, `bcov.test` returns a `htest` class object
containing the following components:

- `statistic`:

  Ball Covariance statistic.

- `p.value`:

  the p-value for the test.

- `replicates`:

  permutation replications of the test statistic.

- `size`:

  sample size.

- `complete.info`:

  a `list` mainly containing two vectors, the first vector is the Ball
  Covariance statistics with different weights, the second is the
  \\p\\-values of weighted Ball Covariance tests.

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating what type of test was performed.

- `data.name`:

  description of data.

If `num.permutations = 0`, `bcov.test` returns a statistic value.

## Details

`bcov.test` is non-parametric tests of independence in Banach spaces. It
can detect the dependence between two random objects (variables) and the
mutual dependence among at least three random objects (variables).

If two samples are pass to arguments `x` and `y`, the sample sizes (i.e.
number of rows or length of the vector) of the two variables must agree.
If a [`list`](https://rdrr.io/r/base/list.html) object is passed to `x`,
this `list` must contain at least two numeric vectors, matrices, or
data.frames, and each element of this `list` must with the same sample
size. Moreover, data pass to `x` or `y` must not contain missing or
infinite values. If `distance = TRUE`, `x` is considered as a distance
matrix or a list containing distance matrices, and `y` is considered as
a distance matrix; otherwise, these arguments are treated as data.

`bcov.test` utilizes the Ball Covariance statistics (see
[`bcov`](https://mamba413.github.io/Ball/reference/bcov.md)) to measure
dependence and derives a \\p\\-value via replicating the random
permutation `num.permutations` times.

See Pan et al 2018 for theoretical properties of the test, including
statistical consistency.

## Note

Actually, `bcov.test` simultaneously computing Ball Covariance
statistics with `"constant"`, `"probability"`, `"chisquare"`, `"rbf"`
weights. Users can get other Ball Covariance statistics with different
weight and their corresponding \\p\\-values in the `complete.info`
element of output. We give a quick example below to illustrate.

## References

Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019)
Ball Covariance: A Generic Measure of Dependence in Banach Space,
Journal of the American Statistical Association, DOI:
10.1080/01621459.2018.1543600

Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2021). Ball: An R
Package for Detecting Distribution Difference and Association in Metric
Spaces, Journal of Statistical Software, Vol.97(6), doi:
10.18637/jss.v097.i06.

## See also

[`bcov`](https://mamba413.github.io/Ball/reference/bcov.md),
[`bcor`](https://mamba413.github.io/Ball/reference/bcov.md)

## Author

Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu, Jin Zhu

## Examples

``` r
set.seed(1)

################# Quick Start #################
noise <- runif(50, min = -0.3, max = 0.3)
x <- runif(50, 0, 4*pi)
y <- cos(x) + noise
# plot(x, y)
res <- bcov.test(x, y)
res
#> 
#>  Ball Covariance test of independence (Permutation)
#> 
#> data:  x and y
#> number of observations = 50
#> replicates = 99, weight: constant
#> bcov.constant = 0.0021965, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 
## get all Ball Covariance statistics:
res[["complete.info"]][["statistic"]]
#>    bcov.constant bcov.probability   bcov.chisquare         bcov.rbf 
#>      0.002196459      0.052685015      0.102469543      0.000000000 
## get all test result:
res[["complete.info"]][["p.value"]]
#>    bcov.constant.pvalue bcov.probability.pvalue   bcov.chisquare.pvalue 
#>                    0.01                    0.01                    0.01 
#>         bcov.rbf.pvalue 
#>                    1.00 

################# Quick Start #################
x <- matrix(runif(50 * 2, -pi, pi), nrow = 50, ncol = 2)
noise <- runif(50, min = -0.1, max = 0.1)
y <- sin(x[,1] + x[,2]) + noise
bcov.test(x = x, y = y, weight = "prob")
#> 
#>  Ball Covariance test of independence (Permutation)
#> 
#> data:  x and y
#> number of observations = 50
#> replicates = 99, weight: probability
#> bcov.probability = 0.035357, p-value = 0.1
#> alternative hypothesis: random variables are dependent
#> 

################# Ball Covariance Test for Non-Hilbert Data #################
# load data:
data("ArcticLake")
# Distance matrix between y:
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
# Distance matrix between x:
Dx <- dist(ArcticLake[["depth"]])
# hypothesis test with BCov:
bcov.test(x = Dx, y = Dy, distance = TRUE)
#> 
#>  Ball Covariance test of independence (Permutation)
#> 
#> data:  Dx and Dy
#> number of observations = 39
#> replicates = 99, weight: constant
#> bcov.constant = 0.0083848, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 

################  Weighted Ball Covariance Test  #################
data("ArcticLake")
Dy <- nhdist(ArcticLake[["x"]], method = "compositional")
Dx <- dist(ArcticLake[["depth"]])
# hypothesis test with weighted BCov:
bcov.test(x = Dx, y = Dy, distance = TRUE, weight = "prob")
#> 
#>  Ball Covariance test of independence (Permutation)
#> 
#> data:  Dx and Dy
#> number of observations = 39
#> replicates = 99, weight: probability
#> bcov.probability = 0.088597, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 

################# Mutual Independence Test #################
x <- rnorm(50)
y <- (x > 0) * x + rnorm(50)
z <- (x <= 0) * x + rnorm(50)
data_list <- list(x, y, z)
bcov.test(data_list)
#> 
#>  Ball Covariance test of mutual independence (Permutation)
#> 
#> data:  data_list
#> number of observations = 50
#> replicates = 99, weight: constant
#> bcov.constant = 0.0018673, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 
data_list <- lapply(data_list, function(x) {
  as.matrix(dist(x))
})
bcov.test(data_list, distance = TRUE)
#> 
#>  Ball Covariance test of mutual independence (Permutation)
#> 
#> data:  data_list
#> number of observations = 50
#> replicates = 99, weight: constant
#> bcov.constant = 0.0018673, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 
bcov.test(data_list, distance = FALSE, weight = "chi")
#> 
#>  Ball Covariance test of mutual independence (Permutation)
#> 
#> data:  data_list
#> number of observations = 50
#> replicates = 99, weight: chisquare
#> bcov.chisquare = 1.7707, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 

################# Mutual Independence Test for Meteorology data #################
data("meteorology")
bcov.test(meteorology)
#> 
#>  Ball Covariance test of mutual independence (Permutation)
#> 
#> data:  meteorology
#> number of observations = 46
#> replicates = 99, weight: constant
#> bcov.constant = 0.0064547, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 

################  Testing via approximate limit distribution  #################
if (FALSE) { # \dontrun{
set.seed(1)
n <- 2000
x <- rnorm(n)
y <- rnorm(n)
bcov.test(x, y, method = "limit")
bcov.test(x, y)
} # }

################  Formula interface  ################
## independence test:
bcov.test(~ CONT + INTG, data = USJudgeRatings)
#> 
#>  Ball Covariance test of independence (Permutation)
#> 
#> data:  CONT and INTG
#> number of observations = 43
#> replicates = 99, weight: constant
#> bcov.constant = 0.00079459, p-value = 0.39
#> alternative hypothesis: random variables are dependent
#> 
## independence test with chisquare weight:
bcov.test(~ CONT + INTG, data = USJudgeRatings, weight = "chi")
#> 
#>  Ball Covariance test of independence (Permutation)
#> 
#> data:  CONT and INTG
#> number of observations = 43
#> replicates = 99, weight: chisquare
#> bcov.chisquare = 0.038503, p-value = 0.39
#> alternative hypothesis: random variables are dependent
#> 
## mutual independence test:
bcov.test(~ CONT + INTG + DMNR, data = USJudgeRatings)
#> 
#>  Ball Covariance test of mutual independence (Permutation)
#> 
#> data:  CONT and INTG and DMNR
#> number of observations = 43
#> replicates = 99, weight: constant
#> bcov.constant = 0.0078353, p-value = 0.01
#> alternative hypothesis: random variables are dependent
#> 
```
