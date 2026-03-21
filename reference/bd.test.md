# Ball Divergence based Equality of Distributions Test

Performs the nonparametric two-sample or \\K\\-sample Ball Divergence
test for equality of multivariate distributions

## Usage

``` r
bd.test(x, ...)

# Default S3 method
bd.test(
  x,
  y = NULL,
  num.permutations = 99,
  method = c("permutation", "limit"),
  distance = FALSE,
  size = NULL,
  seed = 1,
  num.threads = 0,
  kbd.type = c("sum", "maxsum", "max"),
  weight = c("constant", "variance", "rbf"),
  ...
)

# S3 method for class 'formula'
bd.test(formula, data, subset, na.action, ...)
```

## Arguments

- x:

  a numeric vector, matrix, data.frame, or a list containing at least
  two numeric vectors, matrices, or data.frames.

- ...:

  further arguments to be passed to or from methods.

- y:

  a numeric vector, matrix, data.frame.

- num.permutations:

  the number of permutation replications. When `num.permutations = 0`,
  the function just returns the Ball Divergence statistic. Default:
  `num.permutations = 99`.

- method:

  if `method = "permutation"`, a permutation procedure is carried out to
  compute the \\p\\-value; if ` method = "limit"`, an approximate null
  distribution is used when `weight = "constant"`. Any unambiguous
  substring can be given. Default `method = "permutation"`.

- distance:

  if `distance = TRUE`, the elements of `x` will be considered as a
  distance matrix. Default: `distance = FALSE`.

- size:

  a vector recording sample size of each group.

- seed:

  the random seed. Default `seed = 1`.

- num.threads:

  number of threads. If `num.threads = 0`, then all of available cores
  will be used. Default `num.threads = 0`.

- kbd.type:

  a character string specifying the \\K\\-sample Ball Divergence test
  statistic, must be one of `"sum"`, `"summax"`, or `"max"`. Any
  unambiguous substring can be given. Default `kbd.type = "sum"`.

- weight:

  a character string specifying the weight form of Ball Divergence
  statistic. It must be one of `"constant"` or `"variance"`. Any
  unambiguous substring can be given. Default: `weight = "constant"`.

- formula:

  a formula of the form `response ~ group` where `response` gives the
  data values and `group` a vector or factor of the corresponding
  groups.

- data:

  an optional matrix or data frame (or similar: see `model.frame`)
  containing the variables in the formula `formula`. By default the
  variables are taken from `environment(formula)`.

- subset:

  an optional vector specifying a subset of observations to be used.

- na.action:

  a function which indicates what should happen when the data contain
  `NA`s. Defaults to `getOption("na.action")`.

## Value

If `num.permutations > 0`, `bd.test` returns a `htest` class object
containing the following components:

- `statistic`:

  Ball Divergence statistic.

- `p.value`:

  the \\p\\-value for the test.

- `replicates`:

  permutation replications of the test statistic.

- `size`:

  sample sizes.

- `complete.info`:

  a `list` mainly containing two vectors, the first vector is the Ball
  Divergence statistics with different aggregation strategy and weight,
  the second vector is the \\p\\-values of tests.

- `alternative`:

  a character string describing the alternative hypothesis.

- `method`:

  a character string indicating what type of test was performed.

- `data.name`:

  description of data.

If `num.permutations = 0`, `bd.test` returns a statistic value.

## Details

`bd.test` is nonparametric test for the two-sample or \\K\\-sample
problem. It can detect distribution difference between \\K(K \geq 2)\\
sample even though sample size are imbalanced. This test can cope well
multivariate dataset or complex dataset.

If only `x` is given, the statistic is computed from the original pooled
samples, stacked in matrix where each row is a multivariate observation,
or from the distance matrix when `distance = TRUE`. The first `sizes[1]`
rows of `x` are the first sample, the next `sizes[2]` rows of `x` are
the second sample, etc. If `x` is a `list`, its elements are taken as
the samples to be compared, and hence, this `list` must contain at least
two numeric data vectors, matrices or data.frames.

`bd.test` utilizes the Ball Divergence statistics (see
[`bd`](https://mamba413.github.io/Ball/reference/bd.md)) to measure
dispersion and derives a \\p\\-value via replicating the random
permutation `num.permutations` times. The function simply returns the
test statistic when `num.permutations = 0`.

The time complexity of `bd.test` is around \\O(R \times n^2)\\, where
\\R\\ = `num.permutations` and \\n\\ is sample size.

## Note

Actually, `bd.test` simultaneously computing `"sum"`, `"summax"`, and
`"max"` Ball Divergence statistics when \\K \geq 3\\. Users can get
other Ball Divergence statistics and their corresponding \\p\\-values in
the `complete.info` element of output. We give a quick example below to
illustrate.

## References

Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. Ball Divergence:
Nonparametric two sample test. Annals of Statistics. 46 (2018), no. 3,
1109–1137. doi:10.1214/17-AOS1579.
https://projecteuclid.org/euclid.aos/1525313077

Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2021). Ball: An R
Package for Detecting Distribution Difference and Association in Metric
Spaces, Journal of Statistical Software, Vol.97(6), doi:
10.18637/jss.v097.i06.

## See also

[`bd`](https://mamba413.github.io/Ball/reference/bd.md)

## Author

Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang, Jin Zhu

## Examples

``` r
################# Quick Start #################
set.seed(1)
x <- rnorm(50)
y <- rnorm(50, mean = 1)
# plot(density(x))
# lines(density(y), col = "red")
bd.test(x = x, y = y)
#> 
#>  2-sample Ball Divergence Test (Permutation)
#> 
#> data:  x and y 
#> number of observations = 100, group sizes: 50 50
#> replicates = 99, weight: constant
#> bd.constant = 0.092215, p-value = 0.01
#> alternative hypothesis: distributions of samples are distinct
#> 

################# Quick Start #################
x <- matrix(rnorm(100), nrow = 50, ncol = 2)
y <- matrix(rnorm(100, mean = 3), nrow = 50, ncol = 2)
# Hypothesis test with Standard Ball Divergence:
bd.test(x = x, y = y)
#> 
#>  2-sample Ball Divergence Test (Permutation)
#> 
#> data:  x and y 
#> number of observations = 100, group sizes: 50 50
#> replicates = 99, weight: constant
#> bd.constant = 0.5653, p-value = 0.01
#> alternative hypothesis: distributions of samples are distinct
#> 

################# Simlated Non-Hilbert data #################
data("bdvmf")
if (FALSE) { # \dontrun{
library(scatterplot3d)
scatterplot3d(bdvmf[["x"]], color = bdvmf[["group"]], 
              xlab = "X1", ylab = "X2", zlab = "X3")
} # }
# calculate geodesic distance between sample:
Dmat <- nhdist(bdvmf[["x"]], method = "geodesic")
# hypothesis test with BD :
bd.test(x = Dmat, size = c(150, 150), num.permutations = 99, distance = TRUE)
#> 
#>  2-sample Ball Divergence Test (Permutation)
#> 
#> data:  Dmat 
#> number of observations = 300, group sizes: 150 150
#> replicates = 99, weight: constant
#> bd.constant = 0.14483, p-value = 0.01
#> alternative hypothesis: distributions of samples are distinct
#> 

################# Non-Hilbert Real Data #################
# load data:
data("macaques")
# number of femala and male Macaca fascicularis:
table(macaques[["group"]])
#> 
#> f m 
#> 9 9 
# calculate Riemannian shape distance matrix:
Dmat <- nhdist(macaques[["x"]], method = "riemann")
# hypothesis test with BD:
bd.test(x = Dmat, num.permutations = 99, size = c(9, 9), distance = TRUE)
#> 
#>  2-sample Ball Divergence Test (Permutation)
#> 
#> data:  Dmat 
#> number of observations = 18, group sizes: 9 9
#> replicates = 99, weight: constant
#> bd.constant = 0.1922, p-value = 0.03
#> alternative hypothesis: distributions of samples are distinct
#> 

################  K-sample Test  #################
n <- 150
bd.test(rnorm(n), size = c(40, 50, 60))
#> 
#>  3-sample Ball Divergence Test (Permutation)
#> 
#> data:  rnorm(n) 
#> number of observations = 150, group sizes: 40 50 60
#> replicates = 99, weight: constant, kbd.type: sum
#> kbd.sum.constant = 0.056299, p-value = 0.24
#> alternative hypothesis: distributions of samples are distinct
#> 
# alternative input method:
x <- lapply(c(40, 50, 60), rnorm)
res <- bd.test(x)
res
#> 
#>  3-sample Ball Divergence Test (Permutation)
#> 
#> data:  x 
#> number of observations = 150, group sizes: 40 50 60
#> replicates = 99, weight: constant, kbd.type: sum
#> kbd.sum.constant = 0.021721, p-value = 0.87
#> alternative hypothesis: distributions of samples are distinct
#> 
## get all Ball Divergence statistics:
res[["complete.info"]][["statistic"]]
#>    kbd.sum.constant    kbd.sum.variance    kbd.max.constant    kbd.max.variance 
#>          0.02172147          0.02172147          0.01584354          0.01584354 
#> kbd.maxsum.constant kbd.maxsum.variance 
#>          0.01584354          0.01584354 
## get all test result:
res[["complete.info"]][["p.value"]]
#>    kbd.sum.constant.pvalue    kbd.sum.variance.pvalue 
#>                       0.87                       0.87 
#>    kbd.max.constant.pvalue    kbd.max.variance.pvalue 
#>                       0.89                       0.89 
#> kbd.maxsum.constant.pvalue kbd.maxsum.variance.pvalue 
#>                       0.89                       0.89 

################  Testing via approximate limit distribution  #################
if (FALSE) { # \dontrun{
set.seed(1)
n <- 1000
x <- rnorm(n)
y <- rnorm(n)
res <- bd.test(x, y, method = "limit")
bd.test(x, y)
} # }

################  Formula interface  ################
## Two-sample test
bd.test(extra ~ group, data = sleep)
#> 
#>  2-sample Ball Divergence Test (Permutation)
#> 
#> data:  extra by group
#> number of observations = 20, group sizes: 10 10
#> replicates = 99, weight: constant
#> bd.constant = 0.0895, p-value = 0.32
#> alternative hypothesis: distributions of samples are distinct
#> 
## K-sample test
bd.test(Sepal.Width ~ Species, data = iris)
#> 
#>  3-sample Ball Divergence Test (Permutation)
#> 
#> data:  Sepal.Width by Species
#> number of observations = 150, group sizes: 50 50 50
#> replicates = 99, weight: constant, kbd.type: sum
#> kbd.sum.constant = 0.48817, p-value = 0.01
#> alternative hypothesis: distributions of samples are distinct
#> 
bd.test(Sepal.Width ~ Species, data = iris, kbd.type = "max")
#> 
#>  3-sample Ball Divergence Test (Permutation)
#> 
#> data:  Sepal.Width by Species
#> number of observations = 150, group sizes: 50 50 50
#> replicates = 99, weight: constant, kbd.type: max
#> kbd.max.constant = 0.4618, p-value = 0.01
#> alternative hypothesis: distributions of samples are distinct
#> 
```
