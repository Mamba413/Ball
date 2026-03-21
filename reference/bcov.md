# Ball Covariance and Correlation Statistics

Computes Ball Covariance and Ball Correlation statistics, which are
generic dependence measures in Banach spaces.

## Usage

``` r
bcor(x, y, distance = FALSE, weight = FALSE)

bcov(x, y, distance = FALSE, weight = FALSE)
```

## Arguments

- x:

  a numeric vector, matrix, data.frame, or a list containing at least
  two numeric vectors, matrices, or data.frames.

- y:

  a numeric vector, matrix, or data.frame.

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

## Value

- `bcor `:

  Ball Correlation statistic.

- `bcov `:

  Ball Covariance statistic.

## Details

The sample sizes of the two variables must agree, and samples must not
contain missing and infinite values. If we set `distance = TRUE`,
arguments `x`, `y` can be a `dist` object or a symmetric numeric matrix
recording distance between samples; otherwise, these arguments are
treated as data.

`bcov` and `bcor` compute Ball Covariance and Ball Correlation
statistics.

Ball Covariance statistics is a generic dependence measure in Banach
spaces. It enjoys the following properties:

- It is nonnegative and it is equal to zero if and only if variables are
  unassociated;

- It is highly robust;

- It is distribution-free and model-free;

- it is interesting that the HHG is a special case of Ball Covariance
  statistics.

Ball correlation statistics, a normalized version of Ball Covariance
statistics, generalizes Pearson correlation in two fundamental ways:

- It is well-defined for random variables in arbitrary dimension in
  Banach spaces

- BCor is equal to zero implies random variables are unassociated.

The definitions of the Ball Covariance and Ball Correlation statistics
between two random variables are as follows. Suppose, we are given pairs
of independent observations \\\\(x_1, y_1),...,(x_n,y_n)\\\\, where
\\x_i\\ and \\y_i\\ can be of any dimension and the dimensionality of
\\x_i\\ and \\y_i\\ need not be the same. Then, we define sample version
Ball Covariance as: \$\$\mathbf{BCov}\_{\omega, n}^{2}(X,
Y)=\frac{1}{n^{2}}\sum\_{i,j=1}^{n}{(\Delta\_{ij,n}^{X,Y}-\Delta\_{ij,n}^{X}\Delta\_{ij,n}^{Y})^{2}}
\$\$ where: \$\$
\Delta\_{ij,n}^{X,Y}=\frac{1}{n}\sum\_{k=1}^{n}{\delta\_{ij,k}^{X}
\delta\_{ij,k}^{Y}},
\Delta\_{ij,n}^{X}=\frac{1}{n}\sum\_{k=1}^{n}{\delta\_{ij,k}^{X}},
\Delta\_{ij,n}^{Y}=\frac{1}{n}\sum\_{k=1}^{n}{\delta\_{ij,k}^{Y}} \$\$
\$\$\delta\_{ij,k}^{X} = I(x\_{k} \in \bar{B}(x\_{i}, \rho(x\_{i},
x\_{j}))), \delta\_{ij,k}^{Y} = I(y\_{k} \in \bar{B}(y\_{i},
\rho(y\_{i}, y\_{j})))\$\$ Among them, \\\bar{B}(x\_{i}, \rho(x\_{i},
x\_{j}))\\ is a closed ball with center \\x\_{i}\\ and radius
\\\rho(x\_{i}, x\_{j})\\. Similarly, we can define \\
\mathbf{BCov}\_{\omega,n}^2(\mathbf{X},\mathbf{X}) \\ and \\
\mathbf{BCov}\_{\omega,n}^2(\mathbf{Y},\mathbf{Y}) \\. We define Ball
Correlation statistic as follows.
\$\$\mathbf{BCor}\_{\omega,n}^2(\mathbf{X},\mathbf{Y})=
\mathbf{BCov}\_{\omega,n}^2(\mathbf{X},\mathbf{Y})/\sqrt{\mathbf{BCov}\_{\omega,n}^2(\mathbf{X},\mathbf{X})\mathbf{BCov}\_{\omega,n}^2(\mathbf{Y},\mathbf{Y})}
\$\$

We can extend \\\mathbf{BCov}\_{\omega,n}\\ to measure the mutual
independence between \\K\\ random variables:
\$\$\frac{1}{n^{2}}\sum\_{i,j=1}^{n}{\left\[ (\Delta\_{ij,n}^{X\_{1},
...,
X\_{K}}-\prod\_{k=1}^{K}\Delta\_{ij,n}^{X\_{k}})^{2}\prod\_{k=1}^{K}{\hat{\omega}\_{k}(X\_{ki},X\_{kj})}
\right\]}\$\$ where \\X\_{k}(k=1,\ldots,K)\\ are random variables and
\\X\_{ki}\\ is the \\i\\-th observations of \\X\_{k}\\.

See
[`bcov.test`](https://mamba413.github.io/Ball/reference/bcov.test.md)
for a test of independence based on the Ball Covariance statistic.

## References

Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019)
Ball Covariance: A Generic Measure of Dependence in Banach Space,
Journal of the American Statistical Association, DOI:
10.1080/01621459.2018.1543600

Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic
Sure Independence Screening Procedure, Journal of the American
Statistical Association, DOI: 10.1080/01621459.2018.1462709

Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2021). Ball: An R
Package for Detecting Distribution Difference and Association in Metric
Spaces, Journal of Statistical Software, Vol.97(6), doi:
10.18637/jss.v097.i06.

## See also

[`bcov.test`](https://mamba413.github.io/Ball/reference/bcov.test.md),
[`bcorsis`](https://mamba413.github.io/Ball/reference/bcorsis.md)

## Examples

``` r
############# Ball Correlation #############
num <- 50
x <- 1:num
y <- 1:num
bcor(x, y)
#> bcor.constant 
#>             1 
bcor(x, y, weight = "prob")
#> bcor.probability 
#>                1 
bcor(x, y, weight = "chisq")
#> bcor.chisquare 
#>              1 
############# Ball Covariance #############
num <- 50
x <- rnorm(num)
y <- rnorm(num)
bcov(x, y)
#> bcov.constant 
#>  0.0007144357 
bcov(x, y, weight = "prob")
#> bcov.probability 
#>        0.0292523 
bcov(x, y, weight = "chisq")
#> bcov.chisquare 
#>     0.05019143 
```
