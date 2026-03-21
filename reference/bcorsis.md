# Ball Correlation based Sure Independence Screening (BCor-SIS)

Generic non-parametric sure independence screening (SIS) procedure based
on Ball Correlation. Ball correlation is a generic measure of dependence
in Banach spaces.

## Usage

``` r
bcorsis(
  x,
  y,
  d = "small",
  weight = c("constant", "probability", "chisquare"),
  method = "standard",
  distance = FALSE,
  category = FALSE,
  parms = list(d1 = 5, d2 = 5, df = 3),
  num.threads = 0
)
```

## Arguments

- x:

  a numeric matrix or data.frame included \\n\\ rows and \\p\\ columns.
  Each row is an observation vector and each column corresponding to a
  explanatory variable, generally \\p \>\> n\\.

- y:

  a numeric vector, matrix, or data.frame.

- d:

  the hard cutoff rule suggests selecting \\d\\ variables. Setting
  `d = "large"` or `d = "small"` means `n - 1` or `floor(n/log(n))`
  variables are selected. If `d` is a integer, `d` variables are
  selected. Default: `d = "small"`.

- weight:

  a logical or character string used to choose the weight form of Ball
  Covariance statistic.. If input is a character string, it must be one
  of `"constant"`, `"probability"`, `"chisquare"`, `"rbf"`. Any
  unambiguous substring can be given. If input is a logical value, it is
  equivalent to `weight = "probability"` if `weight = TRUE` while
  equivalent to `weight = "constant"` if `weight = FALSE`. Default:
  `weight = FALSE`.

- method:

  specific method for the BCor-SIS procedure. It must be one of
  `"standard"`, `"lm"`, `"gam"`, `"interaction"`, or `"survival"`.
  Setting `method = "standard"` means performing standard SIS procedure
  while the options `"lm"` and `"gam"` mean carrying out iterative SIS
  procedure with ordinary linear regression and generalized additive
  models, respectively. The options `"interaction"` and `"survival"` are
  designed for detecting variables with potential linear interaction and
  associated with left censored responses, respectively. Any unambiguous
  substring can be given. Default: `method = "standard"`.

- distance:

  if `distance = TRUE`, `y` will be considered as a distance matrix.
  Arguments only available when `method = "standard"` and
  `method = "interaction"`. Default: `distance = FALSE`.

- category:

  a logical value or integer vector indicating columns to be selected as
  categorical variables. If `category` is an integer vector, the
  positive/negative integers select/discard the corresponding columns;
  If `category` is a logical value, `category = TRUE` select all
  columns, `category = FALSE` select none column. Default:
  `category = FALSE`.

- parms:

  parameters list only available when `method = "lm"` or `"gam"`. It
  contains three parameters: `d1`, `d2`, and `df`. `d1` is the number of
  initially selected variables, `d2` is the number of variables added in
  each iteration. `df` is a degree freedom of basis in generalized
  additive models playing a role only when `method = "gam"`. Default:
  `parms = list(d1 = 5, d2 = 5, df = 3)`.

- num.threads:

  number of threads. If `num.threads = 0`, then all of available cores
  will be used. Default `num.threads = 0`.

## Value

- `ix `:

  the indices vector corresponding to variables selected by BCor-SIS.

- `method `:

  the method used.

- `weight `:

  the weight used.

- `complete.info `:

  a `list` mainly containing a \\p \times 3\\ matrix, where each row is
  a variable and each column is a weight Ball Correlation statistic. If
  `method = "gam"` or `method = "lm"`, `complete.info` is an empty list.

## Details

`bcorsis` performs a model-free generic sure independence screening
procedure, BCor-SIS, to pick out variables from `x` which are
potentially associated with `y`. BCor-SIS relies on Ball correlation, a
universal dependence measure in Banach spaces. Ball correlation (BCor)
ranges from 0 to 1. A larger BCor implies they are likely to be
associated while Bcor is equal to 0 implies they are unassociated. (See
[`bcor`](https://mamba413.github.io/Ball/reference/bcov.md) for
details.) Consequently, BCor-SIS pick out variables with larger Bcor
values with `y`.

Theory and numerical result indicate that BCor-SIS has following
advantages:

- BCor-SIS can retain the efficient variables even when the
  dimensionality (i.e., `ncol(x)`) is an exponential order of the sample
  size (i.e., `exp(nrow(x))`);

- It is distribution-free and model-free;

- It is very robust;

- It is works well for complex data, such as shape and survival data;

If `x` is a matrix, the sample sizes of `x` and `y` must agree. If `x`
is a [`list`](https://rdrr.io/r/base/list.html) object, each element of
this `list` must with the same sample size. `x` and `y` must not contain
missing or infinite values.

When `method = "survival"`, the matrix or data.frame pass to `y` must
have exactly two columns, where the first column is event (failure) time
while the second column is a dichotomous censored status.

## Note

`bcorsis` simultaneously computing Ball Correlation statistics with
`"constant"`, `"probability"`, and `"chisquare"` weights. Users can get
other Ball Correlation statistics with different weight in the
`complete.info` element of output. We give a quick example below to
illustrate.

## References

Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) A Generic
Sure Independence Screening Procedure, Journal of the American
Statistical Association, DOI: 10.1080/01621459.2018.1462709

## See also

[`bcor`](https://mamba413.github.io/Ball/reference/bcov.md)

## Author

Wenliang Pan, Weinan Xiao, Xueqin Wang, Hongtu Zhu, Jin Zhu

## Examples

``` r
if (FALSE) { # \dontrun{

############### Quick Start for bcorsis function ###############
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n)
y <- 3 * x[, 1] + 5 * (x[, 3])^2 + eps
res <- bcorsis(y = y, x = x)
head(res[["ix"]])
head(res[["complete.info"]][["statistic"]])

############### BCor-SIS: Censored Data Example ###############
data("genlung")
result <- bcorsis(x = genlung[["covariate"]], y = genlung[["survival"]], 
                  method = "survival")
index <- result[["ix"]]
top_gene <- colnames(genlung[["covariate"]])[index]
head(top_gene, n = 1)


############### BCor-SIS: Interaction Pursuing ###############
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n)
y <- 3 * x[, 1] * x[, 5] * x[, 10] + eps
res <- bcorsis(y = y, x = x, method = "interaction")
head(res[["ix"]])

############### BCor-SIS: Iterative Method ###############
library(mvtnorm)
set.seed(1)
n <- 150
p <- 3000
sigma_mat <- matrix(0.5, nrow = p, ncol = p)
diag(sigma_mat) <- 1
x <- rmvnorm(n = n, sigma = sigma_mat)
eps <- rnorm(n)
rm(sigma_mat); gc(reset = TRUE)
y <- 3 * (x[, 1])^2 + 5 * (x[, 2])^2 + 5 * x[, 8] - 8 * x[, 16] + eps
res <- bcorsis(y = y, x = x, method = "lm", d = 15)
res <- bcorsis(y = y, x = x, method = "gam", d = 15)
res[["ix"]]

############### Weighted BCor-SIS: Probability weight ###############
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
eps <- rnorm(n)
y <- 3 * x[, 1] + 5 * (x[, 3])^2 + eps
res <- bcorsis(y = y, x = x, weight = "prob")
head(res[["ix"]])
# Alternative, chisq weight:
res <- bcorsis(y = y, x = x, weight = "chisq")
head(res[["ix"]])

############### BCor-SIS: GWAS data ###############
set.seed(1)
n <- 150
p <- 3000
x <- sapply(1:p, function(i) {
  sample(0:2, size = n, replace = TRUE)
})
eps <- rnorm(n)
y <- 6 * x[, 1] - 7 * x[, 2] + 5 * x[, 3] + eps
res <- bcorsis(x = x, y = y, category = TRUE)
head(res[["ix"]])
head(res[["complete.info"]][["statistic"]])

x <- cbind(matrix(rnorm(n * 2), ncol = 2), x)
# remove the first two columns:
res <- bcorsis(x = x, y = y, category = c(-1, -2))
head(res[["ix"]])

x <- cbind(x[, 3:5], matrix(rnorm(n * p), ncol = p))
res <- bcorsis(x = x, y = y, category = 1:3)
head(res[["ix"]], n = 10)
} # }
```
