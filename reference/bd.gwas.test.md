# Fast K-sample Ball Divergence Test for GWAS Data

Fast K-sample Ball Divergence Test for GWAS Data

## Usage

``` r
bd.gwas.test(
  x,
  snp,
  screening.method = c("permute", "spectrum"),
  refine = TRUE,
  num.permutations,
  distance = FALSE,
  alpha,
  screening.result = NULL,
  verbose = TRUE,
  seed = 1,
  num.threads = 0,
  ...
)
```

## Arguments

- x:

  a numeric vector, matrix, data.frame, dist object.

- snp:

  a numeric matrix recording the values of single nucleotide
  polymorphism (SNP). Each column must be an integer vector.

- screening.method:

  if `screening.method = "spectrum"`, the spectrum method is applied to
  screening the candidate SNPs, or otherwise, the permutation method is
  applied. Default: `screening.method = "permute"`.

- refine:

  a logical value. If `refine = TRUE`, a \\p\\-values refining process
  is applied to the SNPs which passes the pre-screening process.
  Default: `refine = TRUE` (At present, `refine = FALSE` is not
  available).

- num.permutations:

  the number of permutation replications. When `num.permutations = 0`,
  the function just returns the Ball Divergence statistic. Default:
  `num.permutations = 100 * ncol(snp)`

- distance:

  if `distance = TRUE`, the elements of `x` will be considered as a
  distance matrix. Default: `distance = FALSE`.

- alpha:

  the significance level. Default: `0.05 / ncol(snp)`.

- screening.result:

  A object return by `bd.gwas.test` that preserving the pre-screening
  result. It works only if the pre-screening is available. Default:
  `screening.result = NULL`.

- verbose:

  Show computation status and estimated runtimes. Default:
  `verbose = FALSE`.

- seed:

  the random seed. Default `seed = 1`.

- num.threads:

  number of threads. If `num.threads = 0`, then all of available cores
  will be used. Default `num.threads = 0`.

- ...:

  further arguments to be passed to or from methods.

## Value

bd.gwas.test returns a list containing the following components:

- `statistic`:

  ball divergence statistics vector.

- `permuted.statistic`:

  a data.frame containing permuted ball divergence statistic for
  pre-screening SNPs. If `refine = FALSE`, it takes value `NULL`.

- `eigenvalue`:

  the eigenvalue of spectrum decomposition. If `refine = TRUE`, it takes
  value `NULL`.

- `p.value`:

  the p-values of ball divergence test.

- `refined.snp`:

  the SNPs have been refined.

- `refined.p.value`:

  the refined \\p\\-value of significant snp.

- `refined.permuted.statistic`:

  a data.frame containing permuted ball divergence statistics for
  refining \\p\\-values.

- `screening.result`:

  a list containing the result of screening.

## References

Yue Hu, Haizhu Tan, Cai Li, and Heping Zhang. (2021). Identifying
genetic risk variants associated with brain volumetric phenotypes via
K-sample Ball Divergence method. Genetic Epidemiology, 1–11.
https://doi.org/10.1002/gepi.22423

## See also

[`bd`](https://mamba413.github.io/Ball/reference/bd.md),
[`bd.test`](https://mamba413.github.io/Ball/reference/bd.test.md)

## Author

Jin Zhu

## Examples

``` r
# \donttest{
library(Ball)
set.seed(1234)
num <- 200
snp_num <- 500
p <- 5
x <- matrix(rnorm(num * p), nrow = num)
snp <- sapply(1:snp_num, function(i) {
  sample(0:2, size = num, replace = TRUE)
})
snp1 <- sapply(1:snp_num, function(i) {
  sample(1:2, size = num, replace = TRUE)
})
snp <- cbind(snp, snp1)
res <- Ball::bd.gwas.test(x = x, snp = snp)
#> =========== Pre-screening SNPs ===========
#> [1] "None of SNP pass the pre-screening process!"
mean(res[["p.value"]] < 0.05)
#> [1] 0.039
mean(res[["p.value"]] < 0.005)
#> [1] 0.005

## only return the test statistics;
res <- Ball::bd.gwas.test(x = x, snp = snp, num.permutation = 0)

## save pre-screening process results:
x <- matrix(rnorm(num * p), nrow = num)
snp <- sapply(1:snp_num, function(i) {
  sample(0:2, size = num, replace = TRUE, prob = c(1/2, 1/4, 1/4))
})
snp_screening <- Ball::bd.gwas.test(x = x, snp = snp,
                                    alpha = 5*10^-4, 
                                    num.permutations = 19999)
#> =========== Pre-screening SNPs ===========
#> [1] "None of SNP pass the pre-screening process!"
mean(res[["p.value"]] < 0.05)
#> [1] 0
mean(res[["p.value"]] < 0.005)
#> [1] 0
mean(res[["p.value"]] < 0.0005)
#> [1] 0
## refine p-value according to the pre-screening process result:
res <- Ball::bd.gwas.test(x = x, snp = snp, alpha = 5*10^-4,
                          num.permutations = 19999,
                          screening.result = snp_screening[["screening.result"]])
# }
```
