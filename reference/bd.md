# Ball Divergence statistic

Compute Ball Divergence statistic, which is a generic dispersion measure
in Banach spaces.

## Usage

``` r
bd(
  x,
  y = NULL,
  distance = FALSE,
  size = NULL,
  num.threads = 1,
  kbd.type = c("sum", "maxsum", "max")
)
```

## Arguments

- x:

  a numeric vector, matrix, data.frame, or a list containing at least
  two numeric vectors, matrices, or data.frames.

- y:

  a numeric vector, matrix, data.frame.

- distance:

  if `distance = TRUE`, the elements of `x` will be considered as a
  distance matrix. Default: `distance = FALSE`.

- size:

  a vector recording sample size of each group.

- num.threads:

  number of threads. If `num.threads = 0`, then all of available cores
  will be used. Default `num.threads = 0`.

- kbd.type:

  a character string specifying the \\K\\-sample Ball Divergence test
  statistic, must be one of `"sum"`, `"summax"`, or `"max"`. Any
  unambiguous substring can be given. Default `kbd.type = "sum"`.

## Value

- `bd `:

  Ball Divergence statistic

## Details

Given the samples not containing missing values, `bd` returns Ball
Divergence statistics. If we set `distance = TRUE`, arguments `x`, `y`
can be a `dist` object or a symmetric numeric matrix recording distance
between samples; otherwise, these arguments are treated as data.

Ball divergence statistic measure the distribution difference of two
datasets in Banach spaces. The Ball divergence statistic is proven to be
zero if and only if two datasets are identical.

The definition of the Ball Divergence statistics is as follows. Given
two independent samples \\ \\x\_{1}, \ldots, x\_{n}\\ \\ with the
associated probability measure \\\mu\\ and \\ \\y\_{1}, \ldots, y\_{m}\\
\\ with \\\nu\\, where the observations in each sample are *i.i.d*. Let
\\\delta(x,y,z)=I(z\in \bar{B}(x, \rho(x,y)))\\, where \\\delta(x,y,z)\\
indicates whether \\z\\ is located in the closed ball \\\bar{B}(x,
\rho(x,y))\\ with center \\x\\ and radius \\\rho(x, y)\\. We denote:
\$\$ A\_{ij}^{X}=\frac{1}{n}\sum\_{u=1}^{n}{\delta(X_i,X_j,X_u)}, \quad
A\_{ij}^{Y}=\frac{1}{m}\sum\_{v=1}^{m}{\delta(X_i,X_j,Y_v)}, \$\$ \$\$
C\_{kl}^{X}=\frac{1}{n}\sum\_{u=1}^{n}{\delta(Y_k,Y_l,X_u)}, \quad
C\_{kl}^{Y}=\frac{1}{m}\sum\_{v=1}^{m}{\delta(Y_k,Y_l,Y_v)}. \$\$
\\A\_{ij}^X\\ represents the proportion of samples \\ \\x\_{1}, \ldots,
x\_{n}\\ \\ located in the ball \\\bar{B}(X_i,\rho(X_i,X_j))\\ and
\\A\_{ij}^Y\\ represents the proportion of samples \\ \\y\_{1}, \ldots,
y\_{m}\\ \\ located in the ball \\\bar{B}(X_i,\rho(X_i,X_j))\\.
Meanwhile, \\C\_{kl}^X\\ and \\C\_{kl}^Y\\ represent the corresponding
proportions located in the ball \\\bar{B}(Y_k,\rho(Y_k,Y_l))\\. The Ball
Divergence statistic is defined as: \$\$D\_{n,m}=A\_{n,m}+C\_{n,m}\$\$

Ball Divergence can be generalized to the *K*-sample test problem.
Suppose we have \\K\\ group samples, each group include \\n\_{k}\\
samples. The definition of \\K\\-sample Ball Divergence statistic could
be to directly sum up the two-sample Ball Divergence statistics of all
sample pairs (`kbd.type = "sum"`) \$\$\sum\_{1 \leq k \< l \leq
K}{D\_{n\_{k},n\_{l}}},\$\$ or to find one sample with the largest
difference to the others (`kbd.type = "maxsum"`)
\$\$\max\_{t}{\sum\_{s=1, s \neq t}^{K}{D\_{n\_{s}, n\_{t}}},}\$\$ to
aggregate the \\K-1\\ most significant different two-sample Ball
Divergence statistics (`kbd.type = "max"`)
\$\$\sum\_{k=1}^{K-1}{D\_{(k)}},\$\$ where \\D\_{(1)}, \ldots,
D\_{(K-1)}\\ are the largest \\K-1\\ two-sample Ball Divergence
statistics among \\\\D\_{n_s, n_t}\| 1 \leq s \< t \leq K\\\\. When
\\K=2\\, the three types of Ball Divergence statistics degenerate into
two-sample Ball Divergence statistic.

See [`bd.test`](https://mamba413.github.io/Ball/reference/bd.test.md)
for a test of distribution equality based on the Ball Divergence.

## References

Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang. Ball Divergence:
Nonparametric two sample test. Ann. Statist. 46 (2018), no. 3,
1109–1137. doi:10.1214/17-AOS1579.
https://projecteuclid.org/euclid.aos/1525313077

## See also

[`bd.test`](https://mamba413.github.io/Ball/reference/bd.test.md)

## Author

Wenliang Pan, Yuan Tian, Xueqin Wang, Heping Zhang

## Examples

``` r
############# Ball Divergence #############
x <- rnorm(50)
y <- rnorm(50)
bd(x, y)
#> bd.constant 
#>  0.00733456 
```
