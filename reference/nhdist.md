# Distance Matrix Computation for Non-Hilbert Data

This function computes and returns the numeric distance matrix computed
by using the specified distance measure to compute the distances between
the rows of a data matrix.

## Usage

``` r
nhdist(x, method = "geodesic")
```

## Arguments

- x:

  a numeric matrix, data frame or numeric array of dimension \\k \times
  m \times n\\ containing \\n\\ samples in \\k \times m\\ dimension.

- method:

  the distance measure to be used. This must be one of `"geodesic"`,
  `"compositional"`, or `"riemann"`. Any unambiguous substring can be
  given.

## Value

\\n \times n\\ numeric distance matrix

## Details

Available distance measures are geodesic, compositional and riemann.
Denoting any two sample in the dataset as \\x\\ and \\y\\, we give the
definition of distance measures as follows.

geodesic:

The shortest route between two points on the Earth's surface, namely, a
segment of a great circle. \$\$\arccos(x^{T}y), \\x\\\_{2} = \\y\\\_{2}
= 1\$\$

compositional:

First, we apply scale transformation to it, i.e., \\(x\_{i1}/t, ...,
x\_{ip}/t\_{i}), t\_{i} = \sum\_{d=1}^{p}{x\_{d}}\\ . Then, apply the
square root transformation to data and calculate the geodesic distance
between samples.

riemann:

\\k \times m \times n\\ array where \\k\\ = number of landmarks, \\m\\ =
number of dimensions and \\n\\ = sample size. Detail about riemannian
shape distance was given in Kendall, D. G. (1984).

## References

Kendall, D. G. (1984). Shape manifolds, Procrustean metrics and complex
projective spaces, Bulletin of the London Mathematical Society, 16,
81-121.

## Examples

``` r
data('bdvmf')
Dmat <- nhdist(bdvmf[['x']], method = "geodesic")

data("ArcticLake")
Dmat <- nhdist(ArcticLake[['x']], method = "compositional")

data("macaques")
Dmat <- nhdist(macaques[["x"]], method = "riemann")

# unambiguous substring also available:
Dmat <- nhdist(macaques[["x"]], method = "rie")
```
