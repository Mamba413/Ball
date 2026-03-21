# Changelog

## Ball 1.3.12

CRAN release: 2021-09-20

- Add a fast K-sample BD test procedure for category variables,
  especially the GWAS datasets.

## Ball 1.3.11

CRAN release: 2021-03-17

- Update citation information and reference.

## Ball 1.3.10

CRAN release: 2021-01-09

- Add variance (Chi-square) weight for K-sample Ball Divergence test.
- Refine distance computation

## Ball 1.3.9

CRAN release: 2019-12-17

- *bcorsis* can conveniently analyze category variables (e.g., GWAS
  datasets).
- Add permutation-free method for the Ball Divergence and Ball
  Covariance test (both only for constant weight now).
- Add variance (chi-square) weight for Ball Divergence test.
- Fix potential bugs in R version 4.0

## Ball 1.3.8

CRAN release: 2019-05-31

- Improve the reproducibility and efficiency of multi-threading and
  change the default number of threads as the number of CPUs available.
- Improve the computational efficiency of the KBD test.

## Ball 1.3.7

CRAN release: 2018-12-14

- Modify ambiguous arguments of bcov.test, bd.test, bcorsis
- Modify document

## Ball 1.3.6

CRAN release: 2018-12-03

- Formula interface for bd.test and bcov.test
- Optimize the package dependency

## Ball 1.3.5

CRAN release: 2018-11-02

- Faster implementation of mutual independence test
- Multi-thread support for the test of mutual independence
- Modify document

## Ball 1.3.0

CRAN release: 2018-09-30

- Add a KBD statistic designed for detecting the distribution
  distinction when a part of group distributions are identical. (setting
  *kbd.type = “maxsum”*)
- OPENMP based Multi-thread support for KBD
- Optimized OPENMP parallelism

## Ball 1.2.0

CRAN release: 2018-06-17

- Speed up feature screening.
- Speed up mutual independence test.
- Another K-sample test statistic (setting *kbd.type = “max”*) is
  implemented. It is good at detecting the distribution distinction when
  a part of group distributions are identical.
- Add *complete.info* component in the output list of *bcov.test,
  bd.test* and *bcorsis* such that user not need to re-run them when
  user want to obtain the test of result of different statistics.

## Ball 1.1.0

CRAN release: 2018-05-12

- Bug fix
- OPENMP based Multi-thread support for *bd.test* and *bcov.test*
- Speed up feature screening for survival data
- Implement Angular metric for compositional data

## Ball 1.0.0

CRAN release: 2018-01-03

- Initial CRAN version
