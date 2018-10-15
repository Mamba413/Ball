# Ball 1.3.0
* Add a KBD statistic designed for detecting the distribution distinction when a part of group distributions are identical. (setting *kbd.type = "maxsum"*)
* OPENMP based Multi-thread support for KBD
* Optimized OPENMP parallelism

# Ball 1.2.0
* Speed up feature screening.
* Speed up mutual independence test.
* Another K-sample test statistic (setting *kbd.type = "max"*) is implemented. It is good at detecting the distribution distinction when a part of group distributions are identical.
* Add *complete.info* component in the output list of *bcov.test, bd.test* and *bcorsis* such that user not need to re-run them when user want to obtain the test of result of different statistics.

# Ball 1.1.0
* Bug fix
* OPENMP based Multi-thread support for *bd.test* and *bcov.test*
* Speed up feature screening for survival data
* Implement Angular metric for compostional data 

# Ball 1.0.0

* Initial CRAN version



