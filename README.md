<img src=https://github.com/Mamba413/git_picture/blob/master/scrcss.jpg width=135/> Ball Statistics
===========

[![Travis Build Status](https://travis-ci.org/Mamba413/Ball.svg?branch=master)](https://travis-ci.org/Mamba413/Ball)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Mamba413/Ball?branch=master&svg=true)](https://ci.appveyor.com/project/Mamba413/Ball)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/Ball)](https://CRAN.R-project.org/package=Ball)
[![PyPI version](https://badge.fury.io/py/Ball.svg)](https://pypi.python.org/pypi/Ball/)

Introdution
----------
The fundamental problems for data mining and statistical analysis are:

- Whether distributions of two samples are distinct?

- Whether two random variables are dependent?

Ball Statistics provide solutions for these issues. Moreover, a variable screening (or feature screening) procedure is also implemented to tackle ultra high dimensional data. 

Ball statistic have several advantages:

- It's applicable to univariate and multivariate data in Banach space.

- There is no need for moment assumption, which means that outliers and heavy-tail data are no longer a problem.

- They perform well in many setting without complex adjustments for parameters.
 
Particularly, for two-sample or K-sample problem, Ball Divergence has been proved to cope well for imbalanced data, and Ball Covariance and Ball Correlation work well for detecting the relationship between complex responses and/or predictors, such as shape, compositional as well as censored data.     

License
----------
GPL-3

Reference
----------
- Pan, Wenliang; Tian, Yuan; Wang, Xueqin; Zhang, Heping. [Ball Divergence: Nonparametric two sample test](https://projecteuclid.org/euclid.aos/1525313077). Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. [https://projecteuclid.org/euclid.aos/1525313077](https://projecteuclid.org/euclid.aos/1525313077)
- Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) [A Generic Sure Independence Screening Procedure](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1462709#.WupWaoiFM2x), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
- Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019) [Ball Covariance: A Generic Measure of Dependence in Banach Space](https://doi.org/10.1080/01621459.2018.1543600), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1543600
- Jin, Z., Wenliang P., Wei Z., and Xueqin W. (2018). Ball: An R package for detecting distribution difference and association in metric spaces. arXiv preprint arXiv:1811.03750. URL http://arxiv.org/abs/1811.03750.
