<img src=https://github.com/Mamba413/git_picture/blob/master/scrcss.jpg width=135/> Ball Statistics
===========

[![Travis Build Status](https://travis-ci.org/Mamba413/Ball.svg?branch=master)](https://travis-ci.org/Mamba413/Ball)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Mamba413/Ball?branch=master&svg=true)](https://ci.appveyor.com/project/Mamba413/Ball)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/Ball)](https://CRAN.R-project.org/package=Ball)
[![PyPI version](https://badge.fury.io/py/Ball.svg)](https://pypi.python.org/pypi/Ball/)

Introdution
----------
The fundamental problems for data mining, statistical analysis, and machine learning are:
- whether several distributions are different?
- whether random variables are dependent?
- how to pick out useful variables/features from a high-dimensional data?

These issues can be tackled by Ball statistics, which enjoy following admirable advantages:
- available for most of datasets (e.g., traditional tabular data, brain shape, functional connectome, wind direction and so on)
- insensitive to outliers, distribution-free and model-free;
- theoretically guaranteed and computationally efficient.

Softwares
----------
### R package
Install the **Ball** package from CRAN:        
```R
install.packages("Ball")
```
Compared with selective R packages available for datasets in metric spaces:

|                                   | [fastmit](https://cran.r-project.org/web/packages/fastmit) | [energy](https://cran.r-project.org/web/packages/energy) | [HHG](https://cran.r-project.org/web/packages/HHG) | [Ball](https://cran.r-project.org/web/packages/Ball) |
| :-------------------------------- | :----------------------------------------------------------: | :--------------------------------------------------------: | :--------------------------------------------------: | :----------------------------------------------------: |
| Test of equal distributions       | :x:                                                        | :heavy_check_mark:                                       | :heavy_check_mark:                                 | :heavy_check_mark:                                   |
| Test of independence              | :heavy_check_mark:                                         | :heavy_check_mark:                                       | :heavy_check_mark:                                 | :heavy_check_mark:                                   |
| Test of joint independence        | :x:                                                        | :x:                                                      | :x:                                                | :heavy_check_mark:                                   |
| Feature screening/SIS[^sis] | :x:                                                        | :x:                                                      | :x:                                                | :heavy_check_mark:                                   |
| Iterative Feature screening/ISIS[^isis]                     | :x:                                                        | :x:                                                      | :x:                                                | :heavy_check_mark:                                   |
| Datasets in metric spaces         | :heavy_check_mark:                                         | SNT[^snt]                                    | :heavy_check_mark:                                 | :heavy_check_mark:                                   |
| Robustness                        | :heavy_check_mark:                                         | :x:                                                      | :heavy_check_mark:                                 | :heavy_check_mark:                                   |
| Parallel programming              | :x:                                                        | :x:                                                      | :heavy_check_mark:                                 | :heavy_check_mark:                                   |
| Computational efficiency          | :running::running::running:                                | :running::running::running:                              | :running::running:                                 | :running::running::walking:                          |
[^sis]: Sure Independence Screening (SIS)
[^isis]: Iterative Sure Independence Screening (SIS)
[^snt]: strong negative type (SNT)

See the following documents for more details about the **[Ball](https://cran.r-project.org/web/packages/Ball)** package:
- [github page](https://github.com/Mamba413/Ball/tree/master/R-package) (short)
- [vignette](https://cran.r-project.org/web/packages/Ball/vignettes/Ball.html) (moderate)
- [JSS paper](https://arxiv.org/abs/1811.03750) (detailed)

### Python package
Install the **Ball** package from PyPI:        
```shell
pip install Ball
```

References
----------
- Pan, Wenliang; Tian, Yuan; Wang, Xueqin; Zhang, Heping. [Ball Divergence: Nonparametric two sample test](https://projecteuclid.org/euclid.aos/1525313077). Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. [https://projecteuclid.org/euclid.aos/1525313077](https://projecteuclid.org/euclid.aos/1525313077)
- Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) [A Generic Sure Independence Screening Procedure](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1462709#.WupWaoiFM2x), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
- Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019) [Ball Covariance: A Generic Measure of Dependence in Banach Space](https://doi.org/10.1080/01621459.2018.1543600), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1543600
- Jin, Z., Wenliang P., Wei Z., and Xueqin W. (2018). Ball: An R package for detecting distribution difference and association in metric spaces. arXiv preprint arXiv:1811.03750. URL http://arxiv.org/abs/1811.03750.

Bug report
----------
Open an [issue](https://github.com/Mamba413/Ball/issues) or send an email to Jin Zhu at zhuj37@mail2.sysu.edu.cn

Footnotes
----------