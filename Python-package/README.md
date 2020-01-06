<img src=https://github.com/Mamba413/git_picture/blob/master/scrcss.jpg width=135/> Ball Statistics
===========

[![Travis Build Status](https://travis-ci.org/Mamba413/Ball.svg?branch=master)](https://travis-ci.org/Mamba413/Ball)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Mamba413/Ball?branch=master&svg=true)](https://ci.appveyor.com/project/Mamba413/Ball)
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

Installation
------------
### Pypi version         
To install the Ball Python package from Pypi, just run:        
```
pip install Ball
```

Overview: **Ball** package
----------
Three most importance functions in **Ball**:		

|                       |         **bd.test**         |     **bcov.test**      |    **bcorsis**     |
| --------------------- | :-------------------------: | :--------------------: | :----------------: |
| Feature               |       Hypothesis test       |    Hypothesis test     | Feature screening  |
| Type                  | Test of equal distributions |  Test of independence  |    SIS and ISIS    |
| Optional weight       |     :heavy_check_mark:      |   :heavy_check_mark:   | :heavy_check_mark: |
| Parallel programming  |     :heavy_check_mark:      |   :heavy_check_mark:   | :heavy_check_mark: |
| p-value               |     :heavy_check_mark:      |   :heavy_check_mark:   |        :x:         |
| Limit distribution    |    Two-sample test only     | Independence test only |        :x:         |
| Censored data         |             :x:             |          :x:           |    Comming soon    |
| Interaction screening |             :x:             |          :x:           | :heavy_check_mark: |
| GWAS optimization     |             :x:             |          :x:           |    Comming soon    |

- *SIS: Sure Independence Screening*
- *ISIS: Iterative Sure Independence Screening (SIS)*
- *GWAS: Genome-Wide Association Study*

Reference
----------
- Pan, Wenliang; Tian, Yuan; Wang, Xueqin; Zhang, Heping. [Ball Divergence: Nonparametric two sample test](https://projecteuclid.org/euclid.aos/1525313077). Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. [https://projecteuclid.org/euclid.aos/1525313077](https://projecteuclid.org/euclid.aos/1525313077)
- Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) [A Generic Sure Independence Screening Procedure](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1462709#.WupWaoiFM2x), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
- Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019) [Ball Covariance: A Generic Measure of Dependence in Banach Space](https://doi.org/10.1080/01621459.2018.1543600), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1543600
- Jin, Z., Wenliang P., Wei Z., and Xueqin W. (2018). Ball: An R package for detecting distribution difference and association in metric spaces. arXiv preprint arXiv:1811.03750. URL http://arxiv.org/abs/1811.03750.

Bug report
----------
If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. Open an [issue](https://github.com/Mamba413/Ball/issues) or send an email to Jin Zhu at zhuj37@mail2.sysu.edu.cn
