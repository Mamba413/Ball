Ball Statistics in R <img src='https://raw.githubusercontent.com/Mamba413/git_picture/master/ball_logo.png' align="right" height="120" />
===========

<!-- [![Travis Build Status](https://travis-ci.org/Mamba413/Ball.svg?branch=master)](https://travis-ci.org/Mamba413/Ball) -->
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Mamba413/Ball?branch=master&svg=true)](https://ci.appveyor.com/project/Mamba413/Ball)
[![CRAN Status Badge](http://www.r-pkg.org/badges/version/Ball)](https://CRAN.R-project.org/package=Ball)

Introduction
----------
The fundamental problems for data mining, statistical analysis, and machine learning are:
- whether several distributions are different?
- whether random variables are dependent?
- how to pick out useful variables/features from a high-dimensional data?

These issues can be tackled by using **bd.test**, **bcov.test**, and **bcorsis** functions in the **Ball** package, respectively. They enjoy following admirable advantages:
- available for most of datasets (e.g., traditional tabular data, brain shape, functional connectome, wind direction and so on)
- insensitive to outliers, distribution-free and model-free;
- theoretically guaranteed and computationally efficient.

Installation
----------
### CRAN version         
To install the Ball R package from CRAN, just run:        
```R
install.packages("Ball")
```

### Github version       
To install the development version from GitHub, run:      
```R
library(devtools)
install_github("Mamba413/Ball/R-package", build_vignettes = TRUE)
```
*Windows* user will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.       

Overview: **Ball** package
----------
Three most importance functions in **Ball**:		

|                          |         **bd.test**         |        **bcov.test**         |        **bcorsis**        |
| ------------------------ | :-------------------------: | :--------------------------: | :-----------------------: |
| Feature                  |       Hypothesis test       |       Hypothesis test        |     Feature screening     |
| Type                     | Test of equal distributions | Test of (joint) independence | SIS and ISIS |
| Optional weight          |     :heavy_check_mark:      |      :heavy_check_mark:      |    :heavy_check_mark:     |
| Parallel programming     |     :heavy_check_mark:      |      :heavy_check_mark:      |    :heavy_check_mark:     |
| p-value                  |     :heavy_check_mark:      |      :heavy_check_mark:      |            :x:            |
| Limit distribution       |    Two-sample test only     |    Independence test only    |            :x:            |
| Censored data            |             :x:             |             :x:              |    :heavy_check_mark:     |
| Interaction screening    |             :x:             |             :x:              |    :heavy_check_mark:     |
| GWAS optimization |             :x:             |             :x:              |    :heavy_check_mark:     |

- *SIS: Sure Independence Screening*
- *ISIS: Iterative Sure Independence Screening (SIS)*
- *GWAS: Genome-Wide Association Study*

Quick examples
----------
Take *iris* dataset as an example to illustrate how to use **bd.test** and **bcov.test** to deal with the fundamental problems mentioned above.

#### **bd.test**              
```R
virginica <- iris[iris$Species == "virginica", "Sepal.Length"]
versicolor <- iris[iris$Species == "versicolor", "Sepal.Length"]
bd.test(virginica, versicolor)
```

In this example, **bd.test** examines the assumption that Sepal.Length distributions of versicolor and virginica are equal.

If the assumption invalid, the *p*-value of the **bd.test**  will be under 0.05.

In this example, the result is:

```
	2-sample Ball Divergence Test (Permutation)

data:  virginica and versicolor 
number of observations = 100, group sizes: 50 50
replicates = 99, weight: constant
bd.constant = 0.11171, p-value = 0.01
alternative hypothesis: distributions of samples are distinct
```

The R output shows that *p*-value is under 0.05. Consequently, we can conclude that the Sepal.Length distribution of versicolor and virginica are distinct.

#### **bcov.test**        

```R
sepal <- iris[, c("Sepal.Width", "Sepal.Length")]
petal <- iris[, c("Petal.Width", "Petal.Length")]
bcov.test(sepal, petal)
```

In this example, **bcov.test** investigates whether width or length of petal is associated with width and length of sepal. If the dependency really exists, the *p*-value of the **bcov.test** will be under 0.05. In this example, the result is show to be:

```
	Ball Covariance test of independence (Permutation)

data:  sepal and petal
number of observations = 150
replicates = 99, weight: constant
bcov.constant = 0.0081472, p-value = 0.01
alternative hypothesis: random variables are dependent
```
Therefore, the relationship between width and length of sepal and petal exists.

#### **bcorsis**                   
We generate a dataset and demonstrate the usage of **bcorsis** function as follow.

```{r}
## simulate a ultra high dimensional dataset:
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
error <- rnorm(n)
y <- 3 * x[, 1] + 5 * (x[, 3])^2 + error

## BCor-SIS procedure:
res <- bcorsis(y = y, x = x)
head(res[["ix"]], n = 5)
```
In this example, the result is:
```{r}
# [1]    3    1 1601   20  429
```
The **bcorsis** result shows that the first and the third variable are the two most 
important variables in 3000 explanatory variables which is consistent to the simulation settings.

<!-- If you use **Ball** or reference our blog post in a presentation or publication, we would appreciate citations of our package. Here is the corresponding Bibtex entry:
```
@article{zhu2018ball,
  title={Ball: An R package for detecting distribution difference and association in metric spaces},
  author={Zhu, Jin and Pan, Wenliang and Zheng, Wei and Wang, Xueqin},
  journal={arXiv preprint arXiv:1811.03750},
  year={2018}
}
``` -->

Citation
----------
If you use Ball or reference our vignettes in a presentation or publication, we would appreciate citations of our package.
> Zhu J, Pan W, Zheng W, Wang X (2021). “Ball: An R Package for Detecting Distribution Difference and Association in Metric Spaces.” Journal of Statistical Software, 97(6), 1–31. doi: 10.18637/jss.v097.i06.

Here is the corresponding Bibtex entry
```
@Article{,
  title = {{Ball}: An {R} Package for Detecting Distribution Difference and Association in Metric Spaces},
  author = {Jin Zhu and Wenliang Pan and Wei Zheng and Xueqin Wang},
  journal = {Journal of Statistical Software},
  year = {2021},
  volume = {97},
  number = {6},
  pages = {1--31},
  doi = {10.18637/jss.v097.i06},
}
```

Reference
----------
- Pan, Wenliang; Tian, Yuan; Wang, Xueqin; Zhang, Heping. [Ball Divergence: Nonparametric two sample test](https://projecteuclid.org/euclid.aos/1525313077). Ann. Statist. 46 (2018), no. 3, 1109--1137. doi:10.1214/17-AOS1579. [https://projecteuclid.org/euclid.aos/1525313077](https://projecteuclid.org/euclid.aos/1525313077)
- Wenliang Pan, Xueqin Wang, Weinan Xiao & Hongtu Zhu (2018) [A Generic Sure Independence Screening Procedure](https://amstat.tandfonline.com/doi/full/10.1080/01621459.2018.1462709#.WupWaoiFM2x), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1462709
- Wenliang Pan, Xueqin Wang, Heping Zhang, Hongtu Zhu & Jin Zhu (2019) [Ball Covariance: A Generic Measure of Dependence in Banach Space](https://doi.org/10.1080/01621459.2018.1543600), Journal of the American Statistical Association, DOI: 10.1080/01621459.2018.1543600
- Jin Zhu, Wenliang Pan, Wei Zheng, and Xueqin Wang (2021). [Ball: An R Package for Detecting Distribution Difference and Association in Metric Spaces](https://www.jstatsoft.org/article/view/v097i06), Journal of Statistical Software, Vol.97(6), doi: 10.18637/jss.v097.i06.

Bug report
----------
If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. Open an [issue](https://github.com/Mamba413/Ball/issues) or send an email to Jin Zhu at zhuj37@mail2.sysu.edu.cn