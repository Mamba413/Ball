## Ball: Statistical Inference and Sure Independence Screening via Ball Statistics
<!-- [Chinese version](https://gitlab.com/mamba413/Ball/blob/develop_KS_V1.1/README_CN.md) -->

### Introduction
The fundamental problems for data mining and statistical analysis are:

- Whether distributions of two samples are distinct?

- Whether two random variables are dependent?

**Ball** package provides solution for these issues. Moreover, a variables (features) screening procedure also implemented to tackle ultra high dimensional data. The core functions in **Ball** package are: **bd.test**, **bcov.test**, and **bcorsis**.

These functions based on ball statistic have several advantages:

- It's applicable to univariate and multivariate data in banach space.

- There is no need for moment assumption, which means that outliers and heavy-tail data are no longer a problem.

- There is no need for complex adjustments for parameters.
 
Specially, for two-sample or K-sample test, **bd.test** has is proved to cope well for imbalanced data while **bcov.test** and **bcorsis** works well for complex responses and/or predictors, such as censored data and compositional data.     


### Installation        
#### CRAN version         
To install the Ball R package from CRAN, just run:        
```R
install.packages("Ball")
```

#### Github version       
To install the development version from GitHub, run:      
```R
library(devtools)
install_github("Mamba413/Ball", build_vignettes = TRUE)
```
*Windows* user will need to install [Rtools](https://cran.r-project.org/bin/windows/Rtools/) first.       


### Usage         
Take *iris* dataset as an example to illustrate how to use **bd.test** and **bcov.test** to 
deal with the fundamental problems mentioned above.

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
	Nonparametric 2-Samples Ball Divergence Test

data:  virginica and versicolor 
number of observations = 100, group sizes: 50 50
replicates = 99
bd = 0.32912, p-value = 0.01
alternative hypothesis: distributions of samples are distinct
```

The R output shows that *p*-value is under 0.05. Consequently, we can conclude that the Sepal.Length distribution of versicolor and virginica are distinct.

#### **bcov.test**        

```R
sepal <- iris[, c("Sepal.Width", "Sepal.Length")]
petal <- iris[, c("Petal.Width", "Petal.Length")]
bcov.test(sepal, petal)
```

In this example, **bcov.test** investigates whether width or length of petal is associated with width and length of sepal. If the dependence really exists, the *p*-value of the **bcov.test** will be under 0.05.

In this example, the result is:

```
	Ball Covariance test of independence

data:  sepal and petal
number of observations = 150
replicates = 99, Weighted Ball Covariance = FALSE
bcov = 0.0081472, p-value = 0.01
alternative hypothesis: random variables are dependent
```
Therefore, the relationship between width and length of sepal and petal is exists.

#### **bcorsis**                   

We simulate a dataset and demonstrate the usage of **bcorsis** function as follow.

```{r}
## simulate a ultra high dimensional dataset:
set.seed(1)
n <- 150
p <- 3000
x <- matrix(rnorm(n * p), nrow = n)
error <- rnorm(n)
y <- 3*x[, 1] + 5*(x[, 3])^2 + error

## BCor-SIS procedure:
res <- bcorsis(y = y, x = x)
head(res[["ix"]], n = 5)
```
      
In this example, the result is :

```
# [1]    3    1 1601   20  429
```
          
The **bcorsis** result shows that the first and the third variable are the two most 
important variables in 3000 explanatory variables which is consistent to simulation settings.

If you find any bugs, or if you experience any crashes, please report to us. If you have any questions just ask, we won't bite. 

<!-- Please cite our paper if you use Ball. -->

### License
GPL-3