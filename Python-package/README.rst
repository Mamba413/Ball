 Ball Statistics
================


Introduction
~~~~~~~~~~~~

The fundamental problems for data mining and statistical analysis are:

-  Whether distributions of two samples are distinct?

-  Whether two random variables are dependent?

**Ball** package provides solutions for these issues. Moreover, a
variable screening (or feature screening) procedure is also implemented
to tackle ultra high dimensional data. The core functions in **Ball**
package are **bd.test**, **bcov.test**, and **bcorsis**.

These functions based on ball statistic have several advantages:

-  It's applicable to univariate and multivariate data in Banach space.

-  There is no need for moment assumption, which means that outliers and
   heavy-tail data are no longer a problem.

-  They perform well in many setting without complex adjustments for
   parameters.

Particularly, for two-sample or K-sample problem, **bd.test** has been
proved to cope well for imbalanced data, and **bcov.test** and
**bcorsis** work well for detecting the relationship between complex
responses and/or predictors, such as shape, compositional as well as
censored data.

