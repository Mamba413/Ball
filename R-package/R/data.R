#' @title meteorological data
#' @name meteorology
#' @docType data
#' @description 
#' A meteorological data include 46 records about air, soil, humidity, wind and evaporation.
#' 
#' @details 
#' This meteorological data containing 46 observations on five groups of variables: air temperature, soil temperature, 
#' relative humidity, wind speed as well as evaporation. 
#' Among them, maximum, minimum and average value for air temperature, soil temperature, and 
#' relative humidity are recorded. As regards to wind speed and evaporation, 
#' there are univariate numerical variables. We desire to test the independence of these five 
#' groups of variables. 
#' 
#' @format 
#' \itemize{
#' \code{meteorology$air}: A data.frame containing 3 variables: maximum, minimum and average daily air temperature
#' 
#' \code{meteorology$soil}: A data.frame containing 3 covariates: maximum, minimum and average daily soil temperature
#' 
#' \code{meteorology$humidity}: A data.frame containing 3 covariates: maximum, minimum and average daily humidity temperature, 
#' 
#' \code{meteorology$wind}: a vector object record total wind, measured in miles per day
#' 
#' \code{meteorology$evaporation}: a vector object record evaporation
#' }
#' 
NULL

#' @title Lung cancer genomic data
#' @name genlung
#' @docType data
#' @description Publicly available lung cancer genomic data from the Chemores Cohort Study, containing 
#' the expression levels of mRNA, miRNA, artificial noise variables as well as clinical variables and response.
#' 
#' @details 
#' Tissue samples were analysed from a cohort of 123 patients, who underwent complete surgical resection at the Institut Mutualiste
#' Montsouris (Paris, France) between 30 January 2002 and 26 June 2006. The studied outcome was the "Disease-Free Survival Time".
#' Patients were followed until the first relapse occurred or administrative censoring. In this genomic dataset, 
#' the expression levels of Agilent miRNA probes (\eqn{p=939}) were included from the \eqn{n=123} cohort samples. 
#' The miRNA data contains normalized expression levels. See below the paper by Lazar et al. (2013) and Array Express 
#' data repository for the complete description of the samples, tissue preparation, Agilent array technology, and data normalization. 
#' In addition to the genomic data, five clinical variables, also evaluated on the cohort samples, are included as 
#' continuous variable ('Age') and nominal variables ('Type','KRAS.status','EGFR.status','P53.status'). 
#' See Lazar et al. (2013) for more details. Moreover, we add 1056 standard gaussian variables 
#' which are independent with the censored response as noise covariates. This dataset represents a situation where the number of 
#' covariates dominates the number of complete observations or \eqn{p >> n} case.
#' 
#' @format 
#' \itemize{
#' \code{genlung$survival}: A data.frame containing \eqn{n=123} complete observations. 
#' The first column is disease-free survival time and the 
#' second column is censoring status. 
#' 
#' \code{genlung$covariate}: A data.frame containing \eqn{p=2000} covariates.
#' }
#' 
#' @references Lazar V. et al. (2013). Integrated molecular portrait of non-small cell lung cancers. BMC Medical Genomics 6:53-65.
NULL

#' @title Arctic lake sediment samples of different water depth
#' @name ArcticLake
#' @docType data
#' @description Sand, silt and clay compositions of 39 sediment samples of different water 
#' depth in an Arctic lake.
#' 
#' @details Sand, silt and clay compositions of 39 sediment samples at different water 
#' depth (in meters) in an Arctic lake. The additional feature is a concomitant variable or 
#' covariate, water depth, which may account for some of the variation in the compositions. 
#' In statistical terminology, we have a multivariate regression problem with sediment 
#' composition as predictors and water depth as a response. All row percentage sums to 100, 
#' except for rounding errors.
#' 
#' @format 
#' \itemize{
#' \code{ArcticLake$depth}: water depth (in meters). 
#' 
#' \code{ArcticLake$x}: compositions of three covariates: sand, silt, and clay.
#' }
#' 
#' @references Aitchison: The Statistical Analysis of Compositional Data, 1986, Data 5, pp5.
#' @source Aitchison: CODA microcomputer statistical package, 1986, the file name ARCTIC.DAT, here included under the GNU Public Library Licence Version 2 or newer.
#' @note Courtesy of J. Aitchison
NULL


#' @title Male and Female macaque data
#' @name macaques
#' @docType data
#' @description Male and female macaque skull data. 7 landmarks in 3 dimensions, 
#' 18 individuals (9 males, 9 females)
#' 
#' @details In an investigation into sex differences in the crania of a species of 
#' Macaca fascicularis (a type of monkey), random samples of 9 male and 9 female 
#' skulls were obtained by Paul O’Higgins (Hull-York Medical School) (Dryden and Mardia 1993). 
#' A subset of seven anatomical landmarks was located on each cranium and the three-dimensional (3D) 
#' coordinates of each point were recorded.
#' 
#' @format 
#' \itemize{
#' \code{macaques$x}: An array of dimension \eqn{7 \times 3 \times 18}
#' 
#' \code{macaques$group}: A factor indicating the sex ('m' for male and 'f' for female)
#' }
#' 
#' @note Dryden, I.L. and Mardia, K.V. (1998). Statistical Shape Analysis, Wiley, Chichester.
#' @references Dryden, I. L. and Mardia, K. V. (1993). Multivariate shape analysis. Sankhya Series A, 55, 460-480.
NULL


#' @title Simulated von Mises-Fisher Data
#' @name bdvmf
#' @docType data
#' 
#' @description Simulated random vectors following the von Mises-Fisher distribution 
#' with mean direction \eqn{\mu_{x}=(1, 0, 0)} and \eqn{\mu_{y}=(1, 1, 1)}, 
#' and concentration parameter is \eqn{\kappa = 3}.
#' 
#' @details In directional statistics, the von Mises–Fisher distribution 
#' (named after Ronald Fisher and Richard von Mises), is a probability distribution
#' on the \eqn{(p-1)}-dimensional sphere in \eqn{R^{p}}
#' 
#' The parameters \eqn{\mu}, and \eqn{\kappa}, are called the mean direction and concentration
#' parameter, respectively. The greater the value of \eqn{\kappa}, 
#' the higher the concentration of the distribution around the mean 
#' direction \eqn{\mu},. The distribution is unimodal for \eqn{\kappa}, 
#' and is uniform on the sphere for \eqn{\kappa=0}.
#' 
#' @format 
#' \itemize{
#' \code{bdvmf$x}: A \eqn{300 \times 3} numeric matrix containing simulated von Mises-Fisher data. 
#' 
#' \code{bdvmf$group}: A group index vector.
#' }
#' 
#' @references Embleton, N. I. Fisher, T. Lewis, B. J. J. (1993). Statistical analysis of spherical data (1st pbk. ed.). Cambridge: Cambridge University Press. pp. 115–116. ISBN 0-521-45699-1.
#' 
NULL


#' @title Distribution of BD when Null hypothesis is correct.
#' @name BDTestNullDistribution
#' @description Distribution of BD when Null hypothesis is correct.
#' @noRd
NULL

#' @title Distribution of BCor when Null hypothesis is correct.
#' @name BITestNullDistribution
#' @description Distribution of BCor when Null hypothesis is correct.
#' @noRd
NULL


# Various imports
#' @importFrom stats dist
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats as.formula
#' @importFrom stats rnorm
#' @importFrom utils data
#' @importFrom utils head
#' @importFrom utils memory.limit
#' @importFrom survival survfit
#' @importFrom stats lm
#' @importFrom gam gam
#' @importFrom gam s
#' @importFrom stats residuals
#' @importFrom mvtnorm rmvnorm
NULL



#' @title Generate example demostarted data in vignettes
#' @param n The desired sample size.
#' @param type This argument is used to select a specific dataset.
#' @description Three examples (bd-type1error, bd-unvariate, bd-multivariate) generate data suitable for demonstrating the 2-sample test for unvariate or multivariate variables.
#' The remaining examples (bcov-w, bcor-multivariate) are simple univariate and multivariate dependence structures used to demonstrate the tests of independece.
#' 
#' @details Choose a example type and set the size of the example.
#' \code{ball.example} will return an example list to demonstate
#' two-sample test or test of independece. Given data in examples(bd-type1error, bd-unvariate, bd-multivariate), we would like to test whether the two variables are identically distributed using \code{bd.test}.For data in examples(bcov-w, bcor-multivariate), we would like to test whether the two variables are statistically independent using \code{bcov.test}. Noted that they are both dependent.
#' 
#' @return A list contains two elements \code{x} and \code{y} is returned.
#' @noRd
ball.example <- function(n = 50, type) {
  set.seed(1)
  error <- runif(n, min = -0.3, max = 0.3)
  if(type %in% c('bd-type1error', 'bcov-type1error', 'bcor-type1error')) {
    x <- runif(n)
    y <- runif(n)
  }
  if(type %in% c("bcov-unvariate", "bcov-unvariate")) {
    x <- runif(n, 0, 4*pi)
    y <- cos(x) + error
  }
  if(type %in% c("bcor-multivariate", "bcov-multivariate")) {
    p <- 2
    x <- matrix(runif(n, -pi, pi), nrow = n, ncol = p)
    y <- (sin((x[,1])^2 + x[,2])) + error
  }
  if(type %in% c("bd-unvariate")) {
    x <- rnorm(n/2)
    y <- rnorm(n/2, mean = 1)
  }
  if(type %in% c("bd-multivariate")) {
    p <- 2
    x <- matrix(rnorm(n/2*p), nrow = n/2, ncol = p)
    y <- matrix(rnorm(n/2*p, mean = 3), nrow = n/2, ncol = p)
  }
  # if(type %in% c("bd-vmf")) {
  #   data("bdvmf", package = 'Ball', envir = environment(), verbose = FALSE)
  #   return(bdvmf)
  # }
  # if(type %in% c("bd-macaques")) {
  #   data("macaques", package = 'Ball', envir = environment(), verbose = FALSE)
  #   return(macaques)
  # }
  # if(type %in% c("bcov-arcticlake")) {
  #   data("ArcticLake", package = 'Ball', envir = environment(), verbose = FALSE)
  #   return(ArcticLake)
  # }
  list('y' = y, 'x' = x)
}


