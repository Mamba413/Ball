# meteorological data

A meteorological data include 46 records about air, soil, humidity, wind
and evaporation.

## Format

- `meteorology$air`: A data.frame containing 3 variables: maximum,
  minimum and average daily air temperature

  `meteorology$soil`: A data.frame containing 3 covariates: maximum,
  minimum and average daily soil temperature

  `meteorology$humidity`: A data.frame containing 3 covariates: maximum,
  minimum and average daily humidity temperature,

  `meteorology$wind`: a vector object record total wind, measured in
  miles per day

  `meteorology$evaporation`: a vector object record evaporation

## Details

This meteorological data containing 46 observations on five groups of
variables: air temperature, soil temperature, relative humidity, wind
speed as well as evaporation. Among them, maximum, minimum and average
value for air temperature, soil temperature, and relative humidity are
recorded. As regards to wind speed and evaporation, there are univariate
numerical variables. We desire to test the independence of these five
groups of variables.
