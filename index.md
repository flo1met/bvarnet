# bvarnet

## Bayesian Estimation of Multilevel Vector Autoregressive Networks using STAN

The ‘bvarnet’ package allowes user to estimate bayesian multilevel
Vector Auto Regressive (VAR) models for binary, ordinal and continuous
outcome variables. Missing data is handled through listwise deletion and
a skip-lag mechanism, which skips the estimation of the temporal
structure when there is a gap between two timepoints. Further, we
provide functionality to perform hypothesis test and perform
predictions.

## Installation

You can install the development version of bvarnet from
[GitHub](https://github.com/) with:

``` r
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("flo1met/bvarnet")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bvarnet)
## basic example code
```
