
# bvarnet

<!-- badges: start -->
[![R-CMD-check](https://github.com/flo1met/bvarnet/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/flo1met/bvarnet/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/flo1met/bvarnet/graph/badge.svg)](https://codecov.io/gh/flo1met/bvarnet)
[![pkgdown](https://github.com/flo1met/bvarnet/actions/workflows/pkgdown.yml/badge.svg)](https://flo1met.github.io/bvarnet/)
<!-- badges: end -->

## Bayesian Estimation of Multilevel Vector Autoregressive Networks using STAN

The `{bvarnet}` package allows user to estimate Bayesian multilevel Vector Auto Regressive (VAR) models for binary, ordinal and continuous outcome variables. Missing data is handled through listwise deletion and a skip-lag mechanism, which skips the estimation of the temporal structure when there is a gap between two timepoints.
Further, we provide functionality to conduct hypothesis test and perform predictions.

## Installation

You can install the development version of `{bvarnet}` from [GitHub](https://github.com/flo1met/bvarnet) with:

``` r
if(!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("flo1met/bvarnet")
```

## Getting Started

The best place to start learning how to use this package to estimate Bayesian (multilevel) Vector Autoregression is the [Getting Started Vignette]{https://flo1met.github.io/bvarnet/articles/bvarnet.html}.
This vignette covers the basic model syntax, how to specify priors and how to extract the relevant parameters.

## Feature Requests and Contributions

- Cross-sectional Networks
- Correlated Random Effects
- Hierarchical Prior Distributions
- Performance Optimisation

## Roadmap

`{bvarnet}` is actively being developed. While the core functionality is stable, we have several exciting features planned for future releases. For a granular look at our progress, known bugs, and feature requests, please visit our [Issue Tracker](https://github.com/username/packagename/issues).

