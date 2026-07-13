
# bvarnet

<!-- badges: start -->
[![R-CMD-check](https://github.com/flo1met/bvarnet/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/flo1met/bvarnet/actions/workflows/R-CMD-check.yml)
[![codecov](https://codecov.io/gh/flo1met/bvarnet/graph/badge.svg)](https://app.codecov.io/gh/flo1met/bvarnet)
[![pkgdown](https://github.com/flo1met/bvarnet/actions/workflows/pkgdown.yml/badge.svg)](https://flo1met.github.io/bvarnet/)
<!-- badges: end -->

## Bayesian Estimation of Multilevel Vector Autoregressive Networks using Stan

The `bvarnet` package allows user to estimate Bayesian multilevel Vector Auto Regressive (VAR) models for binary, ordinal and continuous outcome variables. Missing data is handled through listwise deletion and a skip-lag mechanism, which skips the estimation of the temporal structure when there is a gap between two timepoints.
Further, we provide functionality to conduct hypothesis test.

## Installation

`bvar()` fits its models with Stan, so it needs a compiled `Stan` executable for each of the three outcome families. Installing is a **two-step flow**, and works the same way whether or not you
have a C++ toolchain:

```r
# Step 1: install bvarnet from CRAN
install.packages("bvarnet")

# Step 2: set up the Stan models
bvarnet::bvarnet_setup_models()
```

`bvarnet_setup_models()` offers to either **download precompiled model binaries** for your platform (*recommended*) or **compile them locally** if you already have a working CmdStan installation. You only need to set this up once; re-run it only after updating `bvarnet` to a new version.

### No toolchain? (most users, recommended)

Just run `bvarnet::bvarnet_setup_models()` and choose the download option when prompted. This fetches a small, platform-specific set of precompiled Stan executables.

### Have a C++ toolchain? (advanced)

For this option you need to have [RTools](https://cran.r-project.org/bin/windows/Rtools/) (Windows) or [Xcode](https://developer.apple.com/xcode/) (Mac) installed. Further, you need to have cmdstanr installed and the C++ toolchain set up. If you don't have CmdStan yet, after installing RTools/Xcode install it using:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan(cores = 2)
```

If you run into any problems, see the
[Getting started with CmdStanR](https://mc-stan.org/cmdstanr/articles/cmdstanr.html) guide.

Then you can use `bvarnet::bvarnet_setup_models()` to compile the models locally. 

Alternatively, `install.packages("bvarnet", type = "source")` compiles the models at install time (requires CmdStan to already be set up). This works equivalent to the two-step flow above, just compiled into the package's install tree instead of a user cache directory.



### Development version

You can install the development version of `bvarnet` from [GitHub](https://github.com/flo1met/bvarnet). If you use the development version, you will have to compile the models yourself!

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("flo1met/bvarnet")
bvarnet::bvarnet_setup_models()
```

## Getting Started

The best place to start learning how to use this package to estimate Bayesian (multilevel) Vector Autoregression is the [Getting Started Vignette](https://flo1met.github.io/bvarnet/articles/bvarnet.html).
This vignette covers the basic model syntax, how to specify priors and how to extract the relevant parameters.

## Feature Requests and Contributions

- Predictions
- Cross-sectional Networks
- Correlated Random Effects
- Hierarchical Prior Distributions
- Performance Optimisation

## Roadmap

`bvarnet` is actively being developed. While the core functionality is stable, we have several features planned for future releases. For bug reports or feature request, please visit our [Issue Tracker](https://github.com/flo1met/bvarnet/issues).

