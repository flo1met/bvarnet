# Build a prior specification object for `bvar()`

Returns a `bvarnet_priors` object containing a `bvarnet_prior` for every
model parameter type. Any argument left as `NULL` uses the package
default. Available prior distributions are:

- normal(loc, scale)

- student_t(loc, scale, df)

- cauchy(loc, scale) For standart deviations and random effects, the
  prior is automatically converted to a half-prior (truncated at `loc`)
  in the Stan code, so the printed format reflects this.

## Usage

``` r
set_priors(
  intercept = NULL,
  beta = NULL,
  phi = NULL,
  sd_u = NULL,
  kappa = NULL,
  sigma = NULL
)
```

## Arguments

- intercept:

  Prior for the intercept. Only applies to gaussian and bernoulli
  models; for ordinal models the intercept is absorbed into the kappa
  (threshold parameter).

- beta:

  Prior for fixed-effect regression coefficients (slopes).

- phi:

  Prior for lag coefficients.

- sd_u:

  Prior for random-effect standard deviations (half-prior).

- kappa:

  Prior for ordinal cut-points (ordinal models only).

- sigma:

  Prior for residual standard deviation (gaussian models only;
  half-prior).

## Value

A `bvarnet_priors` S3 object.
