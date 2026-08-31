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

## Automatic scaling of Gaussian defaults

Gaussian outcomes are modelled on their raw scale, so
[`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md) widens
the default `intercept`, `beta` and `sigma` scales by the outcome SD
before passing them to Stan. A default `beta ~ Normal(0, 1)` on data
with `sd(y) = 17.8` therefore becomes `Normal(0, 17.8)`. This keeps the
unit-scale defaults weakly informative whatever the units of `y`, but it
means the default you see here is not the prior that was used.

A prior you pass explicitly is taken at face value and never rescaled.
The priors actually used are recorded on the fitted object as
`priors_effective` (see
[`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md)), are
reported by [`print()`](https://rdrr.io/r/base/print.html), and are what
[`bf_table`](https://flo1met.github.io/bvarnet/reference/bf_table.md)
divides by.
