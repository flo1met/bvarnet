# Extract labelled parameter summaries from a fitted bvarnet model

Returns a single flat data frame with posterior summaries (mean, median,
and an equal-tailed credible interval) and convergence diagnostics
(Rhat, ESS) for all model parameters.

## Usage

``` r
extract_param(
  object,
  bayes_factor = FALSE,
  null_value = 0,
  type = NULL,
  ci_level = 0.95
)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- bayes_factor:

  Logical; if `TRUE`, append `BF01` and `BF10` columns computed via the
  Savage-Dickey density ratio for beta and phi parameters. Default
  `FALSE`.

- null_value:

  Numeric scalar; the null hypothesis value for Bayes factor computation
  (default 0). Only used when `bayes_factor = TRUE`.

- type:

  Character vector or `NULL` (default). If supplied, only rows matching
  the given type(s) are returned. Valid values are: `"Intercept"`,
  `"Fixed Effect"`, `"Autoregressive"`, `"Cross-lagged"`,
  `"Random Effect SD"`, `"Residual SD"`, `"Threshold"`.

- ci_level:

  Numeric scalar strictly between 0 and 1; the mass of the equal-tailed
  credible interval reported in `ci_lower` and `ci_upper`. Default
  `0.95`.

## Value

A data frame with columns: `type`, `predictor`, `outcome`, `mean`,
`median`, `ci_lower`, `ci_upper`, `rhat`, `ess_bulk`, `ess_tail`, and
optionally `BF01`, `BF10`.
