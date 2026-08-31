# Extract temporal (VAR lag) effects

Returns a data frame of autoregressive and/or cross-lagged parameter
summaries with convergence diagnostics, filtered by lag and effect type.

## Usage

``` r
extract_temporal(
  object,
  lag = NULL,
  effect = c("all", "ar", "cl"),
  bayes_factor = FALSE,
  null_value = 0,
  ci_level = 0.95
)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- lag:

  Integer or `NULL`. If specified, only effects from this lag are
  returned. Default `NULL` (all lags).

- effect:

  Character. One of `"all"` (default), `"ar"` (autoregressive only), or
  `"cl"` (cross-lagged only).

- bayes_factor:

  Logical; if `TRUE`, append BF columns. Default `FALSE`.

- null_value:

  Numeric; null hypothesis for BF. Default 0.

- ci_level:

  Numeric scalar strictly between 0 and 1; the mass of the equal-tailed
  credible interval reported in `ci_lower` and `ci_upper`. Default
  `0.95`.

## Value

A data frame with columns `type`, `predictor`, `outcome`, `mean`,
`median`, `ci_lower`, `ci_upper`, `rhat`, `ess_bulk`, `ess_tail`, and
optionally `BF01`, `BF10`.
