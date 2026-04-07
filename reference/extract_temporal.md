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
  null_value = 0
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

## Value

A data frame with columns `type`, `predictor`, `outcome`, `mean`,
`median`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`, and optionally
`BF01`, `BF10`.
