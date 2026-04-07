# Extract labelled parameter summaries from a fitted bvarnet model

Returns a single flat data frame with posterior summaries (mean, median,
5th/95th percentiles) and convergence diagnostics (Rhat, ESS) for all
model parameters.

## Usage

``` r
extract_param(object, bayes_factor = FALSE, null_value = 0)
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

## Value

A data frame with columns: `type`, `predictor`, `outcome`, `mean`,
`median`, `q5`, `q95`, `rhat`, `ess_bulk`, `ess_tail`, and optionally
`BF01`, `BF10`.
