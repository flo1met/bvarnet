# Summary method for bvarnet objects

Returns a labelled posterior summary table grouped by parameter type,
with convergence diagnostics and optional Bayes factors. Wraps
[`extract_param`](https://flo1met.github.io/bvarnet/reference/extract_param.md).

## Usage

``` r
# S3 method for class 'bvarnet'
summary(object, bayes_factor = FALSE, null_value = 0, ..., ci_level = 0.95)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- bayes_factor:

  Logical; if `TRUE`, append Savage-Dickey BF columns. Default `FALSE`.

- null_value:

  Numeric scalar; null hypothesis value for BF computation. Default 0.

- ...:

  Ignored.

- ci_level:

  Numeric scalar strictly between 0 and 1; the mass of the equal-tailed
  credible interval reported in `ci_lower` and `ci_upper`. Default
  `0.95`. Must be supplied by name.

## Value

An object of class `"summary.bvarnet"` (a list) with elements:

- table:

  Data frame from
  [`extract_param()`](https://flo1met.github.io/bvarnet/reference/extract_param.md).

- family:

  Model family.

- p:

  Number of outcome variables.

- K:

  AR order.

- n:

  Number of observations.

- rhat_max:

  Maximum Rhat across all parameters.

- n_divergences:

  Total divergent transitions.

- ci_level:

  Credible-interval mass used for the table.
