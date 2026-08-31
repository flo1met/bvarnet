# Compute Savege-Dickey Bayes factors

Computes Savage-Dickey density ratio Bayes factors for each (requested
set of) parameter in the model. By default, all applicable parameters
are tested and returned in a tidy data frame. The `type` argument
controls which parameter groups are included; the `variable` argument
can be used to filter to effects involving specific variables. The
`log_BF10` argument allows including the natural log of the Bayes factor
in the output, and `round` controls numeric rounding of the results.

## Usage

``` r
bf_table(
  object,
  type = "all",
  lag = 1L,
  null_value = 0,
  variable = NULL,
  log_BF10 = FALSE,
  round = 5L
)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- type:

  Character vector specifying which parameter groups to test. Use
  `"all"` (default) to include all applicable groups automatically.
  Available options:

  `"ar"`

  :   Autoregressive effects (self-loops). Per-cell BFs for the lag
      specified by `lag`, plus a joint BF.

  `"cl"`

  :   Cross-lagged effects. Same structure as `"ar"`.

  `"intercepts"`

  :   Intercept parameters. Skipped automatically for ordinal outcomes.

  `"fe"`

  :   Non-intercept fixed effects (covariates).

  `"lag_fe"`

  :   Joint BFs for lag \\\times\\ covariate interaction terms. Only
      available when the model was fitted with `fe_interactions`
      containing lag terms.

  `"temporal"`

  :   Joint BF for the entire temporal structure (all AR + CL parameters
      across all lags). When lag \\\times\\ covariate interactions are
      present, additional omnibus rows are included.

- lag:

  Integer; which lag block to use (default 1). Applies to `"ar"` and
  `"cl"` types.

- null_value:

  Numeric scalar; the null hypothesis value (default 0).

- variable:

  Character vector or `NULL` (default). One or more variable names. When
  set, only effects involving these variables are included.

- log_BF10:

  Logical; if `TRUE`, an additional `log_BF10` column (natural log of
  `BF10`) is appended to the output. Default is `FALSE`.

- round:

  Integer or `NULL`; number of decimal places to round numeric output
  columns. Default is `5`. Set to `NULL` to disable rounding.

## Value

A data frame with columns: `type`, `predictor`, `outcome`, `BF10` (and
optionally `log_BF10`).

## Which prior the Bayes factor divides by

The Savage-Dickey density ratio is only valid against the prior the
model was actually fitted with, so `bf_table()` evaluates the prior
density using `object$priors_effective` — resolved per outcome, since
Gaussian nodes scale their default priors by the outcome SD (see
[`set_priors`](https://flo1met.github.io/bvarnet/reference/set_priors.md)).
For a Gaussian outcome with default priors this makes `intercepts` and
`fe` Bayes factors differ from what the unscaled `object$priors` would
give, by roughly a factor of `sd(y)`. `ar`, `cl` and `temporal` tests
are unaffected: `phi` priors are never scaled.
