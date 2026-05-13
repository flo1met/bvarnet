# Computes Savage-Dickey density ratio Bayes factors for each parameter in the requested subset and returns a tidy data frame.

For `type = "fe"` and `"intercepts"`, the table contains three levels:
per-cell (logspline), per-predictor joint (MVN), and a global joint-all
(MVN). For `type = "ar"` and `"cl"`, the existing two-level structure
(per-cell + per-type joint) is unchanged.

## Usage

``` r
bf_table(object, type = "all", lag = 1L, null_value = 0, variable = NULL)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- type:

  Character vector or `"all"` (default). Which parameter groups to test.
  Options: `"ar"` (autoregressive), `"cl"` (cross-lagged),
  `"intercepts"`, `"fe"` (non-intercept fixed effects), `"lag_fe"` (lag
  × predictor interaction joint tests), `"temporal"` (joint test of all
  phi parameters across all lags, i.e. the entire temporal structure
  AR + CL, excluding covariates; additionally emits separate joint rows
  for AR-only and CL-only components; when lag × covariate interactions
  are present, additional rows are emitted for per-interaction-term and
  AR-only / CL-only interaction sub-tests, plus a full temporal +
  interactions omnibus). `"all"` auto-selects all applicable types
  (skips `"intercepts"` for ordinal models and `"lag_fe"` when no lag
  interactions exist). When `variable` is set, `"all"` also skips
  `"intercepts"` and `"fe"`. Per-cell `"ar"` and `"cl"` rows respect the
  `lag` argument; `"temporal"` always covers all lags via joint tests.

- lag:

  Integer; which lag block to use (default 1). Applies to `"ar"` and
  `"cl"` types.

- null_value:

  Numeric scalar; the null hypothesis value (default 0).

- variable:

  Character vector or `NULL` (default). One or more variable names —
  either network variables (from `colnames(standata$Y)`) or covariates
  (from `x_cols`, excluding `"Intercept"` and interaction columns). When
  set, only effects involving these variables are included: network
  variables filter phi rows (effects **from** the variable as lagged
  predictor); covariate names filter fixed-effect rows and lag ×
  covariate interaction rows (effects **of** that covariate). Both types
  can be combined in a single call. Cannot be combined with
  `type = "intercepts"`.

## Value

A data frame with columns: `type`, `predictor`, `outcome`, `BF01`,
`BF10`, `log_BF01`, `post_density`, `prior_density`, `method`.

## Details

`type = "lag_fe"` emits only grouped joint rows for lag × predictor
interaction terms: per-lag-block and full-term omnibus. Per-cell rows
for these parameters are already included when `type = "fe"` is
requested.

When `variable` is non-NULL, only effects involving the named
variable(s) are included. Network variables (from `Y`) filter phi rows
and lag × covariate interaction rows by the lagged predictor. Covariate
names (from `X`) filter fixed-effect rows and interaction terms. Both
types can be combined. `variable` is combinable with `type` and `lag`;
when `type = "all"` and `variable` is set, the auto-selected types
depend on what was requested (network variables → `"ar"`, `"cl"`,
`"temporal"`; covariates → `"fe"`, `"lag_fe"`, `"temporal"` interaction
rows).
