# Build a Stan data list from a long-format data frame

Internal function called by
[`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md).
Constructs the list passed to the Stan model: design matrices Y, X, B,
Z; dimensions; and prior hyperparameters derived from a `bvarnet_priors`
object.

## Usage

``` r
to_stan_data(
  data,
  family,
  id_col,
  time_col,
  y_cols,
  x_cols,
  center_x = F,
  fe_interactions = NULL,
  re_interactions = NULL,
  re_cols = character(0),
  re_temporal = FALSE,
  K,
  na_action = c("listwise"),
  skip_lag = TRUE,
  priors = set_priors()
)
```

## Arguments

- data:

  Data frame in long format.

- family:

  Character. `"bernoulli"`, `"ordinal"`, or `"gaussian"`.

- id_col:

  Character. Subject/group identifier column name.

- time_col:

  Character. Time column name.

- y_cols:

  Character vector. Outcome column names.

- x_cols:

  Character vector. Covariate column names.

- center_x:

  Logical. Grand-mean centre X. Default `FALSE`.

- fe_interactions:

  List or NULL. Fixed-effect interaction terms.

- re_interactions:

  List or NULL. Random-effect interaction terms.

- re_cols:

  Character vector. Random-slope columns.

- re_temporal:

  Logical. Random slopes on lag predictors. Default FALSE.

- K:

  Integer. AR order.

- na_action:

  Character. Only `"listwise"` currently supported.

- skip_lag:

  Logical. Zero-fill lags across irregular time gaps.

- priors:

  A `bvarnet_priors` object. Defaults to
  [`set_priors()`](https://flo1met.github.io/bvarnet/reference/set_priors.md).

## Value

A named list ready to pass to `CmdStanModel$sample()`.
