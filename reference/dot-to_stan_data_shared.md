# Build shared design matrices (family-agnostic)

Always includes the Intercept column in X. Does NOT pack priors or do
family-specific type casting. Used by both
[`to_stan_data()`](https://flo1met.github.io/bvarnet/reference/to_stan_data.md)
(joint path) and `.bvar_nodewise()` (mixed path).

## Usage

``` r
.to_stan_data_shared(
  data,
  id_col,
  time_col,
  y_cols,
  x_cols,
  center_x = FALSE,
  fe_interactions = NULL,
  re_interactions = NULL,
  re_cols = character(0),
  re_temporal = FALSE,
  K,
  na_action = "listwise",
  skip_lag = TRUE
)
```

## Value

A list with p, J, K, n_obs, n_fe, n_re, id, Y, X, B, Z, id_levels,
x_center_means, row_map, n_rows_data, design_spec, fe_interaction_terms,
fe_interaction_colnames.
