# Filter lag-interaction indices to effects FROM specific variables (internal)

Filter lag-interaction indices to effects FROM specific variables
(internal)

## Usage

``` r
.filter_lag_interaction_by_variable(sd, var_idx)
```

## Arguments

- sd:

  The `standata` list.

- var_idx:

  Integer vector of variable column positions.

## Value

Same structure as
[`get_lag_interaction_indices_by_term()`](https://flo1met.github.io/bvarnet/reference/get_lag_interaction_indices_by_term.md)
but with only parameters where the lagged variable is in `var_idx`.
