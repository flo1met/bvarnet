# Filter lag-interaction term groups to a specific set of covariates (internal)

Filter lag-interaction term groups to a specific set of covariates
(internal)

## Usage

``` r
.filter_lag_interaction_by_covariate(sd, x_names)
```

## Arguments

- sd:

  The `standata` list.

- x_names:

  Character vector of covariate names to keep.

## Value

Filtered version of
[`get_lag_interaction_indices_by_term()`](https://flo1met.github.io/bvarnet/reference/get_lag_interaction_indices_by_term.md),
retaining only terms whose non-lag suffix matches `x_names`.
