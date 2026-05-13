# Filter FE by_pred list to specific covariates (internal)

Filter FE by_pred list to specific covariates (internal)

## Usage

``` r
.filter_fe_by_covariate(by_pred, x_names)
```

## Arguments

- by_pred:

  Named list from
  [`get_beta_indices_by_predictor()`](https://flo1met.github.io/bvarnet/reference/get_beta_indices_by_predictor.md).

- x_names:

  Character vector of covariate names to keep.

## Value

Filtered list retaining only predictors in `x_names`.
