# Parse lag×predictor interaction columns from X colnames (dev fallback)

Temporary internal helper for pre-implementation fit objects that lack
`fe_interaction_terms` metadata. Scheduled for removal before alpha.

## Usage

``` r
.parse_lag_interaction_cols(sd)
```

## Arguments

- sd:

  The `standata` list.

## Value

Same structure as
[`get_lag_interaction_indices_by_term()`](https://flo1met.github.io/bvarnet/reference/get_lag_interaction_indices_by_term.md).
