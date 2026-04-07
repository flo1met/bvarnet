# Group Stan beta names by lag×predictor interaction term (internal)

Group Stan beta names by lag×predictor interaction term (internal)

## Usage

``` r
get_lag_interaction_indices_by_term(sd)
```

## Arguments

- sd:

  The `standata` list.

## Value

Named list; each element has `full` (all params), `by_lag` (list of
per-lag-block param vectors), `ar` (AR-like interaction params where
lagged outcome == target outcome), and `cl` (CL-like interaction params
where lagged outcome != target outcome).
