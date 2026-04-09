# Filter phi indices to effects FROM specific variables (internal)

Keeps only `phi[row, col]` entries where the lagged predictor
(`row_within`) is in `var_idx`.

## Usage

``` r
.filter_phi_by_variable(sd, lag, effect, var_idx)
```

## Arguments

- sd:

  The `standata` list.

- lag:

  Integer; which lag block.

- effect:

  One of `"ar"`, `"cl"`, `"all"`.

- var_idx:

  Integer vector of variable column positions.

## Value

Character vector of filtered Stan parameter names.
