# Resolve Stan parameter names for phi sub-matrices (internal)

Returns a character vector of Stan parameter names (e.g., `"phi[1,1]"`)
for autoregressive, cross-lagged, or all effects.

## Usage

``` r
get_phi_indices(sd, lag = 1L, effect = c("ar", "cl", "all"))
```

## Arguments

- sd:

  The `standata` list from a `bvarnet` object.

- lag:

  Integer; which lag block (1 to K). Default 1.

- effect:

  One of `"ar"`, `"cl"`, `"all"`.

## Value

Character vector of Stan parameter names.
