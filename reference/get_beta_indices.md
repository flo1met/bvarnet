# Resolve Stan parameter names for beta sub-groups (internal)

Resolve Stan parameter names for beta sub-groups (internal)

## Usage

``` r
get_beta_indices(sd, type = c("intercepts", "fe"))
```

## Arguments

- sd:

  The `standata` list.

- type:

  One of `"intercepts"` or `"fe"` (non-intercept fixed effects).

## Value

Character vector of Stan parameter names.
