# Group Stan beta parameter names by predictor row (internal)

Returns a named list — one element per predictor row — each element
being a character vector of Stan names `beta[row, 1:p]` for that row.

## Usage

``` r
get_beta_indices_by_predictor(sd, type = c("fe", "intercepts"))
```

## Arguments

- sd:

  The `standata` list.

- type:

  One of `"fe"` or `"intercepts"`.

## Value

Named list of character vectors keyed by predictor name.
