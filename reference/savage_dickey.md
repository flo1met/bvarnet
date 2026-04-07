# Compute the Savage-Dickey density ratio Bayes factor (internal)

This is the internal workhorse called by
[`bf_table()`](https://flo1met.github.io/bvarnet/reference/bf_table.md).
It returns a plain named list for a single test (one parameter or one
joint set of parameters).

## Usage

``` r
savage_dickey(
  object,
  params,
  null_value = 0,
  method = c("auto", "logspline", "mvn")
)
```

## Arguments

- object:

  A `bvarnet` object.

- params:

  Character vector of Stan parameter names.

- null_value:

  Numeric scalar or vector of null values (recycled if scalar).

- method:

  One of `"auto"`, `"logspline"`, `"mvn"`.

## Value

A named list with elements: `BF01`, `BF10`, `log_BF01`, `post_density`,
`prior_density`, `method`, `params`, `null_value`.
