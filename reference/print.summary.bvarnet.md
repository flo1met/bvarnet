# Print a bvarnet summary

Pretty-prints the output of
[`summary.bvarnet`](https://flo1met.github.io/bvarnet/reference/summary.bvarnet.md),
grouping parameters by type and displaying convergence information. Each
group is truncated to `max_rows` rows; use
[`extract_param()`](https://flo1met.github.io/bvarnet/reference/extract_param.md)
or dedicated extractors to see full output.

## Usage

``` r
# S3 method for class 'summary.bvarnet'
print(x, digits = 3, max_rows = 10, ...)
```

## Arguments

- x:

  A `summary.bvarnet` object.

- digits:

  Number of decimal digits for numeric columns. Default 3.

- max_rows:

  Maximum number of rows to print per parameter group. Default 10.

- ...:

  Ignored.

## Value

`x` invisibly.
