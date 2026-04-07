# Format a bvarnet_prior for printing

Format a bvarnet_prior for printing

## Usage

``` r
# S3 method for class 'bvarnet_prior'
format(x, half = FALSE, ...)
```

## Arguments

- x:

  A `bvarnet_prior` object.

- half:

  Logical; if `TRUE` prepends "Half-" to indicate a half-prior (used for
  positive-constrained parameters like sd_u and sigma).

- ...:

  Ignored.

## Value

A character string.
