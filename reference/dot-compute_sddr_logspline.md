# Compute the SDDR using logspline (internal)

Compute the SDDR using logspline (internal)

## Usage

``` r
.compute_sddr_logspline(draws, prior, null = 0)
```

## Arguments

- draws:

  Numeric vector of posterior draws.

- prior:

  A `bvarnet_prior` object.

- null:

  Numeric scalar — null value.

## Value

Named list with elements `BF01`, `post_density`, `prior_density`.
