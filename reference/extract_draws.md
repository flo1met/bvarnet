# Extract raw posterior draws for a single parameter block

Returns an `(iterations * chains)` by `params` matrix with Stan-indexed
column names (e.g. `"beta[1,1]"`, `"phi[2,3]"`).

## Usage

``` r
extract_draws(object, parameter = c("beta", "phi", "sd_u", "sigma", "kappa"))
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- parameter:

  Character. One of `"beta"`, `"phi"`, `"sd_u"`, `"sigma"`, or
  `"kappa"`.

## Value

A numeric matrix with one row per posterior draw and one column per Stan
parameter element.
