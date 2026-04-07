# Get the default prior specification for a given model family

A convenience wrapper around
[`set_priors()`](https://flo1met.github.io/bvarnet/reference/set_priors.md)
for inspecting defaults.

## Usage

``` r
get_default_priors(family)
```

## Arguments

- family:

  One of `"bernoulli"`, `"ordinal"`, `"gaussian"`.

## Value

A `bvarnet_priors` object.
