# Ensure all prior slots required by the model family exist

If slots are missing (e.g. from a filtered
[`get_default_priors()`](https://flo1met.github.io/bvarnet/reference/get_default_priors.md)
object), they are filled with package defaults and a warning is issued.
Called by
[`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md) before
passing priors to
[`to_stan_data()`](https://flo1met.github.io/bvarnet/reference/to_stan_data.md).

## Usage

``` r
.ensure_prior_slots(priors, family_vec)
```

## Arguments

- priors:

  A `bvarnet_priors` object (possibly incomplete).

- family_vec:

  Named character vector of families per node.

## Value

A complete `bvarnet_priors` object.
