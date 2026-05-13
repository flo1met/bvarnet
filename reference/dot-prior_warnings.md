# Emit user-facing warnings about prior usage and return needed-prior names

Warns when user-set priors are not needed by the model, and messages
when the model uses default priors that the user did not explicitly set
(only if the user set at least one prior).

## Usage

``` r
.prior_warnings(priors, family_vec, n_re)
```

## Arguments

- priors:

  A `bvarnet_priors` object.

- family_vec:

  Named character vector of families per node.

- n_re:

  Integer. Number of random-effect columns from the built design.

## Value

Character vector of prior names the model actually uses.
