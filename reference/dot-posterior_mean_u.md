# Posterior mean of subject-level random effects `u`

Returns a 3D array `[node, subject, re]` containing the posterior mean
of each `u[node, subject, re]` element.

## Usage

``` r
.posterior_mean_u(object)
```

## Arguments

- object:

  A `bvarnet` object.

## Value

A 3D array with dimensions `[p, J, n_re]`.
