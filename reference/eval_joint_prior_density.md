# Evaluate the joint prior density for independent priors (internal)

Evaluate the joint prior density for independent priors (internal)

## Usage

``` r
eval_joint_prior_density(prior_list, param_types, null_vec)
```

## Arguments

- prior_list:

  A named list of `bvarnet_prior` objects keyed by parameter type (e.g.,
  `"phi"`, `"beta"`).

- param_types:

  Character vector of prior types for each parameter.

- null_vec:

  Numeric vector of null values (same length as `param_types`).

## Value

Numeric scalar — the product of marginal densities.
