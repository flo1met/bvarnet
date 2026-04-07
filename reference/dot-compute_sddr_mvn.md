# Compute the SDDR using a multivariate normal approximation (internal)

Compute the SDDR using a multivariate normal approximation (internal)

## Usage

``` r
.compute_sddr_mvn(draws_mat, prior_list, param_types, null_vec)
```

## Arguments

- draws_mat:

  Numeric matrix — S rows x d columns of posterior draws.

- prior_list:

  Named list of `bvarnet_prior` objects.

- param_types:

  Character vector of prior types for each column.

- null_vec:

  Numeric vector of null values.

## Value

Named list with elements `BF01`, `post_density`, `prior_density`.
