# Validate variable names and classify into network/covariate (internal)

Validate variable names and classify into network/covariate (internal)

## Usage

``` r
.classify_variable(sd, variable)
```

## Arguments

- sd:

  The `standata` list from a `bvarnet` object.

- variable:

  Character vector of variable names to look up.

## Value

A list with components `y_idx` (integer indices into Y columns, or
`NULL`), `x_names` (character vector of matched covariate names, or
`NULL`), and `has_y`/`has_x` logicals.
