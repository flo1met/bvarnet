# Validate variable names and return column indices (internal)

Validate variable names and return column indices (internal)

## Usage

``` r
.get_variable_indices(sd, variable)
```

## Arguments

- sd:

  The `standata` list from a `bvarnet` object.

- variable:

  Character vector of variable names to look up.

## Value

Integer vector of column positions in `Y`.
