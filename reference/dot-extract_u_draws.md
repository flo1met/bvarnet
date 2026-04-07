# Extract all posterior draws of subject-level random effects `u`

Returns a 4D array with dimensions `[draw, node, subject, re]` matching
the Stan declaration `array[p] matrix[J, n_re] u`.

## Usage

``` r
.extract_u_draws(object)
```

## Arguments

- object:

  A `bvarnet` object.

## Value

A 4D array with dimensions `[S, p, J, n_re]` where
`S = n_iter * n_chains`.

## Details

CmdStan flattens `array[p] matrix[J, n_re] u` as `u[node, subject, re]`
in column-major order within each array element, i.e.
`u[1,1,1], u[1,1,2], ..., u[1,J,n_re], u[2,1,1], ...`.
