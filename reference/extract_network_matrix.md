# Extract a network matrix of temporal coefficients

Returns a named `p x p` matrix of posterior summary statistics for the
VAR lag coefficients at a chosen lag, suitable for network visualisation
(e.g., with igraph or qgraph).

## Usage

``` r
extract_network_matrix(
  object,
  lag = 1L,
  stat = c("mean", "median", "q5", "q95")
)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- lag:

  Integer. Which lag block. Default 1.

- stat:

  Character. Summary statistic to fill the matrix with: `"mean"`
  (default), `"median"`, `"q5"`, or `"q95"`.

## Value

A named `p x p` numeric matrix. Element `[i, j]` gives the effect of
variable `i` (lagged) on variable `j` (outcome). Row and column names
are the outcome variable names.
