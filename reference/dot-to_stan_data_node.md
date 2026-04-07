# Build Stan data for a single node from shared matrices

Build Stan data for a single node from shared matrices

## Usage

``` r
.to_stan_data_node(shared, node, family, priors)
```

## Arguments

- shared:

  List from
  [`.to_stan_data_shared()`](https://flo1met.github.io/bvarnet/reference/dot-to_stan_data_shared.md).

- node:

  Integer, which node (1..p).

- family:

  Character scalar, family for this node.

- priors:

  A `bvarnet_priors` object (original, unmodified).

## Value

Named list ready for `CmdStanModel$sample()`.
