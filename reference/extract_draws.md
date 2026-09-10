# Extract raw posterior draws for one or more parameter blocks

Returns an `(iterations * chains)` by `params` matrix with Stan-indexed
column names (e.g. `"beta[1,1]"`, `"phi[2,3]"`). Several blocks can be
requested at once; their columns are returned side by side in one
matrix.

## Usage

``` r
extract_draws(object, parameter = "all")
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- parameter:

  Character vector naming one or more of `"beta"`, `"phi"`, `"sd_u"`,
  `"sigma"`, `"kappa"`, or `"lp__"`; or the single value `"all"` (the
  default), which returns every block the fitted model actually has
  draws for.

## Value

A numeric matrix with one row per posterior draw and one column per Stan
parameter element, named with the Stan index (e.g. `"phi[1,2]"`). When
several blocks are requested, their columns appear in the order the
blocks were named. Use
[`extract_param`](https://flo1met.github.io/bvarnet/reference/extract_param.md)
for a labelled summary table instead of raw draws.

## Details

The blocks are the fixed effects `beta`, the temporal coefficients
`phi`, the random-effect SDs `sd_u`, the residual SDs `sigma`, the
ordinal thresholds `kappa`, and `lp__`, the sampler's log density.

## Examples

``` r
if (FALSE) { # \dontrun{
# One block
phi <- extract_draws(fit, "phi")

# Several blocks, side by side in one matrix
d <- extract_draws(fit, c("phi", "kappa"))
colnames(d)

# Everything this model has
d <- extract_draws(fit)
} # }
```
