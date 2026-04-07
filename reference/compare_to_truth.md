# Compare fitted model parameters to simulation truth

Extracts posterior summaries from a fitted `bvarnet` object and compares
them to the true parameter values used for data generation.

## Usage

``` r
compare_to_truth(
  fit,
  truth,
  ci_width = 0.9,
  bayes_factor = FALSE,
  null_value = 0
)
```

## Arguments

- fit:

  A fitted `bvarnet` object (output from
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md)).

- truth:

  The `truth` component from
  [`sim_var()`](https://flo1met.github.io/bvarnet/reference/sim_var.md)
  output.

- ci_width:

  Numeric. Width of the credible interval (default 0.90).

- bayes_factor:

  Logical; if `TRUE`, compute Savage-Dickey BFs for beta and phi
  parameters and append `BF01`, `BF10`, and `bf_correct` columns.
  `bf_correct` is `TRUE` when BF01 \> 1 for true null parameters (true
  value == `null_value`) and BF10 \> 1 for true non-null parameters.
  Default `FALSE`.

- null_value:

  Numeric scalar; the null hypothesis value for Bayes factor computation
  (default 0). Only used when `bayes_factor = TRUE`.

## Value

A data frame with columns: parameter, node, true_value, post_mean,
post_sd, ci_lower, ci_upper, covered (logical), and optionally BF01,
BF10, bf_correct.
