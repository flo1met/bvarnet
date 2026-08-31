# Extract random-effect summaries

Returns random-effect standard deviations (group-level variance),
subject-level posterior means, or the full posterior draws of the
subject-level random effects `u`.

## Usage

``` r
extract_random_effects(
  object,
  what = c("sd", "mean_u", "draws_u"),
  ci_level = 0.95
)
```

## Arguments

- object:

  A `bvarnet` object returned by
  [`bvar`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- what:

  Character. What to extract:

  `"sd"`

  :   Data frame of random-effect SD summaries (from `extract_param`).

  `"mean_u"`

  :   3D array `[node, subject, re]` of posterior means of subject-level
      effects.

  `"draws_u"`

  :   4D array `[draw, node, subject, re]` of full posterior draws.

- ci_level:

  Numeric scalar strictly between 0 and 1; the mass of the equal-tailed
  credible interval reported in `ci_lower` and `ci_upper`. Default
  `0.95`. Only used when `what = "sd"`.

## Value

Depends on `what`; see above.

## Array layout

The `"mean_u"` and `"draws_u"` arrays are indexed `[node, subject, re]`
(with a leading `draw` dimension for `"draws_u"`), matching the Stan
declaration `array[p] matrix[J, n_re] u`. `node` is named after the
outcome columns, `subject` after the subject index `1..J`, and `re`
after the random-effect design columns. Elements are placed by parsing
the `u[node, subject, re]` indices from the draws, so the layout is the
same whether the model was fitted through the joint or the nodewise
path.
