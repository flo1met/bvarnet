# Extract random-effect summaries

Returns random-effect standard deviations (group-level variance),
subject-level posterior means, or the full posterior draws of the
subject-level random effects `u`.

## Usage

``` r
extract_random_effects(object, what = c("sd", "mean_u", "draws_u"))
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

## Value

Depends on `what`; see above.
