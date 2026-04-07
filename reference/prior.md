# Construct a single prior distribution

Builds a `bvarnet_prior` object specifying the prior family and its
parameters. Supported families in Phase 1 are `"normal"`, `"student_t"`,
and `"cauchy"`.

## Usage

``` r
prior(family, loc = 0, scale = 1, df = 7)
```

## Arguments

- family:

  Character. One of `"normal"`, `"student_t"`, `"cauchy"`.

- loc:

  Location parameter (default 0).

- scale:

  Scale parameter (default 1). Must be \> 0.

- df:

  Degrees of freedom for `"student_t"` (default 7). Must be0 when
  `family = "student_t"`.

## Value

A `bvarnet_prior` S3 object.
