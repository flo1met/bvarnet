# Get the default prior specification for a given model family

Returns a `bvarnet_priors` object showing the default priors that apply
to a particular model configuration. Parameters irrelevant to the chosen
family or model structure are omitted, so the returned object reflects
what the sampler will actually use.

## Usage

``` r
get_default_priors(family = NULL, has_re = TRUE)
```

## Arguments

- family:

  Character (optional). One of `"bernoulli"`, `"ordinal"`, `"gaussian"`.
  When `NULL` (the default), all parameter priors are shown.

- has_re:

  Logical. Does the model include random effects? Default `TRUE`. When
  `FALSE`, the `sd_u` prior is omitted.

## Value

A `bvarnet_priors` object.
