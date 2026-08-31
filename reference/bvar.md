# Fit a Bayesian multilevel VAR network model

The `bvar` function estimates the posterior distribution of the
specified Bayesian (Multilevel) Vector Autoregression.

## Usage

``` r
bvar(
  id_col,
  time_col,
  y_cols,
  x_cols = NULL,
  center_x = FALSE,
  fe_interactions = NULL,
  re_interactions = NULL,
  re_cols = NULL,
  re_temporal = FALSE,
  K = 1,
  na_action = c("listwise"),
  skip_lag = TRUE,
  data,
  family = c("bernoulli", "ordinal", "gaussian"),
  priors = set_priors(),
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 1,
  seed = NULL,
  adapt_delta = NULL,
  max_treedepth = NULL,
  save_data = FALSE,
  ...
)
```

## Arguments

- id_col:

  Character. Name of the subject/group identifier column.

- time_col:

  Character. Name of the time column. Must be integer-valued (one time
  unit = one lag step); non-integer values error.

- y_cols:

  Character vector. Names of the outcome columns.

- x_cols:

  Character vector or NULL. Names of the covariate columns.

- center_x:

  Logical. Grand-mean centre covariates before fitting? Default `FALSE`.

- fe_interactions:

  List or NULL. Fixed-effect interaction terms to add to the design
  matrix. Each element is a character vector of column names to
  interact, or `c("lag", "x")` to interact all lag columns with a
  covariate.

- re_interactions:

  List or NULL. Random-effect interaction terms.

- re_cols:

  Character vector. Columns from X and/or "Intercept" to include as
  random slopes.

- re_temporal:

  Logical. Include random slopes on lag predictors? Default `FALSE`.

- K:

  Integer. AR order. Default 1.

- na_action:

  Character. Missing-data strategy; currently only `"listwise"`.

- skip_lag:

  Logical. If `TRUE` (default), rows with irregular time gaps have their
  lag predictors zero-filled rather than being dropped. For this,
  `time_col` should be scaled so that each step has a length of one
  (e.g., integer days, weeks, or months). Non-integer values error.

- data:

  Data frame in long format.

- family:

  Character scalar or vector. Observation model per node. A scalar is
  recycled to all `y_cols`. A vector of length `length(y_cols)` (named
  or positional) specifies per-node families. Valid values:
  `"bernoulli"`, `"ordinal"`, `"gaussian"`.

- priors:

  A `bvarnet_priors` object from
  [`set_priors()`](https://flo1met.github.io/bvarnet/reference/set_priors.md).
  Defaults to
  [`set_priors()`](https://flo1met.github.io/bvarnet/reference/set_priors.md)
  (package defaults).

- iter:

  Integer. Number of post-warmup iterations per chain. Default 4000.

- warmup:

  Integer. Number of warmup iterations per chain. Default 1000.

- chains:

  Integer. Number of MCMC chains. Default 4.

- cores:

  Integer. Number of chains to run in parallel. Default 1.

- seed:

  Integer or NULL. RNG seed.

- adapt_delta:

  Numeric in (0, 1). Target average proposal acceptance probability
  during warmup adaptation. Higher values (e.g., 0.95–0.99) reduce
  divergences at the cost of slower sampling. Default `NULL` (CmdStan
  default of 0.8).

- max_treedepth:

  Integer. Maximum depth of the NUTS binary tree. Increasing this allows
  the sampler to take more leapfrog steps per iteration, which can help
  with difficult posteriors (e.g., funnels in hierarchical logistic
  models) but increases computation. Default `NULL` (CmdStan default of
  10).

- save_data:

  Logical. If `TRUE`, store the preprocessed (sorted, listwise-deleted)
  estimation data in the `data_used` slot of the returned object for
  reproducibility and downstream analyses. Default `FALSE`.

- ...:

  Additional named arguments forwarded to the CmdStanR `$sample()`
  method, e.g. `init`, `refresh`, `thin`, `step_size`, or
  `show_messages`. See
  [`?cmdstanr::"model-method-sample"`](https://mc-stan.org/cmdstanr/reference/model-method-sample.html)
  for the full list. Arguments that `bvar()` sets itself (`data`,
  `seed`, `iter_warmup`, `iter_sampling`, `chains`, `parallel_chains`,
  `adapt_delta`, `max_treedepth`, and CmdStanR's deprecated aliases for
  them) are rejected with an error; use the corresponding `bvar()`
  argument instead.

## Value

A `bvarnet` object (a named list) with slots: `draws`, `convergence`,
`diagnostics`, `timing`, `metadata`, `return_codes`, `family`,
`standata`, `priors`, `priors_effective`. If `save_data = TRUE`, also
includes `data_used` (the cleaned estimation data frame).

## Requested versus effective priors

`priors` holds the priors as you specified them. `priors_effective`
holds the priors Stan actually sampled under: one `bvarnet_priors`
object per outcome. The two differ for Gaussian outcomes left at their
default priors, where `intercept`, `beta` and `sigma` scales are
multiplied by the outcome SD so that unit-scale defaults stay weakly
informative on the raw data scale. User-supplied priors are never
rescaled.

Bayes factors
([`bf_table`](https://flo1met.github.io/bvarnet/reference/bf_table.md))
divide by `priors_effective`, as the Savage-Dickey density ratio
requires. Note that the joint path (`family` a single value) scales
every outcome by the mean SD across outcomes, whereas the nodewise path
(mixed `family`) scales each outcome by its own SD, so the same data can
imply slightly different effective priors depending on which path runs.

## See also

[`bvarnet_setup_models`](https://flo1met.github.io/bvarnet/reference/bvarnet_setup_models.md),
which must be run once before the first `bvar()` call to set up the
required Stan models (either by downloading precompiled binaries or
compiling them locally).

## Examples

``` r
if (FALSE) { # \dontrun{
# Run bvar on studentlife data
data(studentlife, package = "bvarnet")
fit <- bvar(
  id_col = "id",
  time_col = "day",
  y_cols = c("anxious", "calm", "conventional", "critical", "dependable"),
  re_temporal = TRUE,
  K = 1,
  data = studentlife,
  family = "ordinal",
  priors = set_priors(),
  seed = 1337)

summary(fit)
} # }
```
