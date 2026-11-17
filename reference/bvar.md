# Fit a Bayesian multilevel VAR network model

The `bvar` function estimates the posterior distribution of the
specified Bayesian (Multilevel) Vector Autoregression.

## Usage

``` r
bvar(
  id_col,
  time_col,
  beep_col = NULL,
  y_cols,
  x_cols = NULL,
  center_x = FALSE,
  fe_interactions = NULL,
  re_interactions = NULL,
  re_cols = NULL,
  re_temporal = FALSE,
  K = 1,
  na_action = c("listwise"),
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
  save_data = TRUE,
  ...
)
```

## Arguments

- id_col:

  Character. Name of the subject/group identifier column.

- time_col:

  Character scalar naming the time-index column in `data`. Its meaning
  depends on `beep_col`:

  - If `beep_col = NULL` (default), this column contains the scheduled
    occasion index for one sequence per subject. Values must be finite
    and integer-valued, with one unit representing one scheduled lag
    step. Indices must continue across the subject's sequence rather
    than restart each day. For lag k, the predecessor must have index
    `time - k`. No separate day or overnight boundary is inferred in
    this mode.

  - If `beep_col` names a column, `time_col` identifies the day or
    session. `beep_col` then contains finite integer-valued beep indices
    within each day/session, with one unit representing one scheduled
    lag step. Beep indices may restart in each day/session. For lag k,
    the predecessor must have index `beep - k` in the same subject and
    day/session.

  Preserve the original scheduled labels when observations are missing.
  Rows without a complete, correctly spaced lag window are excluded as
  responses from the likelihood; their observed outcomes may still be
  used as predictors for later eligible responses.

- beep_col:

  NULL (default), or a character scalar naming the scheduled occasion
  index within each day/session, specified in the `time_col` column. If
  NULL, `time_col` itself supplies the occasion index for one sequence
  per subject, and no day/session splitting is performed. If supplied,
  `time_col` must identify the day/session and this column must contain
  finite integer-valued beep indices, with one unit representing one
  scheduled lag step. Beep indices may restart in each day/session. For
  lag k, the predecessor must have index `beep - k` in the same subject
  and day/session.

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

  Logical. Keep the Stan data used for estimation in the returned
  object's `standata` slot? Default `TRUE`. Set to `FALSE` to remove
  these data from the fitted model before saving or sharing it, for
  example when working with sensitive data. Model dimensions, variable
  labels, posterior draws and effective priors are retained; this option
  does not anonymize posterior results.

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
`model_info`, `priors`, and `priors_effective`. The `standata` slot
contains the estimation data when `save_data = TRUE` (default) and is
`NULL` when `save_data = FALSE`. The `model_info` slot retains only
dimensions and variable/interaction labels needed for summaries and
parameter extraction. No additional data frame is stored.

## Scheduled time and missing observations

Index columns identify eligible transitions; they do not automatically
add time or day effects as covariates. A daily series with one
observation per day can use `time_col = "day"` alone: consecutive days
are ordinary lag steps. With multiple observations per day, supply both
day and beep columns to exclude overnight transitions. Day/session
labels may be numeric, character, factor or Date; occasion indices must
be numeric integers. Retained keys must be unique within subject (and
day/session when supplied), with no missing keys.

Only complete lag windows enter the likelihood. For indices
`1, 2, 4, 5`, a VAR(1) uses `1 -> 2` and `4 -> 5`, excluding `2 -> 4`.
With day/beep indices, day 1 beep 5 never predicts day 2 beep 1, while
day 2 beep 1 can predict day 2 beep 2. Subject random effects remain
shared across that subject's days. Occasion-level listwise deletion is
applied to all selected outcomes and covariates before constructing
lags. This is not a general correction for missingness under MAR or
MNAR.

The former `skip_lag` argument and zero-filled gap likelihood were
removed because they could bias temporal and non-temporal parameters.
Passing the removed argument is an error. Fits may use fewer responses
and give different estimates than the former default.

## Internal predictor centering

Stan samples intercepts (Gaussian/Bernoulli) or thresholds (ordinal)
using centered completed fixed-effect and lag columns. Reported `beta`,
`kappa`, priors and Bayes-factor hypotheses retain the reference of the
public design produced with the requested `center_x` setting.
Random-effect designs and priors are unchanged. This coordinate change
can improve sampling, but does not guarantee a speed improvement.

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
# Scheduled indices retain missing occasions rather than compressing time.
single_index <- data.frame(id = 1, occasion = c(1, 2, 4, 5), y = c(1, 2, 3, 4))
# time_col = "occasion": transitions 1 -> 2 and 4 -> 5 are eligible.
day_beep <- data.frame(id = 1, day = c(1, 1, 1, 2, 2),
                       beep = c(1, 2, 4, 1, 2), y = c(1, 2, 3, 4, 5))
# time_col = "day", beep_col = "beep": only within-day beeps 1 -> 2.
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
