# Simulate data from a multilevel VAR model

Generates data from the generative model implied by each Stan model
family. Useful for testing parameter recovery and model validation.

## Usage

``` r
sim_var(
  N,
  T_obs,
  p,
  K = 1L,
  family = c("bernoulli", "ordinal", "gaussian"),
  alpha = NULL,
  gamma = NULL,
  Phi = NULL,
  sigma = NULL,
  kappa = NULL,
  q = 0L,
  x_gen = NULL,
  sd_alpha = 0.5,
  sd_phi = 0.2,
  sd_gamma = NULL,
  re_temporal = FALSE,
  C = 5L,
  burnin = 500L,
  seed = NULL
)
```

## Arguments

- N:

  Integer. Number of subjects (groups).

- T_obs:

  Integer. Number of time points per subject.

- p:

  Integer. Number of outcome nodes.

- K:

  Integer. AR order (default 1).

- family:

  Character. One of `"bernoulli"`, `"ordinal"`, `"gaussian"`.

- alpha:

  Numeric vector of length `p`. Population intercepts (on logit scale
  for bernoulli, identity for gaussian). For ordinal, this is absorbed
  into kappa and should be left NULL. Generated if NULL.

- gamma:

  Matrix `q x p`. Population covariate effects. Generated if NULL and
  `q > 0`.

- Phi:

  Matrix `(p*K) x p`. Population lag coefficients. Generated if NULL.

- sigma:

  Numeric vector of length `p`. Residual SD per node (gaussian only).
  Generated if NULL.

- kappa:

  List of `p` ordered vectors, each of length `C-1`. Cutpoints per node
  (ordinal only). Generated if NULL.

- q:

  Integer. Number of covariates (default 0).

- x_gen:

  Function `f(N, T_obs)` returning an `N x T_obs x q` array of
  covariates. If NULL, default generation is used.

- sd_alpha:

  Numeric. SD of random intercepts (scalar or p-vector). Default 0.5.
  Set to 0 to simulate a fixed-effects-only model with no between-person
  variation in intercepts.

- sd_phi:

  Numeric. SD of random lag coefficients (scalar or matrix). Default
  0.2.

- sd_gamma:

  Numeric or NULL. SD of random covariate slopes. NULL means no random
  slopes on covariates.

- re_temporal:

  Logical. Include random slopes on lag predictors? Default FALSE.

- C:

  Integer. Number of ordinal categories (ordinal only, default 5).

- burnin:

  Integer. Number of time points to discard as warmup before recording
  data (default 500). The VAR process is simulated for `burnin + T_obs`
  time points per subject, and the first `burnin` are discarded. This
  allows the process to reach its stationary distribution before data
  collection begins.

- seed:

  Integer or NULL. RNG seed.

## Value

A list with two components:

- data:

  A long-format data frame with columns `id`, `t`, `y_1`, ..., `y_p`,
  and optionally `x_1`, ..., `x_q`.

- truth:

  A list of true generating parameters.

## Details

To simulate a VAR without any random effects (i.e. all subjects share
identical parameters), set `sd_alpha = 0`, `re_temporal = FALSE` (the
default), and `sd_gamma = NULL` (the default).
