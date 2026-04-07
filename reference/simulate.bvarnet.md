# Simulate new trajectories from a fitted bvarnet model

Generates new multilevel VAR trajectories using the posterior
parameters. Reuses the existing response kernels from
[`sim_var()`](https://flo1met.github.io/bvarnet/reference/sim_var.md).

## Usage

``` r
# S3 method for class 'bvarnet'
simulate(
  object,
  nsim = 20L,
  seed = NULL,
  method = c("posterior-mean", "posterior-sample"),
  ndraws = 10L,
  N = NULL,
  burnin = 200L,
  x_gen = NULL,
  subject_re = c("sample", "zero"),
  ...
)
```

## Arguments

- object:

  A `bvarnet` object from
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- nsim:

  Integer. Number of time points per subject. Default 20.

- seed:

  Integer or NULL. RNG seed.

- method:

  Character. `"posterior-mean"` returns one data frame.
  `"posterior-sample"` returns a list of data frames (one per draw).

- ndraws:

  Integer. Number of posterior draws when `method = "posterior-sample"`.
  Default 10.

- N:

  Integer. Number of subjects. Defaults to J from the fitted model.

- burnin:

  Integer. Burn-in time points to discard. Default 200.

- x_gen:

  Function `f(N, T)` returning an `N x T x q` array of covariates, or
  NULL for zero covariates.

- subject_re:

  Character. `"sample"` draws REs from `N(0, sd_u)`. `"zero"` sets all
  REs to 0.

- ...:

  Ignored.

## Value

For `method = "posterior-mean"`: a data frame with columns `id`, `t`,
`y_1`, ..., `y_p` (and optional `x_1`, ..., `x_q`). For
`"posterior-sample"`: a list of such data frames, one per draw.

## Examples

``` r
if (instantiate::stan_cmdstan_exists()) {
  sim <- sim_var(N = 5, T_obs = 30, p = 2, family = "gaussian", seed = 1)
  fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "gaussian",
              iter = 200, warmup = 100, chains = 2, seed = 1)
  sim_data <- simulate(fit, nsim = 20, seed = 42)
}
```
