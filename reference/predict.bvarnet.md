# Predict from a fitted bvarnet model

Computes one-step-ahead or recursive forecasts for long-format
time-series data. Supports population-level (`subject_re = "zero"`) and
subject-specific (`subject_re = "posterior-mean"`) predictions. Also
serves as the out-of-sample engine: fit on training data, call
`predict(fit, newdata = test_data)`.

## Usage

``` r
# S3 method for class 'bvarnet'
predict(
  object,
  newdata = NULL,
  forecast = c("one-step", "recursive"),
  conditioning_window = NULL,
  type = c("link", "response", "probabilities"),
  method = c("posterior-mean", "posterior-sample"),
  ndraws = NULL,
  seed = NULL,
  subject_re = c("zero", "posterior-mean", "sample"),
  new_subject = c("zero", "sample"),
  ...
)
```

## Arguments

- object:

  A `bvarnet` object from
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md).

- newdata:

  Data frame in long format. If `NULL`, the original training data
  design matrices (stored in `object$standata`) are used for in-sample
  fitted values.

- forecast:

  Character. `"one-step"` (default) uses observed lag values from
  `newdata` for every row. `"recursive"` feeds predicted values back
  into the lag buffer after the conditioning window.

- conditioning_window:

  Integer scalar, named integer vector, or `NULL`. Number of observed
  time points (per subject) used to initialise the lag buffer before
  recursive forecasting begins. Must be `>= K`. If `NULL`, defaults to
  `K` (minimum lag history). Ignored when `forecast = "one-step"`.

- type:

  Character. Output type: `"link"` (linear predictor), `"response"`
  (mean on the outcome scale), or `"probabilities"` (category-level
  probabilities/mean+sd).

- method:

  Character. `"posterior-mean"` uses posterior means only
  (deterministic). `"posterior-sample"` averages over `ndraws` draws and
  returns `attr(,"sd")` with across-draw SD.

- ndraws:

  Integer. Number of posterior draws to use when
  `method = "posterior-sample"`. Defaults to the smaller of 100 and the
  total available draws.

- seed:

  Integer or NULL. RNG seed for draw selection and new-subject sampling.

- subject_re:

  Character. How to handle random effects: `"zero"` (population-level, u
  = 0), `"posterior-mean"` (posterior mean of u for seen subjects), or
  `"sample"` (draw-specific u for seen subjects).

- new_subject:

  Character. Fallback for unseen subjects: `"zero"` (u = 0) or
  `"sample"` (draw from RE distribution).

- ...:

  Ignored.

## Value

For `type = "link"` or `"response"`: a numeric matrix with rows matching
the original data (or `nrow(newdata)`), with `NA` for the first K rows
per subject. For `type = "probabilities"`: a list of `p` matrices. When
`method = "posterior-sample"`, the output carries `attr(,"sd")` and
`attr(,"ndraws")`.

## Examples

``` r
if (instantiate::stan_cmdstan_exists()) {
  sim <- sim_var(N = 5, T_obs = 30, p = 2, family = "gaussian", seed = 1)
  fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "gaussian",
              iter = 200, warmup = 100, chains = 2, seed = 1)
  preds <- predict(fit, type = "response")
}
```
