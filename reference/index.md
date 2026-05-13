# Package index

## Model Fitting

Estimate Bayesian multilevel VAR models

- [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md) : Fit
  a Bayesian multilevel VAR network model

- [`get_default_priors()`](https://flo1met.github.io/bvarnet/reference/get_default_priors.md)
  : Get the default prior specification for a given model family

- [`prior()`](https://flo1met.github.io/bvarnet/reference/prior.md) :
  Construct a single prior distribution

- [`set_priors()`](https://flo1met.github.io/bvarnet/reference/set_priors.md)
  :

  Build a prior specification object for
  [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md)

## Parameter Extraction

Extract and summarize model parameters

- [`extract_param()`](https://flo1met.github.io/bvarnet/reference/extract_param.md)
  : Extract labelled parameter summaries from a fitted bvarnet model
- [`extract_draws()`](https://flo1met.github.io/bvarnet/reference/extract_draws.md)
  : Extract raw posterior draws for a single parameter block
- [`extract_network_matrix()`](https://flo1met.github.io/bvarnet/reference/extract_network_matrix.md)
  : Extract a network matrix of temporal coefficients
- [`extract_temporal()`](https://flo1met.github.io/bvarnet/reference/extract_temporal.md)
  : Extract temporal (VAR lag) effects
- [`extract_random_effects()`](https://flo1met.github.io/bvarnet/reference/extract_random_effects.md)
  : Extract random-effect summaries

## Inference

Bayes factors and model comparison

- [`bf_table()`](https://flo1met.github.io/bvarnet/reference/bf_table.md)
  : Computes Savage-Dickey density ratio Bayes factors for each
  parameter in the requested subset and returns a tidy data frame.

## Simulation & Prediction

Simulate data and generate predictions

- [`sim_var()`](https://flo1met.github.io/bvarnet/reference/sim_var.md)
  : Simulate data from a multilevel VAR model

## Summaries

Print and summarise model objects

- [`summary(`*`<bvarnet>`*`)`](https://flo1met.github.io/bvarnet/reference/summary.bvarnet.md)
  : Summary method for bvarnet objects
- [`print(`*`<bvarnet>`*`)`](https://flo1met.github.io/bvarnet/reference/print.bvarnet.md)
  : Print a bvarnet model object
- [`print(`*`<bvarnet_prior>`*`)`](https://flo1met.github.io/bvarnet/reference/print.bvarnet_prior.md)
  : Print a bvarnet_prior
- [`print(`*`<bvarnet_priors>`*`)`](https://flo1met.github.io/bvarnet/reference/print.bvarnet_priors.md)
  : Print a bvarnet_priors specification
- [`print(`*`<summary.bvarnet>`*`)`](https://flo1met.github.io/bvarnet/reference/print.summary.bvarnet.md)
  : Print a bvarnet summary
- [`format(`*`<bvarnet_prior>`*`)`](https://flo1met.github.io/bvarnet/reference/format.bvarnet_prior.md)
  : Format a bvarnet_prior for printing
- [`compare_to_truth()`](https://flo1met.github.io/bvarnet/reference/compare_to_truth.md)
  : Compare fitted model parameters to simulation truth

## Data

Example datasets

- [`studentlife`](https://flo1met.github.io/bvarnet/reference/studentlife.md)
  : StudentLife Data
