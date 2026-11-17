# Changelog

## bvarnet 1.0.3.9000

- Removed `skip_lag` immediately, including its zero-filled gap
  likelihood, which could bias shared parameters and temporal estimates.
  Passing this argument now errors for any value. All model families use
  complete-transition analysis: responses with any unavailable required
  lag are excluded from the likelihood. Observation counts and estimates
  can therefore differ from the former default.
- Added `beep_col`. With `beep_col = NULL`, `time_col` is the scheduled
  integer occasion index within a subject. With both supplied,
  `time_col` labels days or sessions and `beep_col` indexes scheduled
  occasions within them. Temporal lags never cross subject/day
  boundaries or missing scheduled occasions. Subject random effects
  remain shared across days.
- Switching to a completely centered backend. The Stan models now center
  all fixed effects parameter internally. This can lead to improved
  posterior geometry in complex models. This does not change parameter
  estimates or Bayes factors!
- Internal centered blocks are excluded from stored public
  draws/diagnostics.
- Custom initial values are mapped to fit the transformation.

## bvarnet 1.0.3

CRAN release: 2026-09-10

- Bugfix: Bayes factors for interaction effects were sometimes displayed
  wrong.

## bvarnet 1.0.3

CRAN release: 2026-09-10

- Minor bugfix that let to an error when running a model without
  covariates.
- [`extract_draws()`](https://flo1met.github.io/bvarnet/reference/extract_draws.md)
  now accepts several parameters at once.
- Minor documentation fixes
- Bugfix: Removing order constrained from adjacent category model.

## bvarnet 1.0.2

CRAN release: 2026-08-31

### New features

- Introduce toolchain free installation via
  [`bvarnet_setup_models()`](https://flo1met.github.io/bvarnet/reference/bvarnet_setup_models.md).

### Bug fixes affecting results

- Fixed automatic scaling of Gaussian default priors and related
  misscalulation of BFs
- Fixed that `extract_random_effects(what = "mean_u")` and
  `what = "draws_u"` returned scrambled subject estimates and RE
  mislabeling in mixed-family fits.
- Change default CI to 95% and introduce argument to vary CI width.

### Other changes

- Defined the `%||%` operator internally.
- [`print()`](https://rdrr.io/r/base/print.html) now reports the
  effective scale of default Gaussian priors, which are widened by the
  outcome SD before reaching Stan. Previously this scaling was
  invisible.
- Minor bug and documentation fixes.
- [`bvar()`](https://flo1met.github.io/bvarnet/reference/bvar.md) gains
  `...`, forwarding additional arguments (e.g. `init`, `refresh`,
  `thin`, `step_size`) to CmdStanR’s `$sample()` method.
- Added a safeguard against duplicated `(id, time)` rows, which
  previously produced an ambiguous, silently-contaminated lag design;
  these now error with guidance to deduplicate or aggregate.
- `time_col` must now be integer-valued (one time unit = one lag step).

## bvarnet 1.0.1

CRAN release: 2026-06-13

- Add installation safeguard on CRAN.

## bvarnet 1.0.0

CRAN release: 2026-05-27

- Initial CRAN submission.
