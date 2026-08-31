# bvarnet 1.0.1.9000

## New features
* Introduce toolchain free installation via `bvarnet_setup_models()`.

## Bug fixes affecting results

* Fixed automatic scaling of Gaussian default priors and related misscalulation of BFs 
* Fixed that `extract_random_effects(what = "mean_u")` and `what = "draws_u"` returned scrambled subject estimates and RE mislabeling in mixed-family fits.
* Change default CI to 95% and introduce argument to vary CI width.


## Other changes

* Defined the `%||%` operator internally.
* `print()` now reports the effective scale of default Gaussian priors, which
  are widened by the outcome SD before reaching Stan. Previously this scaling
  was invisible.
* Minor bug and documentation fixes.
* `bvar()` gains `...`, forwarding additional arguments (e.g. `init`, `refresh`,
  `thin`, `step_size`) to CmdStanR's `$sample()` method. 
* Added a safeguard against duplicated `(id, time)` rows, which previously
  produced an ambiguous, silently-contaminated lag design; these now error
  with guidance to deduplicate or aggregate.
* `time_col` must now be integer-valued (one time unit = one lag step).

# bvarnet 1.0.1

* Add installation safeguard on CRAN.

# bvarnet 1.0.0

* Initial CRAN submission.
