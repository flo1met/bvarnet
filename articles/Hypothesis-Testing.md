# Hypothesis-Testing

## Setup

First, lets load the package:

``` r
library(bvarnet)
library(qgraph)
```

## Data

Now, we can load the example data:

``` r
data(studentlife)
```

There is some missing data in the dataset. The models default options
handle this by themselves. For a further elaboration on this, you can
read
[`vignette("Missing-Data")`](https://flo1met.github.io/bvarnet/articles/Missing-Data.md).

## Model Estimation

This time, well simply use the default priors, for the ordinal model
they are

``` r
fit <- bvar(
  id_col = "id",
  time_col = "day",
  y_cols = c("anxious", "calm", "conventional", "critical", "dependable"),
  x_cols = c("sleep_hour"),
  re_cols = NULL,
  re_temporal = FALSE,
  K = 1,
  data = studentlife,
  family = c("ordinal"),
  priors = set_priors(),
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1337
)
```

## Model Output

Again, lets get the model overview and summary:

``` r
print(fit)
#> BVAR Network fit
#> ======================================== 
#> Family:      ordinal
#> Outcomes (p): 5 
#> Lags (K):     1 
#> Fixed eff.:   1 
#> Observations: 111 
#> Rhat max:    1.001
#> Divergences: 2  WARNING: check model/priors.
#> Priors:       beta ~ Normal(0, 1), phi ~ Normal(0, 0.5), kappa ~ Normal(0, 2) (all defaults)
#> Total time:  11.7 sec
#> ========================================
summary(fit)
#> BVAR Network Summary
#> ================================================== 
#> Family: ordinal | p=5 | K=1 | n=111
#> Rhat max: 1.001 | Divergences: 2
#>   WARNING: divergent transitions detected — check model/priors.
#> 
#> --- Fixed Effect ---
#>  predictor  outcome      mean   median q5     q95   rhat ess_bulk ess_tail
#>  sleep_hour anxious      -0.005 -0.005 -0.060 0.049 1    15771.43 12695.05
#>  sleep_hour calm         -0.032 -0.032 -0.088 0.022 1    17766.75 12604.63
#>  sleep_hour conventional -0.016 -0.015 -0.076 0.042 1    17485.60 12642.45
#>  sleep_hour critical     -0.008 -0.008 -0.056 0.038 1    14049.46 12084.18
#>  sleep_hour dependable   -0.003 -0.003 -0.060 0.053 1    18991.13 12742.34
#> 
#> 
#> --- Autoregressive ---
#>  predictor         outcome      mean  median q5     q95   rhat  ess_bulk ess_tail
#>  lag1_anxious      anxious      0.069 0.070  -0.260 0.395 1.000 16772.27 12193.66
#>  lag1_calm         calm         0.165 0.166  -0.145 0.478 1.000 12365.66 11825.79
#>  lag1_conventional conventional 0.391 0.390   0.036 0.750 1.001 15129.26 12680.89
#>  lag1_critical     critical     0.238 0.235  -0.020 0.503 1.000 18369.60 11582.80
#>  lag1_dependable   dependable   0.427 0.425   0.129 0.733 1.000 16335.00 11476.86
#> 
#> 
#> --- Cross-lagged ---
#>  predictor         outcome      mean   median q5     q95   rhat ess_bulk ess_tail
#>  lag1_calm         anxious      -0.249 -0.247 -0.574 0.070 1    12852.77 11468.22
#>  lag1_conventional anxious      -0.028 -0.027 -0.370 0.314 1    13704.43 12188.62
#>  lag1_critical     anxious       0.296  0.295  0.017 0.582 1    19347.83 13118.49
#>  lag1_dependable   anxious       0.085  0.084 -0.194 0.365 1    16984.34 12402.25
#>  lag1_anxious      calm         -0.028 -0.030 -0.343 0.290 1    16803.47 12540.98
#>  lag1_conventional calm         -0.199 -0.200 -0.542 0.146 1    13545.69 11926.93
#>  lag1_critical     calm         -0.173 -0.172 -0.453 0.099 1    19113.18 11669.83
#>  lag1_dependable   calm          0.167  0.164 -0.098 0.438 1    16138.68 11770.26
#>  lag1_anxious      conventional -0.313 -0.313 -0.663 0.032 1    16432.09 12229.12
#>  lag1_calm         conventional -0.063 -0.062 -0.396 0.266 1    13903.03 12044.49
#> 
#> ... 10 more rows. Use extract_temporal(fit, effect = "cl") for full output.
#> 
#> --- Threshold ---
#>  predictor               outcome mean   median q5     q95    rhat  ess_bulk ess_tail
#>  kappa(anxious, c1)      —       -1.392 -1.384 -2.106 -0.710 1.000 10840.42 10882.22
#>  kappa(calm, c1)         —       -1.787 -1.745 -2.665 -1.035 1.000 13231.26 11291.59
#>  kappa(conventional, c1) —       -1.553 -1.527 -2.347 -0.847 1.000 14172.12 10875.45
#>  kappa(critical, c1)     —        0.026  0.035 -0.497  0.524 1.001 12451.60 11045.22
#>  kappa(dependable, c1)   —       -2.043 -1.991 -3.108 -1.139 1.001 11713.93 10367.08
#>  kappa(anxious, c2)      —        0.405  0.406 -0.163  0.966 1.000 16633.22 13007.72
#>  kappa(calm, c2)         —       -1.215 -1.204 -1.855 -0.619 1.000 17831.13 12444.68
#>  kappa(conventional, c2) —       -1.079 -1.075 -1.709 -0.463 1.000 16960.92 12615.81
#>  kappa(critical, c2)     —        0.423  0.423 -0.034  0.889 1.000 16313.25 13806.83
#>  kappa(dependable, c2)   —       -1.082 -1.076 -1.759 -0.435 1.000 16347.72 12825.94
#> 
#> ... 10 more rows. Use extract_param(fit, type = "Threshold") for full output.
#> 
#> ==================================================
#> Use extract_param() or extract_param(fit, type = "...") for the full parameter table.
#> Use extract_network_matrix() for the temporal network matrix.
```

## Hypothesis Testing with Bayes Factors

To estimate Bayes factors for a single parameter, or a set of
parameters, we can use the
[`bf_table()`](https://flo1met.github.io/bvarnet/reference/bf_table.md)
function. Bayes factors are computed via the **Savage–Dickey density
ratio** (SDDR), which, for a point null $H_{0}:\theta = \theta_{0}$
nested in $H_{1}:\theta \neq \theta_{0}$ with shared nuisance
parameters, gives

$${BF}_{01}\; = \;\frac{p\left( \theta = \theta_{0} \mid y \right)}{p\left( \theta = \theta_{0} \right)},$$

i.e. the ratio of the posterior to the prior density evaluated at the
null value.
[`bf_table()`](https://flo1met.github.io/bvarnet/reference/bf_table.md)
reports ${BF}_{10} = 1/{BF}_{01}$ by approximating the posterior density
at $\theta_{0} = 0$ from the MCMC draws.

To perform a hypothesis test on effects involving a specific variable,
we can call the
[`bf_table()`](https://flo1met.github.io/bvarnet/reference/bf_table.md)
function, specifying which variable we want to look at using the
`variable` argument. Lets investigate if the hours of sleep the
individual got in the previous night (`sleep_hour`) has an influence on
the network structure:

``` r
bf_tab <- bf_table(fit, variable = "sleep_hour")
bf_tab
#>                   type  predictor      outcome    BF10
#> 1         Fixed Effect sleep_hour      anxious 0.03386
#> 2         Fixed Effect sleep_hour         calm 0.05187
#> 3         Fixed Effect sleep_hour conventional 0.03760
#> 4         Fixed Effect sleep_hour     critical 0.02951
#> 5         Fixed Effect sleep_hour   dependable 0.03403
#> 6 Fixed Effect (joint) sleep_hour            — 0.00000
```

The output now shows 6 rows. The five “Fixed Effect” rows show if the
the covariate `sleep_hour` has an influence on each node separately. The
“Fixed Effect (join)” row displays a single Bayes factor that in this
case tests the hypothesis

H: Is there a difference in the baseline of the network structure that
is explained by the variable `sleep_hour`

As we can see, the `BF_{10} = bf_tab[]`. Therefore we can conlude that
the variable `sleep_hour` does not influence the baseline level of the
network.

## Joint Bayes Factors

For more elaborate analyses, or to report more results at once, we can
get a table with all single and joint Bayes factors. For this we can
call the `bf_table(fit)` function

``` r
bf_tab <- bf_table(fit)
bf_tab
#>                      type         predictor      outcome     BF10
#> 1          Autoregressive      lag1_anxious      anxious  0.42563
#> 2          Autoregressive         lag1_calm         calm  0.54824
#> 3          Autoregressive lag1_conventional conventional  2.23181
#> 4          Autoregressive     lag1_critical     critical  1.01615
#> 5          Autoregressive   lag1_dependable   dependable  6.07838
#> 6  Autoregressive (joint)            all_ar            —  2.49597
#> 7            Cross-lagged         lag1_calm      anxious  0.85207
#> 8            Cross-lagged lag1_conventional      anxious  0.41647
#> 9            Cross-lagged     lag1_critical      anxious  1.52971
#> 10           Cross-lagged   lag1_dependable      anxious  0.38114
#> 11           Cross-lagged      lag1_anxious         calm  0.38559
#> 12           Cross-lagged lag1_conventional         calm  0.66942
#> 13           Cross-lagged     lag1_critical         calm  0.55497
#> 14           Cross-lagged   lag1_dependable         calm  0.54173
#> 15           Cross-lagged      lag1_anxious conventional  1.28950
#> 16           Cross-lagged         lag1_calm conventional  0.40758
#> 17           Cross-lagged     lag1_critical conventional  0.39803
#> 18           Cross-lagged   lag1_dependable conventional  0.36941
#> 19           Cross-lagged      lag1_anxious     critical  0.38522
#> 20           Cross-lagged         lag1_calm     critical  0.33946
#> 21           Cross-lagged lag1_conventional     critical  0.37277
#> 22           Cross-lagged   lag1_dependable     critical  0.49510
#> 23           Cross-lagged      lag1_anxious   dependable  2.31460
#> 24           Cross-lagged         lag1_calm   dependable  0.46823
#> 25           Cross-lagged lag1_conventional   dependable 24.44204
#> 26           Cross-lagged     lag1_critical   dependable  0.73579
#> 27   Cross-lagged (joint)            all_cl            —  0.00333
#> 28       Temporal (joint)           all_phi            —  0.00000
#> 29           Fixed Effect        sleep_hour      anxious  0.03386
#> 30           Fixed Effect        sleep_hour         calm  0.05187
#> 31           Fixed Effect        sleep_hour conventional  0.03760
#> 32           Fixed Effect        sleep_hour     critical  0.02951
#> 33           Fixed Effect        sleep_hour   dependable  0.03403
#> 34   Fixed Effect (joint)        sleep_hour            —  0.00000
```

The table can be used to easily get a set of Bayes factors. If we have
hypotheses concerning the temporal structure, we can call the
`bf_table(fit)` function using the argument `type = "temporal"`.

``` r
bf_temp <- bf_table(fit, type = "temporal")
bf_temp
#>               type predictor outcome BF10
#> 1 Temporal (joint)   all_phi       —    0
```

The following arguments are available as filters for the `type`
argument: “ar”, “cl”, “intercepts”, “fe”, “lag_fe”, and “temporal”.
