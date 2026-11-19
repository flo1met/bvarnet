# Missing-Data

## Setup

First, lets load the package:

``` r
library(bvarnet)
```

## Data

Now, we can load the example data:

``` r
data(studentlife)
```

## Missing Data in `bvarnet`

Inspecting the dataset, we can see that there are quite a lot of missing
values:

``` r
head(studentlife)
#> # A tibble: 6 × 76
#>      id   day happy happyornot   sad sadornot social_number stress_level sleep_hour sleep_quality difficult_stay_awake anxious
#>   <dbl> <dbl> <dbl>      <int> <dbl>    <dbl>         <dbl>        <dbl>      <dbl>         <dbl>                <int>   <dbl>
#> 1     0     1    NA         NA    NA       NA         NA               3       7.5           1.5                     0      NA
#> 2     0     2    NA         NA    NA       NA          4               2       5.33          2.33                    0      NA
#> 3     0     4    NA         NA    NA       NA          2               4       9             2                       0      NA
#> 4     0     5    NA         NA    NA       NA          2              NA       3             4                       0      NA
#> 5     0     6    NA         NA    NA       NA         NA               1      NA            NA                      NA      NA
#> 6     0     7    NA         NA    NA       NA          2.67            2       5             3                       1      NA
#> # ℹ 64 more variables: calm <dbl>, conventional <dbl>, critical <dbl>, dependable <dbl>, disorganized <dbl>, enthusiastic <dbl>,
#> #   open_to_experiences <dbl>, reserved <dbl>, sympathetic <dbl>, act_running_ep_0 <dbl>, act_running_ep_1 <dbl>,
#> #   act_running_ep_2 <dbl>, act_running_ep_3 <dbl>, act_running_ep_4 <dbl>, act_still_ep_0 <dbl>, act_still_ep_1 <dbl>,
#> #   act_still_ep_2 <dbl>, act_still_ep_3 <dbl>, act_still_ep_4 <dbl>, act_unknown_ep_0 <dbl>, act_unknown_ep_1 <dbl>,
#> #   act_unknown_ep_2 <dbl>, act_unknown_ep_3 <dbl>, act_unknown_ep_4 <dbl>, act_walking_ep_0 <dbl>, act_walking_ep_1 <dbl>,
#> #   act_walking_ep_2 <dbl>, act_walking_ep_3 <dbl>, act_walking_ep_4 <dbl>, act_on_foot_ep_0 <dbl>, act_on_foot_ep_1 <dbl>,
#> #   act_on_foot_ep_2 <dbl>, act_on_foot_ep_3 <dbl>, act_on_foot_ep_4 <dbl>, audio_convo_duration_ep_0 <dbl>, …
```

If we look at the subset of variables that we use in the
`Vignette("bvarnet")` c(“anxious”, “calm”, “conventional”, “critical”,
“dependable”), we see there is approximately $`90\%`$ of missing data:

``` r
sum(is.na(studentlife[c("anxious", "calm", "conventional", "critical", "dependable")])) / (5*nrow(studentlife))
#> [1] 0.8988494
```

## Fitting the Model

We fit the model using the usual call

``` r
fit <- bvar(
  id_col = "id",
  time_col = "day",
  y_cols = c("anxious", "calm", "conventional", "critical", "dependable"),
  x_cols = NULL,
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

## Missing Data Handling

Currently, `bvarnet` is handling missing data by restricting estimation
to completem correctly soaced time-series. Before constructing lagged
predictors, `bvarnet` applies listwise deletion at the
measurement-occasion level. Within each subject, a row is removed if any
selected outcome or covariate is missing. The subject is not removed
completely simply because one or more measurement occasions are
incomplete (complete-case analysis). Entirely absent time points and
rows removed through listwise deletion can therefore create gaps in the
remaining time series. After listwise deletion, the gaps are preserved
rather than compressed. To illustrate this we have the following
example: suppose complete observations are available at occasions $`1`$,
$`2`$, $`4`$, and $`5`$. So the observation at occasion $`3`$ is
missing. A VAR($`1`$) model uses the transitions $`1 \rightarrow 2`$ and
$`4 \rightarrow 5`$. It does not treat $`2 \rightarrow 4`$ as a one-step
transition. Occasion $`4`$ is therefore excluded as a response because
occasion $`3`$ is unavailable, but occasion $`4`$’s observed outcomes
can still predict occasion $`5`$. Excluding an observation as a response
therefore does **not** automatically remove it from subsequent lagged
predictors. The scheduled occasion index must be integer-valued, with
one unit representing one lag step. This index is supplied by
`time_col`, or by `beep_col` when separate day/session labels are used.
Original occasion numbers must be retained when measurements are
missing. Listwise deletion is computationally simple, but it reduces the
available sample size and can produce biased estimates when the retained
observations differ systematically from the omitted observations. This
procedure does not impute missing values or estimate relationships
across incomplete lag windows. It reduces the information available for
estimating both temporal and non-temporal parameters, and retaining
correctly spaced transitions does not itself guarantee unbiased
estimates. Results still depend on the missingness process and on how
the retained observations differ from those omitted.

## Model Output

We can print an overview of the the model we just fitted using the
`print(fit)` functions:

``` r
print(fit)
#> BVAR Network fit
#> ======================================== 
#> Family:      ordinal
#> Outcomes (p): 5 
#> Lags (K):     1 
#> Fixed eff.:   0 
#> Observations: 67 
#> Rhat max:    1.001
#> Divergences: 0
#> Priors:       beta ~ Normal(0, 1), phi ~ Normal(0, 0.5), kappa ~ Normal(0, 2) (all defaults)
#> Total time:  1.8 sec
#> ========================================
```

Here we get information about the number of variables and lags, the
number of observations and some first indications if the model did
converge.

To further inspect the model parameters we can use the `summary(fit)`
function:

``` r
summary(fit)
#> BVAR Network Summary
#> ================================================== 
#> Family: ordinal | p=5 | K=1 | n=67
#> Rhat max: 1.001 | Divergences: 0
#> Credible interval: 95% (equal-tailed)
#> 
#> --- Autoregressive ---
#>  predictor         outcome      mean  median ci_lower ci_upper rhat ess_bulk ess_tail
#>  lag1_anxious      anxious      0.335 0.333  -0.017   0.695    1    24028.15 13285.85
#>  lag1_calm         calm         0.264 0.261  -0.094   0.636    1    21285.17 14037.66
#>  lag1_conventional conventional 0.663 0.657   0.290   1.061    1    20118.54 12641.12
#>  lag1_critical     critical     0.360 0.357   0.081   0.658    1    25754.04 12266.11
#>  lag1_dependable   dependable   0.493 0.489   0.192   0.823    1    21332.14 12422.27
#> 
#> 
#> --- Cross-lagged ---
#>  predictor         outcome      mean   median ci_lower ci_upper rhat ess_bulk ess_tail
#>  lag1_calm         anxious      -0.296 -0.292 -0.666   0.057    1    20896.24 13981.84
#>  lag1_conventional anxious       0.216  0.214 -0.142   0.587    1    22586.54 14071.34
#>  lag1_critical     anxious       0.125  0.125 -0.155   0.410    1    27492.43 11988.09
#>  lag1_dependable   anxious       0.116  0.116 -0.181   0.412    1    24005.22 12934.04
#>  lag1_anxious      calm         -0.113 -0.113 -0.484   0.251    1    22827.71 12798.41
#>  lag1_conventional calm         -0.223 -0.220 -0.614   0.155    1    22023.63 13241.97
#>  lag1_critical     calm         -0.107 -0.106 -0.412   0.192    1    26471.74 12015.22
#>  lag1_dependable   calm          0.173  0.170 -0.123   0.492    1    23787.15 13315.81
#>  lag1_anxious      conventional -0.151 -0.147 -0.524   0.214    1    21837.14 13462.33
#>  lag1_calm         conventional -0.035 -0.034 -0.381   0.311    1    22211.31 13923.50
#> 
#> ... 10 more rows. Use extract_temporal(fit, effect = "cl") for full output.
#> 
#> --- Threshold ---
#>  predictor               outcome mean   median ci_lower ci_upper rhat ess_bulk ess_tail
#>  kappa(anxious, c1)      —       -0.697 -0.684 -2.419   0.957    1    21033.84 12530.14
#>  kappa(calm, c1)         —       -2.043 -2.009 -4.306   0.063    1    22862.97 13749.32
#>  kappa(conventional, c1) —        0.180  0.192 -1.784   2.072    1    21798.18 14399.65
#>  kappa(critical, c1)     —        0.006  0.001 -1.644   1.661    1    22687.71 13188.50
#>  kappa(dependable, c1)   —       -1.229 -1.217 -3.194   0.650    1    24006.53 13586.09
#>  kappa(anxious, c2)      —        1.626  1.631 -0.031   3.261    1    21835.94 13374.13
#>  kappa(calm, c2)         —       -1.083 -1.079 -2.858   0.689    1    22016.99 13225.12
#>  kappa(conventional, c2) —       -0.694 -0.695 -2.498   1.073    1    20883.11 14297.24
#>  kappa(critical, c2)     —        1.290  1.282 -0.398   2.997    1    21588.94 13353.84
#>  kappa(dependable, c2)   —       -0.714 -0.707 -2.471   0.989    1    23092.55 13070.41
#> 
#> ... 10 more rows. Use extract_param(fit, type = "Threshold") for full output.
#> 
#> ==================================================
#> Use extract_param() or extract_param(fit, type = "...") for the full parameter table.
#> Use extract_network_matrix() for the temporal network matrix.
```

### Practical implications

Only observations with a complete and correctly spaced lag window enter
the likelihood. This preserves the intended timing of all retained
lagged relationships, but it may substantially reduce the number of
observations available for estimating both temporal and non-temporal
parameters. This procedure does not impute missing values or recovers
temporal information from incomplete lag windows. If listwise deletion
is inappropriate for the missingness process, users should consider a
suitable imputation procedure or an explicit model for the missing data.
