# Random-Effects

## Setup

``` r
library(bvarnet)

# subject to be removed again...
library(qgraph)
```

## Data

``` r
data(studentlife)
```

There is some missing data in the dataset. The models default options
handle this by themselves. For a further elaboration on this, you can
read
[`vignette("Missing-Data")`](https://flo1met.github.io/bvarnet/articles/Missing-Data.md).

## Random-Effects Arguments

For the following tutorial we will extend the network containing the the
variables `anxious`, `calm`, `conventional`, `critical`, and
`dependable` from the `Vignette(bvarnet)`. To introduce a multilevel
structure to our model, we can use the two arguments: `re_cols` and
`re_temporal`.

The `re_cols` argument accepts all variables that we specify in `x_cols`
and `"Intercept"`, to introduce random effects on the baseline model.
The `re_temporal` argument is a binary TRUE/FALSE indicator to introduce
random effects on the temporal structure.

## Model Estimation

To estimate a model with random effects, we can use the following code,
where we specify the `re_cols = c("Intercept")` to introduce random
effects on the baseline model, and `re_temporal = TRUE` to introduce
random effects on the temporal structure:

``` r
fit <- bvar(
  id_col = "id",
  time_col = "day",
  y_cols = c("anxious", "calm", "conventional", "critical", "dependable"),
  x_cols = NULL,
  re_cols = c("Intercept"),
  re_temporal = TRUE,
  K = 1,
  data = studentlife,
  family = c("ordinal"),
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1337
)
```

## Model Output

``` r
print(fit)
#> BVAR Network fit
#> ======================================== 
#> Family:      ordinal
#> Outcomes (p): 5 
#> Lags (K):     1 
#> Fixed eff.:   0 
#> Random eff.:  6 
#> Observations: 147 
#> Rhat max:    1.002
#> Divergences: 0
#> Priors:       beta ~ Normal(0, 1), phi ~ Normal(0, 0.5), sd_u ~ Half-Normal(0, 1), kappa ~ Normal(0, 2) (all defaults)
#> Total time:  55.4 sec
#> ========================================
summary(fit)
#> BVAR Network Summary
#> ================================================== 
#> Family: ordinal | p=5 | K=1 | n=147
#> Rhat max: 1.002 | Divergences: 0
#> 
#> --- Autoregressive ---
#>  predictor         outcome      mean   median q5     q95   rhat  ess_bulk  ess_tail 
#>  lag1_anxious      anxious      -0.147 -0.145 -0.552 0.246 1.001 14144.797 10285.229
#>  lag1_calm         calm         -0.164 -0.163 -0.510 0.179 1.000 11041.313 11301.026
#>  lag1_conventional conventional  0.100  0.100 -0.328 0.523 1.000 14948.680 11292.342
#>  lag1_critical     critical      0.021  0.043 -0.478 0.450 1.000  7627.365  9737.259
#>  lag1_dependable   dependable    0.433  0.432  0.077 0.798 1.001 12392.768 11687.742
#> 
#> 
#> --- Cross-lagged ---
#>  predictor         outcome      mean   median q5     q95   rhat ess_bulk  ess_tail 
#>  lag1_calm         anxious      -0.055 -0.054 -0.416 0.300 1    11921.220 11122.966
#>  lag1_conventional anxious      -0.151 -0.150 -0.570 0.265 1    12262.803 11589.538
#>  lag1_critical     anxious       0.423  0.429 -0.001 0.833 1     9948.882  9928.583
#>  lag1_dependable   anxious      -0.030 -0.029 -0.377 0.316 1    11970.579 10997.386
#>  lag1_anxious      calm          0.009  0.006 -0.367 0.390 1    14091.704 11395.360
#>  lag1_conventional calm          0.109  0.108 -0.277 0.498 1    13148.537 11706.244
#>  lag1_critical     calm         -0.126 -0.125 -0.462 0.209 1    14628.540 11751.291
#>  lag1_dependable   calm          0.173  0.172 -0.146 0.497 1    13038.248 12165.331
#>  lag1_anxious      conventional -0.176 -0.174 -0.599 0.244 1    14940.034 12389.111
#>  lag1_calm         conventional  0.002  0.004 -0.377 0.379 1    12654.900 10321.554
#> 
#> ... 10 more rows. Use extract_temporal(fit, effect = "cl") for full output.
#> 
#> --- Random Effect SD ---
#>  predictor    outcome      mean  median q5    q95   rhat  ess_bulk ess_tail 
#>  anxious      Intercept    1.466 1.441  0.975 2.049 1.001 5771.638  8999.489
#>  calm         Intercept    1.177 1.159  0.779 1.648 1.001 6702.738  9353.987
#>  conventional Intercept    1.126 1.105  0.730 1.595 1.000 7457.640 10383.186
#>  critical     Intercept    1.357 1.326  0.863 1.952 1.001 6335.629  9610.739
#>  dependable   Intercept    0.959 0.937  0.577 1.415 1.000 5733.278  9413.112
#>  anxious      lag1_anxious 0.201 0.167  0.015 0.506 1.000 7103.922  7853.716
#>  calm         lag1_anxious 0.160 0.129  0.012 0.416 1.000 7432.433  7148.949
#>  conventional lag1_anxious 0.251 0.210  0.019 0.616 1.000 6548.701  7779.383
#>  critical     lag1_anxious 0.299 0.258  0.026 0.711 1.001 6084.399  7761.607
#>  dependable   lag1_anxious 0.297 0.242  0.021 0.766 1.000 5121.613  6943.454
#> 
#> ... 20 more rows. Use extract_random_effects(fit) for full output.
#> 
#> --- Threshold ---
#>  predictor               outcome mean   median q5     q95    rhat ess_bulk  ess_tail 
#>  kappa(anxious, c1)      —       -2.101 -2.088 -2.896 -1.368 1     6916.424  9765.659
#>  kappa(calm, c1)         —       -2.218 -2.181 -3.134 -1.428 1     7725.885 10490.939
#>  kappa(conventional, c1) —       -2.032 -1.994 -2.893 -1.289 1     9459.150 10938.668
#>  kappa(critical, c1)     —       -0.265 -0.262 -0.894  0.354 1     6852.727  8909.718
#>  kappa(dependable, c1)   —       -2.502 -2.467 -3.532 -1.586 1     8625.055  8684.235
#>  kappa(anxious, c2)      —        0.443  0.444 -0.191  1.071 1     7571.908 11243.181
#>  kappa(calm, c2)         —       -1.501 -1.494 -2.162 -0.872 1    10199.428 10299.144
#>  kappa(conventional, c2) —       -1.255 -1.247 -1.859 -0.676 1    11832.378 12090.803
#>  kappa(critical, c2)     —        1.164  1.154  0.548  1.816 1     7963.996 10452.848
#>  kappa(dependable, c2)   —       -1.242 -1.232 -1.848 -0.666 1    10986.504 11362.741
#> 
#> ... 10 more rows. Use extract_param(fit) for full output.
#> 
#> ==================================================
#> Use extract_param() for the full parameter table.
#> Use extract_network_matrix() for the temporal network matrix.
```

## Extracting Random Effects

Additionally to the `extract_*` that we already described in
`Vignette(bvarnet)`, we can use the
[`extract_random_effects()`](https://flo1met.github.io/bvarnet/reference/extract_random_effects.md)
function to only extract the random effects:

``` r
re <- extract_random_effects(fit)
re
#>                type    predictor           outcome      mean    median          q5       q95     rhat ess_bulk  ess_tail
#> 1  Random Effect SD      anxious         Intercept 1.4658813 1.4405375 0.974907417 2.0494664 1.000511 5771.638  8999.489
#> 2  Random Effect SD         calm         Intercept 1.1768493 1.1586542 0.779077144 1.6483890 1.000520 6702.738  9353.987
#> 3  Random Effect SD conventional         Intercept 1.1263609 1.1046796 0.729782739 1.5952215 1.000058 7457.640 10383.186
#> 4  Random Effect SD     critical         Intercept 1.3569946 1.3259133 0.863404840 1.9516987 1.000630 6335.629  9610.739
#> 5  Random Effect SD   dependable         Intercept 0.9587509 0.9369383 0.577423646 1.4146214 1.000467 5733.278  9413.112
#> 6  Random Effect SD      anxious      lag1_anxious 0.2014235 0.1670989 0.015447781 0.5063193 1.000359 7103.922  7853.716
#> 7  Random Effect SD         calm      lag1_anxious 0.1599548 0.1290305 0.012097991 0.4159594 1.000280 7432.433  7148.949
#> 8  Random Effect SD conventional      lag1_anxious 0.2510272 0.2100882 0.019368192 0.6161868 1.000499 6548.701  7779.383
#> 9  Random Effect SD     critical      lag1_anxious 0.2985005 0.2584240 0.026466864 0.7112709 1.000631 6084.399  7761.607
#> 10 Random Effect SD   dependable      lag1_anxious 0.2966100 0.2420921 0.020789221 0.7656729 1.000356 5121.613  6943.454
#> 11 Random Effect SD      anxious         lag1_calm 0.1400964 0.1116320 0.010623852 0.3702453 1.000092 8388.615  8389.123
#> 12 Random Effect SD         calm         lag1_calm 0.1717285 0.1456824 0.014364122 0.4220376 1.000588 5932.985  7207.692
#> 13 Random Effect SD conventional         lag1_calm 0.1616909 0.1303977 0.012481142 0.4176642 1.000347 6438.309  7672.721
#> 14 Random Effect SD     critical         lag1_calm 0.1378683 0.1082645 0.009859052 0.3682334 1.000254 8913.645  8443.741
#> 15 Random Effect SD   dependable         lag1_calm 0.1977470 0.1628938 0.015441451 0.5053782 1.001017 5239.681  7090.355
#> 16 Random Effect SD      anxious lag1_conventional 0.1962377 0.1649739 0.015909338 0.4803717 1.000489 6342.986  7313.237
#> 17 Random Effect SD         calm lag1_conventional 0.1351080 0.1074978 0.008974315 0.3566013 1.000528 7260.226  7255.846
#> 18 Random Effect SD conventional lag1_conventional 0.1897339 0.1571840 0.014196576 0.4806879 1.000518 6041.512  7858.500
#> 19 Random Effect SD     critical lag1_conventional 0.2396792 0.2018236 0.017498089 0.5958192 1.000120 6233.629  8242.834
#> 20 Random Effect SD   dependable lag1_conventional 0.2494178 0.2068788 0.018816474 0.6267372 1.000221 4653.808  6644.290
#> 21 Random Effect SD      anxious     lag1_critical 0.3101418 0.2717926 0.025833968 0.7256860 1.000891 5482.684  7204.193
#> 22 Random Effect SD         calm     lag1_critical 0.1664858 0.1323557 0.011989979 0.4384806 1.000182 8193.372  8237.181
#> 23 Random Effect SD conventional     lag1_critical 0.2696374 0.2321138 0.023865145 0.6409327 1.000172 6659.391  7381.395
#> 24 Random Effect SD     critical     lag1_critical 0.4338374 0.4035013 0.050755481 0.9331531 1.000605 4276.709  5323.753
#> 25 Random Effect SD   dependable     lag1_critical 0.2796911 0.2270858 0.021128453 0.7212001 1.000758 6098.635  8119.354
#> 26 Random Effect SD      anxious   lag1_dependable 0.1942395 0.1655554 0.016816816 0.4720359 1.000271 6012.751  7057.398
#> 27 Random Effect SD         calm   lag1_dependable 0.1319020 0.1083331 0.009981044 0.3357622 1.000104 6782.349  6758.308
#> 28 Random Effect SD conventional   lag1_dependable 0.1473785 0.1175822 0.010304257 0.3863459 1.000707 7119.038  6843.407
#> 29 Random Effect SD     critical   lag1_dependable 0.1564366 0.1250607 0.011192005 0.4120482 1.000184 7757.491  8173.588
#> 30 Random Effect SD   dependable   lag1_dependable 0.2361776 0.1949235 0.018555265 0.5932606 1.000534 4333.254  6605.437
```
