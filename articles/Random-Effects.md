# Random-Effects

## Setup

``` r
library(bvarnet)
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
#> Total time:  42.9 sec
#> ========================================
summary(fit)
#> BVAR Network Summary
#> ================================================== 
#> Family: ordinal | p=5 | K=1 | n=147
#> Rhat max: 1.002 | Divergences: 0
#> 
#> --- Autoregressive ---
#>  predictor         outcome      mean   median q5     q95   rhat  ess_bulk  ess_tail 
#>  lag1_anxious      anxious      -0.150 -0.147 -0.553 0.242 1.000 12127.960 11495.325
#>  lag1_calm         calm         -0.164 -0.165 -0.518 0.187 1.001  9331.043 10426.126
#>  lag1_conventional conventional  0.097  0.096 -0.326 0.518 1.000  9984.664 10959.196
#>  lag1_critical     critical      0.020  0.043 -0.475 0.453 1.001  6178.876  9401.899
#>  lag1_dependable   dependable    0.434  0.434  0.075 0.793 1.000 11031.359 11762.321
#> 
#> 
#> --- Cross-lagged ---
#>  predictor         outcome      mean   median q5     q95   rhat  ess_bulk  ess_tail 
#>  lag1_calm         anxious      -0.060 -0.058 -0.421 0.300 1.000  9860.096 10454.267
#>  lag1_conventional anxious      -0.147 -0.143 -0.580 0.265 1.000 10649.275 11114.092
#>  lag1_critical     anxious       0.429  0.432 -0.003 0.842 1.000  7860.406  9727.691
#>  lag1_dependable   anxious      -0.029 -0.028 -0.380 0.316 1.000  8579.806 10701.049
#>  lag1_anxious      calm          0.009  0.007 -0.358 0.385 1.000 10943.302 10812.556
#>  lag1_conventional calm          0.114  0.113 -0.269 0.504 1.000  9787.858  9782.278
#>  lag1_critical     calm         -0.123 -0.125 -0.461 0.215 1.001 11641.220 11337.815
#>  lag1_dependable   calm          0.171  0.170 -0.147 0.497 1.000 10875.083 11537.899
#>  lag1_anxious      conventional -0.176 -0.178 -0.590 0.243 1.000 11803.823 11360.979
#>  lag1_calm         conventional  0.004  0.004 -0.377 0.380 1.001  8710.923  9393.587
#> 
#> ... 10 more rows. Use extract_temporal(fit, effect = "cl") for full output.
#> 
#> --- Random Effect SD ---
#>  predictor    outcome      mean  median q5    q95   rhat  ess_bulk ess_tail
#>  anxious      Intercept    1.465 1.444  0.978 2.028 1.002 5265.490 7743.781
#>  calm         Intercept    1.178 1.155  0.780 1.658 1.001 5785.705 8819.844
#>  conventional Intercept    1.127 1.106  0.737 1.590 1.001 6009.031 7327.977
#>  critical     Intercept    1.363 1.336  0.858 1.974 1.001 4875.025 8113.773
#>  dependable   Intercept    0.966 0.944  0.583 1.422 1.000 5336.862 8154.476
#>  anxious      lag1_anxious 0.202 0.169  0.016 0.503 1.001 5809.229 6993.924
#>  calm         lag1_anxious 0.161 0.129  0.012 0.419 1.001 6742.384 6608.842
#>  conventional lag1_anxious 0.254 0.212  0.019 0.633 1.000 5586.614 6996.407
#>  critical     lag1_anxious 0.299 0.256  0.023 0.726 1.000 4817.791 6684.910
#>  dependable   lag1_anxious 0.299 0.240  0.023 0.768 1.002 4628.541 7090.355
#> 
#> ... 20 more rows. Use extract_random_effects(fit) for full output.
#> 
#> --- Threshold ---
#>  predictor               outcome mean   median q5     q95    rhat  ess_bulk  ess_tail 
#>  kappa(anxious, c1)      —       -2.099 -2.074 -2.909 -1.371 1.000  5356.457  8025.435
#>  kappa(calm, c1)         —       -2.226 -2.189 -3.144 -1.437 1.001  7422.138 10360.114
#>  kappa(conventional, c1) —       -2.035 -2.006 -2.892 -1.284 1.000  6615.821  9820.399
#>  kappa(critical, c1)     —       -0.265 -0.259 -0.890  0.333 1.000  5610.687  7665.892
#>  kappa(dependable, c1)   —       -2.505 -2.469 -3.507 -1.609 1.000  8289.672  9571.903
#>  kappa(anxious, c2)      —        0.446  0.444 -0.186  1.079 1.001  6012.700  9350.799
#>  kappa(calm, c2)         —       -1.506 -1.497 -2.166 -0.891 1.000  8043.093  9262.475
#>  kappa(conventional, c2) —       -1.258 -1.252 -1.867 -0.678 1.000  8082.262 10988.755
#>  kappa(critical, c2)     —        1.169  1.157  0.561  1.821 1.000  6447.810  9636.228
#>  kappa(dependable, c2)   —       -1.241 -1.230 -1.852 -0.665 1.000 11358.657 12787.251
#> 
#> ... 10 more rows. Use extract_param(fit, type = "Threshold") for full output.
#> 
#> ==================================================
#> Use extract_param() or extract_param(fit, type = "...") for the full parameter table.
#> Use extract_network_matrix() for the temporal network matrix.
```

## Extracting Random Effects

Additionally to the `extract_*` functions that we already described in
`Vignette(bvarnet)`, we can use the
[`extract_random_effects()`](https://flo1met.github.io/bvarnet/reference/extract_random_effects.md)
function to only extract the random effects:

``` r
re <- extract_random_effects(fit)
re
#>                type    predictor           outcome      mean    median          q5       q95      rhat ess_bulk ess_tail
#> 1  Random Effect SD      anxious         Intercept 1.4652817 1.4439362 0.977649331 2.0282597 1.0015907 5265.490 7743.781
#> 2  Random Effect SD         calm         Intercept 1.1783653 1.1549465 0.779576993 1.6583403 1.0009675 5785.705 8819.844
#> 3  Random Effect SD conventional         Intercept 1.1268582 1.1063721 0.736724058 1.5901623 1.0006901 6009.031 7327.977
#> 4  Random Effect SD     critical         Intercept 1.3634879 1.3362195 0.857573136 1.9744304 1.0006644 4875.025 8113.773
#> 5  Random Effect SD   dependable         Intercept 0.9660984 0.9443749 0.582831873 1.4222424 1.0004303 5336.862 8154.476
#> 6  Random Effect SD      anxious      lag1_anxious 0.2020052 0.1685053 0.016233745 0.5032068 1.0006659 5809.229 6993.924
#> 7  Random Effect SD         calm      lag1_anxious 0.1609254 0.1289984 0.012195480 0.4193622 1.0011414 6742.384 6608.842
#> 8  Random Effect SD conventional      lag1_anxious 0.2538834 0.2116468 0.019302953 0.6331248 1.0003903 5586.614 6996.407
#> 9  Random Effect SD     critical      lag1_anxious 0.2986767 0.2556347 0.022967897 0.7264043 1.0003753 4817.791 6684.910
#> 10 Random Effect SD   dependable      lag1_anxious 0.2987636 0.2399345 0.023133795 0.7682914 1.0021471 4628.541 7090.355
#> 11 Random Effect SD      anxious         lag1_calm 0.1394255 0.1112678 0.010881843 0.3615685 1.0002288 6431.470 6540.138
#> 12 Random Effect SD         calm         lag1_calm 0.1711120 0.1428554 0.013783058 0.4292299 1.0004208 4667.061 6146.463
#> 13 Random Effect SD conventional         lag1_calm 0.1648274 0.1344180 0.013319257 0.4176716 1.0006599 5273.898 7165.568
#> 14 Random Effect SD     critical         lag1_calm 0.1421269 0.1124673 0.010267661 0.3735351 1.0001253 7205.369 7004.134
#> 15 Random Effect SD   dependable         lag1_calm 0.1983639 0.1632243 0.015700375 0.5037090 1.0003510 4551.673 5418.571
#> 16 Random Effect SD      anxious lag1_conventional 0.1951214 0.1662170 0.016230538 0.4777358 1.0002119 5611.433 7220.911
#> 17 Random Effect SD         calm lag1_conventional 0.1382549 0.1102174 0.010239865 0.3621581 1.0015230 6263.953 6290.534
#> 18 Random Effect SD conventional lag1_conventional 0.1890291 0.1529822 0.014335658 0.4863923 1.0003510 5346.015 6306.486
#> 19 Random Effect SD     critical lag1_conventional 0.2409120 0.2021247 0.019716436 0.5921124 1.0004442 4726.675 5829.030
#> 20 Random Effect SD   dependable lag1_conventional 0.2495825 0.2088717 0.020474920 0.6229940 1.0006090 4476.032 7940.819
#> 21 Random Effect SD      anxious     lag1_critical 0.3066851 0.2681832 0.027216169 0.7317527 1.0005179 4969.800 6305.771
#> 22 Random Effect SD         calm     lag1_critical 0.1681908 0.1354386 0.012202660 0.4435305 1.0004330 7164.754 7388.659
#> 23 Random Effect SD conventional     lag1_critical 0.2692275 0.2309014 0.020660210 0.6603602 1.0001849 5570.973 5814.862
#> 24 Random Effect SD     critical     lag1_critical 0.4391420 0.4077634 0.053331996 0.9286383 1.0005797 3702.812 3961.537
#> 25 Random Effect SD   dependable     lag1_critical 0.2822168 0.2289254 0.021953689 0.7283563 0.9999966 5497.042 7410.768
#> 26 Random Effect SD      anxious   lag1_dependable 0.1979930 0.1712029 0.017472326 0.4800424 1.0009843 4905.884 7042.269
#> 27 Random Effect SD         calm   lag1_dependable 0.1325212 0.1072447 0.009313681 0.3400509 1.0007978 5465.635 5488.219
#> 28 Random Effect SD conventional   lag1_dependable 0.1451051 0.1148876 0.010273553 0.3802055 1.0003109 6261.696 6855.308
#> 29 Random Effect SD     critical   lag1_dependable 0.1577871 0.1263587 0.011028128 0.4126913 1.0001326 6949.019 6812.282
#> 30 Random Effect SD   dependable   lag1_dependable 0.2270672 0.1850662 0.017049599 0.5886999 1.0005338 4024.393 6203.294
```
