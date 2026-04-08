# bvarnet

This package implements methods to analyse (multilevel) Vector
Autoregression commonly used for the analysis of intensive longitudinal
datasets like EMA data.

In the following vignette we will introduce the basic functionality of
the package. For more elaborated usecases, we refer to the additional
vignettes

- [`vignette("Hypothesis-Testing")`](https://flo1met.github.io/bvarnet/articles/Hypothesis-Testing.md)
- [`vignette("Mixed-Model")`](https://flo1met.github.io/bvarnet/articles/Mixed-Model.md)
- [`vignette("Missing-Data")`](https://flo1met.github.io/bvarnet/articles/Missing-Data.md)
- `vignette("Random-Effect")`
- [`vignette("Prediction")`](https://flo1met.github.io/bvarnet/articles/Prediction.md)

First, we have to load the package:

``` r
library(bvarnet)

# subject to be removed again...
library(openesm)
#> Welcome to openesm!
#> 
#> Find documentation and usage examples at https://openesmdata.org
#> 
#> Get started:
#>   * list_datasets() to browse available datasets
#>   * get_dataset('dataset_id') to download a specific dataset
#>   * ?get_dataset or visit the website for detailed guides
#> 
#> Attaching package: 'openesm'
#> The following object is masked from 'package:utils':
#> 
#>     cite
library(qgraph)
```

There is some missing data in the dataset. The models default options
handle this by themselves. For a further elaboration on this, you can
read
[`vignette("Missing-Data")`](https://flo1met.github.io/bvarnet/articles/Missing-Data.md).

``` r
data <- openesm::get_dataset("0004")
#> ℹ Using cached dataset index (less than 24 hours old).
#> ✔ Loading dataset "0004" version "1.0.0"
#> 
#> ── openESM Dataset: "0004" ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#> • Dataset version: "1.0.0"
#> • Metadata version: "1.0.1"
#> • Authors: Wang et al. (2014)
#> • Paper DOI: https://doi.org/10.1371/journal.pone.0266516
#> • License: CC BY-NC 4.0
#> • Data: A tibble with 49 participants and 64 maximum time points per participant
#> ℹ Use `cite(dataset)` for citation information.
#> ℹ Use `notes(dataset)` for additional information about the dataset.
#> ℹ Please ensure you follow the license terms for this dataset.

df <- data$data

df$anxious <- floor(df$anxious)
df$calm <- floor(df$calm)
df$conventional <- floor(df$conventional)
df$critical <- floor(df$critical)
df$dependable <- floor(df$dependable)
```

``` r
fit <- bvar(
  id_col = "id",
  time_col = "day",
  y_cols = c("anxious", "calm", "conventional", "critical", "dependable"),
  x_cols = NULL,
  re_cols = NULL,
  re_temporal = FALSE,
  K = 1,
  data = df,
  family = c("ordinal"),
  priors = set_priors(),
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  seed = 1337
)
#> Error:
#> ! Stan model file "" not found.
```

Most models will give a warning when initiating sampling, that looks
like this:

``` text
Chain 1 Informational Message: The current Metropolis proposal is about to be rejected because of the following issue:
Chain 1 Exception: categorical_logit_glm_lpmf: Intercept[3] is -inf, but must be finite! (in '/var/folders/n5/38kfmkv55hq8bnd03344d2b40000gn/T/RtmpkjgDq1/model-a60f79e9a37c.stan', line 153, column 8 to column 94)
Chain 1 If this warning occurs sporadically, such as for highly constrained variable types like covariance matrices, then the sampler is fine,
Chain 1 but if this warning occurs often then your model may be either severely ill-conditioned or misspecified.
Chain 1 
```

We can print an overview of the model we just ran using the `print(fit)`
functions:

``` r
print(fit)
#> BVAR Network fit
#> ======================================== 
#> Family:      ordinal
#> Outcomes (p): 5 
#> Lags (K):     1 
#> Fixed eff.:   0 
#> Observations: 147 
#> Rhat max:    1.001
#> Divergences: 2  WARNING: check model/priors.
#> Priors:      beta ~ Normal(0, 1), phi ~ Normal(0, 0.5), sd_u ~ Half-Normal(0, 1), kappa ~ Normal(0, 2)
#> Total time:  17.8 sec
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
#> Family: ordinal | p=5 | K=1 | n=147
#> Rhat max: 1.001 | Divergences: 2
#>   WARNING: divergent transitions detected — check model/priors.
#> 
#> --- Autoregressive ---
#>  predictor         outcome      mean  median q5     q95   rhat ess_bulk ess_tail
#>  lag1_anxious      anxious      0.203 0.203  -0.054 0.460 1    15617.37 12599.58
#>  lag1_calm         calm         0.176 0.174  -0.058 0.415 1    13025.39 11190.21
#>  lag1_conventional conventional 0.558 0.554   0.280 0.843 1    14346.16 11205.97
#>  lag1_critical     critical     0.278 0.277   0.075 0.492 1    18103.80 11283.74
#>  lag1_dependable   dependable   0.476 0.473   0.247 0.714 1    15253.27 12152.74
#> 
#> 
#> --- Cross-lagged ---
#>  predictor         outcome      mean   median q5     q95    rhat ess_bulk ess_tail
#>  lag1_calm         anxious      -0.308 -0.307 -0.556 -0.069 1    12276.01 12256.25
#>  lag1_conventional anxious       0.113  0.112 -0.145  0.375 1    15288.24 12209.56
#>  lag1_critical     anxious       0.065  0.065 -0.142  0.273 1    21335.42 11693.19
#>  lag1_dependable   anxious       0.054  0.052 -0.158  0.267 1    14645.31 12404.76
#>  lag1_anxious      calm         -0.085 -0.085 -0.338  0.168 1    15499.73 13047.83
#>  lag1_conventional calm         -0.162 -0.162 -0.423  0.096 1    14847.92 12191.17
#>  lag1_critical     calm         -0.077 -0.077 -0.286  0.131 1    20946.54 12223.59
#>  lag1_dependable   calm          0.116  0.115 -0.096  0.328 1    15982.24 11988.57
#>  lag1_anxious      conventional -0.184 -0.184 -0.462  0.095 1    14348.25 12403.62
#>  lag1_calm         conventional -0.105 -0.105 -0.349  0.134 1    12989.52 12059.56
#> 
#> ... 10 more rows. Use extract_temporal(fit, effect = "cl") for full output.
#> 
#> --- Threshold ---
#>  predictor               outcome mean   median q5     q95    rhat  ess_bulk ess_tail 
#>  kappa(anxious, c1)      —       -1.033 -1.029 -1.470 -0.613 1.000 11123.68  9566.673
#>  kappa(calm, c1)         —       -1.300 -1.266 -1.869 -0.848 1.000 14564.11 10916.416
#>  kappa(conventional, c1) —       -1.033 -1.012 -1.471 -0.674 1.000 14373.14 11229.523
#>  kappa(critical, c1)     —        0.086  0.095 -0.239  0.379 1.000 12541.34 10745.366
#>  kappa(dependable, c1)   —       -1.726 -1.680 -2.521 -1.081 1.000 10748.70 10274.783
#>  kappa(anxious, c2)      —        0.459  0.464  0.139  0.765 1.000 13762.57 12456.461
#>  kappa(calm, c2)         —       -0.884 -0.875 -1.244 -0.552 1.000 18423.96 13368.980
#>  kappa(conventional, c2) —       -0.707 -0.710 -1.008 -0.397 1.000 15460.49 10840.830
#>  kappa(critical, c2)     —        0.548  0.547  0.298  0.806 1.000 16603.34 12957.191
#>  kappa(dependable, c2)   —       -0.859 -0.860 -1.263 -0.456 1.001 10048.09  6937.794
#> 
#> ... 10 more rows. Use extract_param(fit) for full output.
#> 
#> ==================================================
#> Use extract_param() for the full parameter table.
#> Use extract_network_matrix() for the temporal network matrix.
```

Here we can see, that we can not see all category threshold parameters
($\kappa$). To inspect them completely we have to extract them using
[`extract_param()`](https://flo1met.github.io/bvarnet/reference/extract_param.md):

``` r
params <- extract_param(fit)
params[params$type == "Threshold",] # add seperate extractor function for kappa, sigma..
#>         type               predictor outcome        mean      median         q5        q95      rhat ess_bulk  ess_tail
#> 26 Threshold      kappa(anxious, c1)       — -1.03285109 -1.02852445 -1.4703235 -0.6134608 1.0001894 11123.68  9566.673
#> 27 Threshold         kappa(calm, c1)       — -1.30014533 -1.26565785 -1.8692502 -0.8476075 1.0001328 14564.11 10916.416
#> 28 Threshold kappa(conventional, c1)       — -1.03336729 -1.01205985 -1.4711615 -0.6738048 1.0000887 14373.13 11229.523
#> 29 Threshold     kappa(critical, c1)       —  0.08567509  0.09505895 -0.2386167  0.3793493 1.0003189 12541.34 10745.366
#> 30 Threshold   kappa(dependable, c1)       — -1.72568614 -1.68045880 -2.5214071 -1.0805193 1.0003006 10748.70 10274.783
#> 31 Threshold      kappa(anxious, c2)       —  0.45936818  0.46442090  0.1385890  0.7650862 1.0000192 13762.57 12456.461
#> 32 Threshold         kappa(calm, c2)       — -0.88438434 -0.87529079 -1.2444187 -0.5517570 1.0001721 18423.96 13368.980
#> 33 Threshold kappa(conventional, c2)       — -0.70673963 -0.71005004 -1.0081459 -0.3967633 0.9999131 15460.49 10840.830
#> 34 Threshold     kappa(critical, c2)       —  0.54833491  0.54694036  0.2981471  0.8061212 1.0002716 16603.34 12957.191
#> 35 Threshold   kappa(dependable, c2)       — -0.85949747 -0.85968742 -1.2632601 -0.4558807 1.0010041 10048.09  6937.794
#> 36 Threshold      kappa(anxious, c3)       —  0.95189554  0.93533212  0.5917419  1.3648925 0.9999869 14502.71 11369.641
#> 37 Threshold         kappa(calm, c3)       — -0.36802642 -0.37436518 -0.6744710 -0.0440537 1.0001101 13760.73 12421.296
#> 38 Threshold kappa(conventional, c3)       —  0.62996635  0.62775770  0.2608396  1.0047809 0.9999474 13535.62 10171.390
#> 39 Threshold     kappa(critical, c3)       —  0.76699860  0.75749318  0.5062937  1.0600305 0.9999889 17305.51 12710.418
#> 40 Threshold   kappa(dependable, c3)       —  0.12516613  0.12449102 -0.2202632  0.4724522 1.0004251 11098.90  7811.925
#> 41 Threshold      kappa(anxious, c4)       —  1.92440208  1.88030470  1.1873595  2.8280811 1.0000381 14967.02 10543.057
#> 42 Threshold         kappa(calm, c4)       —  1.35545685  1.34728995  0.9012824  1.8378872 0.9999767 18072.38 11041.277
#> 43 Threshold kappa(conventional, c4)       —  2.57589420  2.55145715  1.8081027  3.4241004 1.0004303 19188.28 10144.966
#> 44 Threshold     kappa(critical, c4)       —  1.15580404  1.10903125  0.7258669  1.7561748 1.0001773 20899.16 12835.077
#> 45 Threshold   kappa(dependable, c4)       —  1.15106984  1.14520455  0.6797434  1.6414832 1.0000701 14456.67 10624.680
```

If we are interested in the temporal network structure, we can inspect
this using either the `extract_temporal(fit)`, or the
`extract_network_matrix(fit)` functions. Here we will use the
`extract_network_matrix(fit)` to plot the network:

``` r
nw_mat <- extract_network_matrix(fit)
qgraph::qgraph(nw_mat)
```

![plot of chunk unnamed-chunk-31](figure/unnamed-chunk-31-1.png)

plot of chunk unnamed-chunk-31
