# ──────────────────────────────────────────────────────────────────────────────
# test-sim_bvarnet.R — simulation function tests
#
# Follows the bgms stochastic-robust testing philosophy:
#   • Range invariants (binary ∈ {0,1}, ordinal ∈ {1,...,C}, gaussian finite)
#   • Dimension consistency (N × T_obs × p layout)
#   • Seed reproducibility
#   • Parameter effect direction (coarse, wide tolerance)
#   • Burnin behaviour
#   • VAR stability enforcement
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 Return structure
# ═══════════════════════════════════════════════════════════════════════════════

test_that("sim_var returns correct structure for all families", {
  for (fam in c("bernoulli", "gaussian", "ordinal")) {
    result <- sim_var(N = 3, T_obs = 20, p = 2, K = 1,
                      family = fam, q = 1, seed = 1, burnin = 10)

    expect_type(result, "list")
    expect_named(result, c("data", "truth"))

    # data is a data.frame
    expect_s3_class(result$data, "data.frame")

    # truth is a list
    expect_type(result$truth, "list")

    # truth always contains these core elements
    expect_true(all(c("alpha", "Phi", "family", "N", "T_obs", "p", "K", "q",
                       "sd_alpha", "alpha_i", "Phi_i", "burnin") %in%
                      names(result$truth)))
  }
})


test_that("truth structure varies correctly by family", {
  gauss <- sim_var(N = 3, T_obs = 20, p = 2, family = "gaussian", seed = 1, burnin = 10)
  expect_true("sigma" %in% names(gauss$truth))
  expect_true(all(!is.na(gauss$truth$sigma)))
  # kappa entries are NULL for non-ordinal

  expect_true(all(vapply(gauss$truth$kappa, is.null, logical(1L))))
  expect_null(gauss$truth$C)

  ordinal <- sim_var(N = 3, T_obs = 20, p = 2, family = "ordinal", seed = 1, burnin = 10)
  expect_true("kappa" %in% names(ordinal$truth))
  expect_true("C" %in% names(ordinal$truth))
  # sigma is NA for ordinal nodes
  expect_true(all(is.na(ordinal$truth$sigma)))

  bern <- sim_var(N = 3, T_obs = 20, p = 2, family = "bernoulli", seed = 1, burnin = 10)
  expect_true(all(is.na(bern$truth$sigma)))
  expect_true(all(vapply(bern$truth$kappa, is.null, logical(1L))))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 Dimension consistency
# ═══════════════════════════════════════════════════════════════════════════════

test_that("output data frame has correct dimensions and columns", {
  N <- 5; T_obs <- 30; p <- 3; q <- 2
  result <- sim_var(N = N, T_obs = T_obs, p = p, q = q,
                    family = "bernoulli", seed = 10, burnin = 50)
  df <- result$data

  expect_equal(nrow(df), N * T_obs)
  expect_true("id" %in% names(df))
  expect_true("t"  %in% names(df))

  y_cols <- paste0("y_", 1:p)
  x_cols <- paste0("x_", 1:q)
  expect_true(all(y_cols %in% names(df)))
  expect_true(all(x_cols %in% names(df)))
})


test_that("truth matrices have correct dimensions", {
  N <- 4; p <- 3; K <- 2; q <- 2
  result <- sim_var(N = N, T_obs = 20, p = p, K = K, q = q,
                    family = "gaussian", seed = 20, burnin = 10)
  truth <- result$truth

  expect_equal(length(truth$alpha), p)
  expect_equal(dim(truth$Phi), c(p * K, p))
  expect_equal(dim(truth$gamma), c(q, p))
  expect_equal(length(truth$sigma), p)
  expect_equal(dim(truth$alpha_i), c(N, p))
  expect_equal(dim(truth$Phi_i), c(N, p * K, p))
})


test_that("ordinal kappa has correct structure", {
  p <- 3; C <- 5
  result <- sim_var(N = 3, T_obs = 20, p = p, family = "ordinal",
                    C = C, seed = 30, burnin = 10)
  kappa <- result$truth$kappa

  expect_type(kappa, "list")
  expect_length(kappa, p)

  for (node in seq_len(p)) {
    expect_length(kappa[[node]], C - 1)
    # cutpoints must be strictly ordered
    expect_true(!is.unsorted(kappa[[node]], strictly = TRUE))
  }
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 Range invariants
# ═══════════════════════════════════════════════════════════════════════════════

test_that("bernoulli responses are 0 or 1", {
  result <- sim_var(N = 10, T_obs = 40, p = 4, family = "bernoulli",
                    seed = 100, burnin = 50)
  y_cols <- paste0("y_", 1:4)
  for (col in y_cols) {
    vals <- result$data[[col]]
    expect_true(all(vals %in% c(0L, 1L)))
  }
})


test_that("ordinal responses are in {1, ..., C}", {
  C <- 5
  result <- sim_var(N = 10, T_obs = 40, p = 4, family = "ordinal",
                    C = C, seed = 101, burnin = 50)
  y_cols <- paste0("y_", 1:4)
  for (col in y_cols) {
    vals <- result$data[[col]]
    expect_true(all(vals %in% 1:C))
  }
})


test_that("gaussian responses are finite", {
  result <- sim_var(N = 10, T_obs = 40, p = 4, family = "gaussian",
                    seed = 102, burnin = 50)
  y_cols <- paste0("y_", 1:4)
  for (col in y_cols) {
    vals <- result$data[[col]]
    expect_true(all(is.finite(vals)))
  }
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 Seed reproducibility
# ═══════════════════════════════════════════════════════════════════════════════

test_that("same seed produces identical output", {
  r1 <- sim_var(N = 5, T_obs = 20, p = 2, family = "bernoulli",
                q = 1, seed = 777, burnin = 10)
  r2 <- sim_var(N = 5, T_obs = 20, p = 2, family = "bernoulli",
                q = 1, seed = 777, burnin = 10)

  expect_identical(r1$data, r2$data)
  expect_identical(r1$truth$Phi, r2$truth$Phi)
  expect_identical(r1$truth$alpha, r2$truth$alpha)
  expect_identical(r1$truth$alpha_i, r2$truth$alpha_i)
})


test_that("different seeds produce different output", {
  r1 <- sim_var(N = 5, T_obs = 20, p = 2, family = "bernoulli",
                q = 1, seed = 777, burnin = 10)
  r2 <- sim_var(N = 5, T_obs = 20, p = 2, family = "bernoulli",
                q = 1, seed = 888, burnin = 10)

  expect_false(identical(r1$data, r2$data))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §5 Burnin behaviour
# ═══════════════════════════════════════════════════════════════════════════════

test_that("burnin = 0 is allowed and produces T_obs rows per subject", {
  result <- sim_var(N = 3, T_obs = 20, p = 2, family = "bernoulli",
                    seed = 50, burnin = 0)
  expect_equal(nrow(result$data), 3 * 20)
  expect_equal(result$truth$burnin, 0L)
})


test_that("increasing burnin does not change output dimensions", {
  r1 <- sim_var(N = 3, T_obs = 20, p = 2, family = "gaussian",
                seed = 51, burnin = 10)
  r2 <- sim_var(N = 3, T_obs = 20, p = 2, family = "gaussian",
                seed = 51, burnin = 500)

  expect_equal(nrow(r1$data), nrow(r2$data))
  expect_equal(r1$truth$T_obs, r2$truth$T_obs)
})


test_that("burnin and non-burnin produce different data (same seed)", {
  r1 <- sim_var(N = 5, T_obs = 30, p = 2, family = "gaussian",
                seed = 52, burnin = 0)
  r2 <- sim_var(N = 5, T_obs = 30, p = 2, family = "gaussian",
                seed = 52, burnin = 500)

  # Same seed but different burnin → different recorded data
  expect_false(identical(r1$data, r2$data))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §6 VAR stability
# ═══════════════════════════════════════════════════════════════════════════════

test_that("auto-generated Phi is always stable", {
  for (trial in 1:5) {
    result <- sim_var(N = 3, T_obs = 20, p = 4, K = 2,
                      family = "gaussian", seed = trial, burnin = 10)
    Phi <- result$truth$Phi
    p <- result$truth$p; K <- result$truth$K

    companion <- bvarnet:::build_companion(Phi, p, K)
    max_ev <- max(abs(eigen(companion, only.values = TRUE)$values))
    expect_lt(max_ev, 1.0)
  }
})


test_that("unstable user-supplied Phi triggers warning and rescaling", {
  # Create an obviously unstable Phi (diagonal >> 1)
  unstable_Phi <- matrix(0, 3, 3)
  diag(unstable_Phi) <- 2.0

  expect_warning(
    result <- sim_var(N = 3, T_obs = 20, p = 3, K = 1,
                      family = "gaussian", Phi = unstable_Phi,
                      seed = 60, burnin = 10),
    "not VAR-stable"
  )

  # After rescaling, stability should hold
  Phi <- result$truth$Phi
  expect_true(bvarnet:::check_var_stability(Phi, 3, 1))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §7 No random effects mode (sd_alpha = 0)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("sd_alpha = 0 produces identical intercepts across subjects", {
  result <- sim_var(N = 10, T_obs = 20, p = 3, family = "gaussian",
                    sd_alpha = 0, seed = 70, burnin = 10)

  alpha_i <- result$truth$alpha_i
  # All rows should equal the population alpha
  for (i in seq_len(nrow(alpha_i))) {
    expect_equal(alpha_i[i, ], result$truth$alpha, tolerance = 1e-12)
  }
})


test_that("re_temporal = TRUE produces person-varying lag coefficients", {
  result <- sim_var(N = 5, T_obs = 20, p = 2, K = 1,
                    family = "gaussian", re_temporal = TRUE,
                    sd_phi = 0.3, seed = 71, burnin = 10)

  Phi_i <- result$truth$Phi_i
  # Check that at least two subjects differ in Phi
  expect_false(identical(Phi_i[1, , ], Phi_i[2, , ]))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §8 Covariate generation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("covariates are included when q > 0", {
  result <- sim_var(N = 3, T_obs = 20, p = 2, q = 3,
                    family = "bernoulli", seed = 80, burnin = 10)
  df <- result$data

  x_cols <- paste0("x_", 1:3)
  expect_true(all(x_cols %in% names(df)))

  # truth should have gamma

  expect_true(!is.null(result$truth$gamma))
  expect_equal(dim(result$truth$gamma), c(3, 2))
})


test_that("no covariate columns when q = 0", {
  result <- sim_var(N = 3, T_obs = 20, p = 2, q = 0,
                    family = "bernoulli", seed = 81, burnin = 10)
  df <- result$data

  expect_false(any(grepl("^x_", names(df))))
  expect_null(result$truth$gamma)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §9 K > 1 (higher-order VAR)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("K = 2 produces correct Phi dimensions and valid output", {
  result <- sim_var(N = 3, T_obs = 30, p = 3, K = 2,
                    family = "gaussian", seed = 90, burnin = 10)

  expect_equal(dim(result$truth$Phi), c(6, 3))  # p*K x p
  expect_equal(nrow(result$data), 3 * 30)
  expect_true(all(is.finite(result$data$y_1)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §10 User-supplied parameters
# ═══════════════════════════════════════════════════════════════════════════════

test_that("user-supplied alpha is stored in truth", {
  my_alpha <- c(0.5, -0.5)
  result <- sim_var(N = 3, T_obs = 20, p = 2, family = "gaussian",
                    alpha = my_alpha, seed = 95, burnin = 10)
  expect_equal(result$truth$alpha, my_alpha)
})


test_that("user-supplied Phi (stable) is preserved", {
  Phi <- matrix(c(0.3, 0.1, -0.1, 0.4), 2, 2)
  result <- sim_var(N = 3, T_obs = 20, p = 2, K = 1,
                    family = "gaussian", Phi = Phi, seed = 96, burnin = 10)
  expect_equal(result$truth$Phi, Phi)
})


test_that("user-supplied sigma is preserved in gaussian truth", {
  my_sigma <- c(1.0, 2.0, 0.5)
  result <- sim_var(N = 3, T_obs = 20, p = 3, family = "gaussian",
                    sigma = my_sigma, seed = 97, burnin = 10)
  expect_equal(result$truth$sigma, my_sigma)
})


test_that("user-supplied kappa is preserved in ordinal truth", {
  my_kappa <- list(c(-2, -1, 0, 1), c(-1.5, -0.5, 0.5, 1.5))
  result <- sim_var(N = 3, T_obs = 20, p = 2, family = "ordinal",
                    C = 5, kappa = my_kappa, seed = 98, burnin = 10)
  expect_equal(result$truth$kappa, my_kappa)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §11 Internal helpers
# ═══════════════════════════════════════════════════════════════════════════════

test_that("build_companion returns correct dimensions", {
  p <- 3; K <- 2
  Phi <- matrix(runif(p * K * p, -0.3, 0.3), p * K, p)
  comp <- bvarnet:::build_companion(Phi, p, K)
  expect_equal(dim(comp), c(p * K, p * K))
})


test_that("build_companion for K=1 equals t(Phi)", {
  Phi <- matrix(c(0.3, 0.1, 0.05, 0.4), 2, 2)
  comp <- bvarnet:::build_companion(Phi, 2, 1)
  expect_equal(comp, t(Phi))
})


test_that("rescale_to_stable makes unstable Phi stable", {
  Phi <- matrix(0, 3, 3)
  diag(Phi) <- 1.5  # unstable
  Phi_stable <- bvarnet:::rescale_to_stable(Phi, 3, 1)
  expect_true(bvarnet:::check_var_stability(Phi_stable, 3, 1))
})


test_that("generate_default_kappa returns ordered cutpoints", {
  for (C in c(2, 3, 5, 10)) {
    kappa <- bvarnet:::generate_default_kappa(C)
    expect_length(kappa, C - 1)
    expect_true(!is.unsorted(kappa, strictly = TRUE))
  }
})


test_that("generate_response_binary returns 0/1 integers", {
  set.seed(200)
  eta <- rnorm(100)
  y <- bvarnet:::generate_response_binary(eta)
  expect_true(all(y %in% c(0L, 1L)))
  expect_length(y, 100)
})


test_that("generate_response_gaussian returns finite values", {
  set.seed(201)
  eta <- rnorm(50)
  sigma <- rep(1, 50)
  y <- bvarnet:::generate_response_gaussian(eta, sigma)
  expect_true(all(is.finite(y)))
  expect_length(y, 50)
})


test_that("generate_response_ordinal returns values in {1, ..., C}", {
  set.seed(202)
  p <- 5; C <- 4
  eta <- rnorm(p)
  kappa <- replicate(p, sort(rnorm(C - 1)), simplify = FALSE)
  y <- bvarnet:::generate_response_ordinal(eta, kappa, C, p)
  expect_true(all(y %in% 1:C))
  expect_length(y, p)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §12 Stochastic property tests (wide tolerance)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("bernoulli: large positive intercept → high marginal P(Y=1)", {
  # With very large intercepts, most responses should be 1
  result <- sim_var(N = 20, T_obs = 100, p = 2, family = "bernoulli",
                    alpha = c(5, 5), sd_alpha = 0, burnin = 100, seed = 300)
  y_cols <- paste0("y_", 1:2)
  for (col in y_cols) {
    prop <- mean(result$data[[col]])
    expect_gt(prop, 0.8)
  }
})


test_that("bernoulli: large negative intercept → low marginal P(Y=1)", {
  result <- sim_var(N = 20, T_obs = 100, p = 2, family = "bernoulli",
                    alpha = c(-5, -5), sd_alpha = 0, burnin = 100, seed = 301)
  y_cols <- paste0("y_", 1:2)
  for (col in y_cols) {
    prop <- mean(result$data[[col]])
    expect_lt(prop, 0.2)
  }
})


test_that("gaussian: zero Phi → means approx equal to intercept", {
  alpha_true <- c(3, -2, 0.5)
  Phi_zero <- matrix(0, 3, 3)
  result <- sim_var(N = 50, T_obs = 100, p = 3, K = 1,
                    family = "gaussian", alpha = alpha_true,
                    Phi = Phi_zero, sd_alpha = 0, burnin = 100, seed = 302)

  for (j in 1:3) {
    col <- paste0("y_", j)
    obs_mean <- mean(result$data[[col]])
    # Wide tolerance: within 0.5 of true intercept
    expect_equal(obs_mean, alpha_true[j], tolerance = 0.5)
  }
})


test_that("ordinal: extreme kappa shifts marginal distribution", {
  # Very negative kappa → most mass on higher categories
  kappa_low <- list(c(-10, -9, -8, -7), c(-10, -9, -8, -7))
  result_low <- sim_var(N = 20, T_obs = 100, p = 2, family = "ordinal",
                        C = 5, kappa = kappa_low, sd_alpha = 0,
                        burnin = 100, seed = 303)

  # Very positive kappa → most mass on lower categories
  kappa_high <- list(c(7, 8, 9, 10), c(7, 8, 9, 10))
  result_high <- sim_var(N = 20, T_obs = 100, p = 2, family = "ordinal",
                         C = 5, kappa = kappa_high, sd_alpha = 0,
                         burnin = 100, seed = 304)

  mean_low  <- mean(result_low$data$y_1)
  mean_high <- mean(result_high$data$y_1)

  # Low kappa → higher category mean, High kappa → lower category mean
  expect_gt(mean_low, mean_high)
})
