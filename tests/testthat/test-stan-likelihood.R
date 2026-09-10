# ──────────────────────────────────────────────────────────────────────────────
# test-stan-likelihood.R — do the .stan files compute the documented density?
#
# Each test optimises the compiled Stan model and an independent R
# implementation of the same log posterior (helper-oracles.R) on identical
# data, and asserts the two modes agree. Because `optimize(jacobian = FALSE)`
# maximises the constrained-space log density, it is directly comparable to a
# plain `optim()` over the natural parameters, and agreement is exact to
# optimiser tolerance rather than approximate.
#
# These are deterministic and take well under a second each once compiled — no
# MCMC, no sampling noise, no tolerance guesswork.
#
# What this class of test catches: sign flips, off-by-one category coding, the
# lagged outcome entering on the wrong scale, transposed phi indexing, priors
# applied to the wrong parameter block, and parameter declarations that
# constrain the space more than the model does (see the ordinal test).
# ──────────────────────────────────────────────────────────────────────────────

expect_stan_matches_oracle <- function(family, model_name, tolerance = 1e-3) {
  skip_if_no_stan_compile(model_name)

  sd_ <- make_oracle_standata(family)
  # These tests fit a fixed-effects model; the oracle assumes n_re == 0.
  expect_identical(sd_$n_re, 0L)

  standata <- sd_[intersect(names(sd_), .stan_data_names(family))]

  stan_est   <- stan_flatten(stan_map(model_name, standata), family, sd_)
  oracle_est <- oracle_flatten(oracle_map(family, sd_), family, sd_)

  expect_equal(stan_est, oracle_est, tolerance = tolerance)
  invisible(list(stan = stan_est, oracle = oracle_est))
}

#' Which standata elements the Stan program actually declares
#' (to_stan_data() also carries R-side prediction metadata)
.stan_data_names <- function(family) {
  base <- c("p", "J", "K", "n_obs", "n_fe", "n_re", "id", "Y", "X", "B", "Z",
            "prior_beta_fam", "beta_loc", "beta_scale", "beta_df",
            "prior_phi_fam", "phi_loc", "phi_scale", "phi_df",
            "prior_sd_fam", "sd_loc", "sd_scale", "sd_df")
  if (family %in% c("bernoulli", "gaussian"))
    base <- c(base, "prior_intercept_fam", "intercept_loc",
              "intercept_scale", "intercept_df")
  if (family == "gaussian")
    base <- c(base, "prior_sigma_fam", "sigma_loc", "sigma_scale", "sigma_df")
  if (family == "ordinal")
    base <- c(base, "C", "prior_kappa_fam", "kappa_loc",
              "kappa_scale", "kappa_df")
  base
}


test_that("model_binary computes the bernoulli-logit log posterior", {
  expect_stan_matches_oracle("bernoulli", "model_binary")
})


test_that("model_gaussian computes the normal-identity log posterior", {
  expect_stan_matches_oracle("gaussian", "model_gaussian")
})


test_that("model_ordinal computes the adjacent-category log posterior", {
  expect_stan_matches_oracle("ordinal", "model_ordinal")
})


test_that("adjacent-category thresholds are not order-constrained", {
  # Data are simulated with kappa = (1.6, 0.2, 1.1) and (0.9, -0.4, 1.2), both
  # non-monotone, so a constrained model cannot reach the true mode.
  skip_if_no_stan_compile("model_ordinal")

  res <- expect_stan_matches_oracle("ordinal", "model_ordinal")

  kappa <- matrix(res$stan[grep("^kappa", names(res$stan))],
                  nrow = 2L, byrow = FALSE)

  # 1. the fitted thresholds are genuinely non-monotone -- if this fails the
  #    fixture has stopped exercising the case the test exists for
  expect_true(any(apply(kappa, 1, is.unsorted, strictly = TRUE)),
              info = "fixture no longer produces non-monotone thresholds")

  # 2. no node has adjacent thresholds pinned together, the signature of an
  #    active ordering constraint
  gaps <- abs(as.vector(t(apply(kappa, 1, diff))))
  expect_true(all(gaps > 1e-4),
              info = "adjacent thresholds are tied -- kappa looks order-constrained")
})


test_that("sim_var round-trips non-monotone kappa", {
  # The simulator must not re-impose the constraint the Stan model dropped.
  kappa <- list(c(1.6, 0.2, 1.1), c(0.9, -0.4, 1.2))
  sim <- sim_var(N = 4, T_obs = 20, p = 2, K = 1, family = "ordinal",
                 C = 4, kappa = kappa, seed = 5, burnin = 0)
  expect_equal(sim$truth$kappa, kappa)
})
