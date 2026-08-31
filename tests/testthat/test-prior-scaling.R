# ──────────────────────────────────────────────────────────────────────────────
# test-prior-scaling.R — the prior a fit records must be the prior Stan used.
#
# Gaussian defaults are widened by the outcome SD before they reach Stan. If the
# object stores the unscaled prior, savage_dickey() divides by the wrong density
# and every default-prior BF on a Gaussian intercept/fixed effect is biased by
# roughly sd(y). These tests pin the two halves of that contract: the builders
# report what they sent, and savage_dickey() looks the prior up per node.
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 The builders report their effective priors
# ═══════════════════════════════════════════════════════════════════════════════

test_that("to_stan_data records the scales it sent to Stan (gaussian)", {
  df <- make_test_df(N = 4, T_obs = 15, p = 2, family = "gaussian")
  df$y_1 <- df$y_1 * 10   # force s_y far from 1

  sd_out <- to_stan_data(df, "gaussian", "id", "t", paste0("y_", 1:2),
                         character(0), K = 1, priors = set_priors())
  eff <- attr(sd_out, "priors_effective")

  expect_s3_class(eff, "bvarnet_priors")
  expect_equal(eff$beta$scale,      sd_out$beta_scale)
  expect_equal(eff$intercept$scale, sd_out$intercept_scale)
  expect_equal(eff$sigma$scale,     sd_out$sigma_scale)

  # And the scaling actually happened, so the test can distinguish the bug.
  expect_gt(sd_out$beta_scale, 1)
})

test_that("to_stan_data leaves non-gaussian and user-supplied priors unscaled", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "bernoulli")
  sd_b <- to_stan_data(df, "bernoulli", "id", "t", paste0("y_", 1:2),
                       character(0), K = 1, priors = set_priors())
  expect_equal(attr(sd_b, "priors_effective")$beta$scale, sd_b$beta_scale)
  expect_equal(sd_b$beta_scale, 1)

  dg <- make_test_df(N = 4, T_obs = 15, p = 2, family = "gaussian")
  dg$y_1 <- dg$y_1 * 10
  user <- set_priors(beta = prior("normal", 0, 3))
  sd_g <- to_stan_data(dg, "gaussian", "id", "t", paste0("y_", 1:2),
                       character(0), K = 1, priors = user)
  # User asked for 3; a user prior is taken at face value and never rescaled.
  expect_equal(sd_g$beta_scale, 3)
  expect_equal(attr(sd_g, "priors_effective")$beta$scale, 3)
  # The default intercept still scales.
  expect_gt(sd_g$intercept_scale, 1)
})

test_that(".to_stan_data_node records per-node scales", {
  set.seed(99)
  df <- data.frame(
    id  = rep(1:3, each = 15),
    t   = rep(1:15, 3),
    y_1 = rnorm(45, sd = 5),
    y_2 = rnorm(45, sd = 0.1)
  )
  shared <- .to_stan_data_shared(df, "id", "t", paste0("y_", 1:2),
                                 character(0), K = 1)
  priors <- set_priors()

  n1 <- .to_stan_data_node(shared, 1L, "gaussian", priors)
  n2 <- .to_stan_data_node(shared, 2L, "gaussian", priors)

  expect_equal(attr(n1, "priors_effective")$beta$scale, n1$beta_scale)
  expect_equal(attr(n2, "priors_effective")$beta$scale, n2$beta_scale)
  # y_1 is ~50x more variable than y_2, so the nodes must not share a scale.
  expect_gt(n1$beta_scale, n2$beta_scale)
})

test_that("effective priors are carried as an attribute, never as Stan data", {
  # bvar() filters standata by name, but .bvar_nodewise() passes its node list
  # to $sample() unfiltered — a named element would reach Stan's JSON.
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "gaussian")
  sd_out <- to_stan_data(df, "gaussian", "id", "t", paste0("y_", 1:2),
                         character(0), K = 1, priors = set_priors())
  expect_false("priors_effective" %in% names(sd_out))
  expect_false(is.null(attr(sd_out, "priors_effective")))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 savage_dickey() divides by the prior Stan actually used
# ═══════════════════════════════════════════════════════════════════════════════

test_that("savage_dickey uses the effective, not the requested, prior scale", {
  obj <- make_mock_bvarnet(family = "gaussian", n_iter = 1000L, n_chains = 1L)

  # As if Stan had been given a beta scale of 10 while the user asked for 1.
  eff <- set_priors()
  eff$beta$scale <- 10
  obj$priors_effective <- rep(list(eff), 2L)

  res <- savage_dickey(obj, "beta[2,1]", method = "logspline")
  expect_equal(res$prior_density, dnorm(0, 0, 10))
  expect_false(isTRUE(all.equal(res$prior_density, dnorm(0, 0, 1))))
})

test_that("savage_dickey resolves the effective prior per node", {
  obj <- make_mock_bvarnet(family = "gaussian", n_iter = 1000L, n_chains = 1L)

  # Mixed/nodewise fits scale each Gaussian node by its own sd(Y_node).
  eff1 <- set_priors(); eff1$beta$scale <- 2
  eff2 <- set_priors(); eff2$beta$scale <- 20
  obj$priors_effective <- list(eff1, eff2)

  expect_equal(savage_dickey(obj, "beta[2,1]")$prior_density, dnorm(0, 0, 2))
  expect_equal(savage_dickey(obj, "beta[2,2]")$prior_density, dnorm(0, 0, 20))
})

test_that("beta row 1 uses the intercept prior on gaussian/bernoulli nodes", {
  eff <- set_priors(intercept = prior("normal", 0, 5))   # beta stays at 1

  for (fam in c("gaussian", "bernoulli")) {
    obj <- make_mock_bvarnet(family = fam, n_iter = 1000L, n_chains = 1L)
    obj$priors_effective <- rep(list(eff), 2L)

    expect_equal(savage_dickey(obj, "beta[1,1]")$prior_density, dnorm(0, 0, 5))
    expect_equal(savage_dickey(obj, "beta[2,1]")$prior_density, dnorm(0, 0, 1))
  }
})

test_that("beta row 1 uses the beta prior on ordinal nodes", {
  # Ordinal strips the Intercept from X — beta[1,j] is a covariate, and the
  # intercept is absorbed by the cutpoints.
  obj <- make_mock_bvarnet(family = "ordinal", n_iter = 1000L, n_chains = 1L)
  eff <- set_priors(intercept = prior("normal", 0, 5))
  obj$priors_effective <- rep(list(eff), 2L)

  expect_equal(savage_dickey(obj, "beta[1,1]")$prior_density, dnorm(0, 0, 1))
})

test_that("savage_dickey warns and falls back on objects with no effective priors", {
  obj <- make_mock_bvarnet(family = "gaussian", n_iter = 1000L, n_chains = 1L)
  obj$priors_effective <- NULL

  expect_warning(res <- savage_dickey(obj, "phi[1,1]"), "predates")
  expect_equal(res$prior_density, dnorm(0, 0, 0.5))   # phi default, never scaled
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 End-to-end against a real CmdStan fit
# ═══════════════════════════════════════════════════════════════════════════════

test_that("a fitted object's effective priors match the scales sent to Stan", {
  skip_if_no_bvarnet_models()

  sim <- sim_var(N = 4, T_obs = 15, p = 2, K = 1, family = "gaussian", seed = 1)
  sim$data$y_1 <- sim$data$y_1 * 10   # force s_y != 1

  fit <- bvar(id_col = "id", time_col = "t", y_cols = c("y_1", "y_2"),
              K = 1, family = "gaussian", data = sim$data, iter = 200, warmup = 150,
              chains = 1, seed = 1)

  expect_equal(fit$priors_effective[[1]]$beta$scale,      fit$standata$beta_scale)
  expect_equal(fit$priors_effective[[1]]$intercept$scale, fit$standata$intercept_scale)
  expect_gt(fit$standata$intercept_scale, 1)

  # What the user asked for is preserved, unscaled.
  expect_equal(fit$priors$intercept$scale, 1)

  # And the SDDR divides by the scaled prior.
  res <- savage_dickey(fit, "beta[1,1]", method = "logspline")
  expect_equal(res$prior_density, dnorm(0, 0, fit$standata$intercept_scale))
})

test_that("phi priors are never scaled, so phi BFs are unaffected", {
  skip_if_no_bvarnet_models()

  sim <- sim_var(N = 4, T_obs = 15, p = 2, K = 1, family = "gaussian", seed = 2)
  sim$data$y_1 <- sim$data$y_1 * 10

  fit <- bvar(id_col = "id", time_col = "t", y_cols = c("y_1", "y_2"),
              K = 1, family = "gaussian", data = sim$data, iter = 200, warmup = 150,
              chains = 1, seed = 2)

  expect_equal(fit$priors_effective[[1]]$phi$scale, fit$priors$phi$scale)
  expect_equal(savage_dickey(fit, "phi[1,1]")$prior_density, dnorm(0, 0, 0.5))
})
