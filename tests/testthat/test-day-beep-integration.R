test_that("installed-model sampling supports day/beep data in every family and mixed fits", {
  skip_if_not_installed("cmdstanr")
  skip_if_not(instantiate::stan_cmdstan_exists())
  d <- expand.grid(beep = 1:8, day = 1:2, id = 1:2)
  d$g <- sin(seq_len(nrow(d)))
  d$b <- seq_len(nrow(d)) %% 2
  d$o <- seq_len(nrow(d)) %% 3 + 1
  d$x <- seq_len(nrow(d))/20
  d$g[d$beep == 4] <- NA_real_
  d$b[d$beep == 4] <- NA_real_
  d$o[d$beep == 4] <- NA_real_
  for (family in list("gaussian", "bernoulli", "ordinal", c("gaussian", "bernoulli", "ordinal"))) {
    cols <- if (length(family) == 3) c("g","b","o") else
      switch(family, gaussian = "g", bernoulli = "b", ordinal = "o")
    fit <- suppressWarnings(suppressMessages(bvar(id_col = "id", time_col = "day",
      beep_col = "beep", y_cols = cols, x_cols = "x", data = d, family = family,
      K = 2, center_x = TRUE, re_cols = "Intercept", re_temporal = TRUE,
      fe_interactions = list(c("lag", "x")),
      chains = 2, cores = 2, iter = 100, warmup = 100, adapt_delta = 0.95,
      init = 0, seed = 111, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)))
    expect_equal(fit$standata$n_obs, 12L)
    expect_equal(fit$standata$J, 2L)
    expect_equal(sort(unique(fit$standata$index_obs$beep)), c(3,7,8))
    expect_true(all(fit$return_codes == 0))
    expect_false(any(grepl("beta_centered|kappa_centered", dimnames(fit$draws)[[3]])))
  }
})
