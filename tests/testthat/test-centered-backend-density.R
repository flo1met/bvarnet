test_that("compiled targets retain the full reference density for every prior family", {
  for (family in c("gaussian", "bernoulli", "ordinal")) {
    skip_if_no_stan_compile(.bvarnet_model_for_family(family))
    skip_if_not_installed("Rcpp")
    for (n_re in c(0L, 2L)) for (prior_family in 1:3) {
      fixture <- centered_fixture(family, n_re, prior_family)
      sd <- fixture$data
      fit <- centered_fixed_fit(sd, fixture$pars, family)
      # Different points exercise a draw-dependent prior shift, not a constant
      # normalization correction that happens to match a single evaluation.
      for (offset in c(0, .13)) {
        pars <- fixture$pars
        pars$phi <- pars$phi + offset
        pars$beta <- pars$beta - offset / 2
        init <- .bvarnet_centered_init(pars, sd, family, 1)[[1]]
        uv <- fit$unconstrain_variables(init)
        expected <- reference_density(pars, sd, family)
        actual <- fit$log_prob(uv, jacobian = FALSE)
        info <- paste(family, n_re, prior_family, offset)
        expect_equal(actual, unname(expected["total"]), tolerance = 1e-9, info = info)
        # With the full-target identity established, isolate the prior too.
        expect_equal(actual - expected[["likelihood"]], expected[["prior"]], tolerance = 1e-9, info = info)
        jacobian <- sum(log(pars$sd_u)) + if (family == "gaussian") sum(log(pars$sigma)) else 0
        expect_equal(fit$log_prob(uv, jacobian = TRUE) - actual, jacobian, tolerance = 1e-9, info = info)
        public <- natural_stan_parameters(fit$constrain_variables(uv), sd)
        expect_equal(public$beta, pars$beta, tolerance = 1e-12, ignore_attr = TRUE)
        expect_equal(public$phi, pars$phi, tolerance = 1e-12, ignore_attr = TRUE)
        if (family == "ordinal") expect_equal(public$kappa, pars$kappa, tolerance = 1e-12, ignore_attr = TRUE)
        if (n_re > 0) for (j in seq_len(sd$p))
          expect_equal(matrix(public$u[j, , ], sd$J, sd$n_re),
                       sweep(matrix(pars$z_u[j, , ], sd$J, sd$n_re), 2, pars$sd_u[j, ], "*"),
                       tolerance = 1e-12, ignore_attr = TRUE)
      }
    }
  }
})

test_that("compiled gradients obey the full coordinate mapping and positivity Jacobians", {
  for (family in c("gaussian", "bernoulli", "ordinal")) {
    skip_if_no_stan_compile(.bvarnet_model_for_family(family))
    skip_if_not_installed("Rcpp")
    fixture <- centered_fixture(family)
    sd <- fixture$data
    fit <- centered_fixed_fit(sd, fixture$pars, family)
    uv <- fit$unconstrain_variables(.bvarnet_centered_init(fixture$pars, sd, family, 1)[[1]])
    for (jacobian in c(FALSE, TRUE)) {
      target <- function(v) {
        sampled <- fit$constrain_variables(v, transformed_parameters = FALSE, generated_quantities = FALSE)
        pars <- reference_from_sampled(sampled, sd, family)
        value <- reference_density(pars, sd, family)[["total"]]
        if (jacobian) value <- value + sum(log(pars$sd_u)) +
          if (family == "gaussian") sum(log(pars$sigma)) else 0
        value
      }
      gradient <- vapply(seq_along(uv), function(i) {
        step <- 1e-5 * max(1, abs(uv[i]))
        plus <- minus <- uv; plus[i] <- plus[i] + step; minus[i] <- minus[i] - step
        (target(plus) - target(minus)) / (2 * step)
      }, numeric(1))
      actual <- fit$grad_log_prob(uv, jacobian = jacobian)
      expect_true(all(abs(actual - gradient) <= 2e-5 * (1 + abs(gradient))), info = family)
    }
  }
})

test_that("location reconstruction has the fixed signs, including ordinal n_fe=0", {
  for (family in c("gaussian", "bernoulli", "ordinal")) {
    skip_if_no_stan_compile(.bvarnet_model_for_family(family))
    skip_if_not_installed("Rcpp")
    fixture <- centered_fixture(family, 0L)
    sd <- fixture$data
    sd$p <- sd$K <- sd$J <- 1L; sd$n_obs <- 2L; sd$id <- c(1L, 1L)
    sd$X <- if (family == "ordinal") matrix(c(8, 12)) else cbind(1, c(8, 12))
    sd$n_fe <- ncol(sd$X)
    sd$B <- matrix(if (family == "bernoulli") c(0, 1) else c(1, 3))
    sd$Z <- matrix(numeric(), 2, 0); sd$Y <- matrix(if (family == "ordinal") c(1L, 4L) else c(0, 1))
    pars <- list(beta = matrix(if (family == "ordinal") .3 else c(if (family == "bernoulli") 1.8 else 1.2, .3), sd$n_fe),
                 phi = matrix(.4), sd_u = matrix(numeric(), 1, 0), z_u = array(numeric(), c(1, 1, 0)))
    if (family == "gaussian") pars$sigma <- 1
    if (family == "ordinal") pars$kappa <- matrix(c(5.4, 4, 4.9), 1)
    fit <- centered_fixed_fit(sd, pars, family)
    values <- natural_stan_parameters(fit$constrain_variables(fit$unconstrain_variables(.bvarnet_centered_init(pars, sd, family, 1)[[1]])), sd)
    if (family == "ordinal") {
      expect_equal(as.vector(values$kappa_centered), c(1.6, .2, 1.1), tolerance = 1e-12)
      expect_equal(as.vector(values$kappa), c(5.4, 4, 4.9), tolerance = 1e-12)
      sd$X <- matrix(numeric(), 2, 0); sd$n_fe <- 0L; pars$beta <- matrix(numeric(), 0, 1)
      empty_fit <- centered_fixed_fit(sd, pars, family)
      uv <- empty_fit$unconstrain_variables(.bvarnet_centered_init(pars, sd, family, 1)[[1]])
      expect_equal(empty_fit$log_prob(uv, jacobian = FALSE), reference_density(pars, sd, family)[["total"]], tolerance = 1e-9)
    } else {
      expect_equal(values$beta_centered[1, 1], 5, tolerance = 1e-12)
      expect_equal(values$beta[1, 1], if (family == "bernoulli") 1.8 else 1.2, tolerance = 1e-12)
    }
  }
})

test_that("modeled-row means align for user centering, interactions, missingness and gaps", {
  skip_if_no_stan_compile("model_gaussian")
  skip_if_not_installed("Rcpp")
  df <- data.frame(id = rep(1:2, each = 8), t = rep(c(1:4, 6:9), 2),
                   y1 = seq_len(16) / 10, y2 = rev(seq_len(16)) / 7,
                   x1 = rep(c(-3, -2, -1, 0, 0, 1, 2, 3), 2), x2 = rep(c(1, 3), 8))
  df$x2[4] <- NA_real_
  for (center_x in c(FALSE, TRUE)) for (skip_lag in c(FALSE, TRUE)) {
    sd <- to_stan_data(df, "gaussian", "id", "t", c("y1", "y2"), c("x1", "x2"),
                       K = 2, center_x = center_x, skip_lag = skip_lag,
                       fe_interactions = list(c("x1", "x2"), c("lag", "x1")),
                       re_cols = "x1", re_temporal = TRUE,
                       re_interactions = list(c("x1", "x2")))
    Z <- sd$Z; X <- sd$X; B <- sd$B; priors <- attr(sd, "priors_effective")
    sd <- sd[intersect(names(sd), .stan_data_names("gaussian"))]
    pars <- list(beta = matrix(seq_len(sd$n_fe * sd$p) / 50, sd$n_fe, sd$p),
                 phi = matrix(seq_len(sd$p * sd$K * sd$p) / 100, sd$p * sd$K, sd$p),
                 sd_u = matrix(.5, sd$p, sd$n_re), z_u = array(.1, c(sd$p, sd$J, sd$n_re)), sigma = c(1, 1.2))
    fit <- centered_fixed_fit(sd, pars, "gaussian")
    mapped <- .bvarnet_centered_init(pars, sd, "gaussian", 1)[[1]]
    uv <- fit$unconstrain_variables(mapped)
    expect_equal(fit$constrain_variables(uv)$beta, pars$beta, tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(fit$log_prob(uv, jacobian = FALSE), reference_density(pars, sd, "gaussian")[["total"]], tolerance = 1e-9)
    expect_identical(sd$Z, Z); expect_identical(sd$X, X); expect_identical(sd$B, B)
    expect_equal(priors$intercept$scale, sd$intercept_scale)
    expect_true(any(colMeans(sd$B) != 0))
  }
})
