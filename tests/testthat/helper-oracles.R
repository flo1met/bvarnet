# ──────────────────────────────────────────────────────────────────────────────
# helper-oracles.R — independent reference implementations of the three
# likelihoods, plus machinery to compile and optimise the .stan sources.
#
# WHY THIS FILE EXISTS
#
# The rest of the suite tests R-level plumbing against mocks: that the right
# numbers reach Stan and that what comes back is reshaped correctly. None of it
# asks the question these tests ask -- whether the .stan files compute the
# density the package documents.
#
# That question needs an answer from outside the package. The simulator, the
# mock fixtures and the validation checks are all built from the same reading
# of the model as the Stan code, so a sim -> fit -> recover test confirms that
# the pieces agree with each other, not that any of them is right. Wherever
# that shared reading is wrong, every such test still passes.
#
# So the oracles here are a second, independent derivation, and the tests are a
# check for agreement between two things that were written separately. Which
# only works while they stay separate:
#
#   *** The densities below are written from the model definitions in the
#   *** package documentation, NOT transcribed from the .stan files. If you
#   *** ever "fix" an oracle by copying from Stan, you have deleted the test.
#
# See tests/testthat/test-stan-likelihood.R for how they are used.
# ──────────────────────────────────────────────────────────────────────────────


# ── Locating and compiling the Stan sources ──────────────────────────────────

#' Path to a bundled .stan source
#'
#' Source checkouts keep them in `src/stan/`; installed packages ship them
#' alongside the executables in `bin/stan/` (see src/install.libs.R). Returns
#' NA_character_ if neither.
#'
#' The checkout is tried first, and the order matters: an installed copy can be
#' stale relative to the working tree, and silently compiling it would test a
#' different file from the one being edited. `test_path()` is used rather than a
#' path relative to the working directory because the two resolve differently
#' depending on how the suite is invoked -- under `R CMD check` neither reaches
#' a checkout, and the installed copy is then both the only and the right
#' answer, since check installs from the same sources.
#' @noRd
bvarnet_stan_source <- function(model_name) {
  file <- paste0(model_name, ".stan")
  root <- tryCatch(testthat::test_path("..", ".."),
                   error = function(e) NA_character_)
  cand <- c(
    if (!is.na(root))
      normalizePath(file.path(root, "src", "stan", file), mustWork = FALSE),
    system.file("bin", "stan", file, package = "bvarnet")
  )
  cand <- cand[nzchar(cand) & file.exists(cand)]
  if (length(cand) == 0L) NA_character_ else cand[1L]
}

#' Skip unless we can compile a Stan source here
#' @noRd
skip_if_no_stan_compile <- function(model_name) {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("cmdstanr")
  testthat::skip_if_not(instantiate::stan_cmdstan_exists(), "CmdStan not found")
  testthat::skip_if(is.na(bvarnet_stan_source(model_name)),
                    paste0("Stan source for ", model_name, " not available"))
}

#' Compile a bundled .stan source, cached for the session
#'
#' Compiles into a session temp dir so the repo working tree stays clean.
#' Deliberately compiles the SOURCE rather than resolving the shipped
#' executable via `.bvarnet_stan_model()`: this test guards the .stan file in
#' the repo, which in a dev checkout is ahead of the installed binary.
#' @noRd
compile_stan_source <- function(model_name) {
  key <- paste0("stanmod_", model_name)
  if (!exists(key, envir = .test_cache)) {
    dir <- file.path(tempdir(), "bvarnet_oracle_models")
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    assign(key,
           cmdstanr::cmdstan_model(bvarnet_stan_source(model_name), dir = dir),
           envir = .test_cache)
  }
  get(key, envir = .test_cache)
}

#' Penalised-MLE (MAP) fit of a compiled Stan model
#'
#' `jacobian = FALSE` maximises the log density in the CONSTRAINED space, with
#' no change-of-variables adjustment — exactly what a plain `optim()` over the
#' natural parameters computes. That correspondence is what makes the
#' comparison exact rather than approximate.
#' @noRd
stan_map <- function(model_name, standata, seed = 1L) {
  fit <- compile_stan_source(model_name)$optimize(
    data = standata, jacobian = FALSE, seed = seed, refresh = 0
  )
  s <- fit$summary()
  stats::setNames(s$estimate, s$variable)
}


# ── Priors (mirrors set_priors(): 1 = normal, 2 = student-t, 3 = cauchy) ──────

#' @noRd
oracle_lprior <- function(x, fam, loc, scale, df) {
  x <- as.vector(x)
  if (length(x) == 0L) return(0)
  switch(as.character(fam),
    "1" = sum(stats::dnorm(x, loc, scale, log = TRUE)),
    "2" = sum(stats::dt((x - loc) / scale, df, log = TRUE) - log(scale)),
    "3" = sum(stats::dcauchy(x, loc, scale, log = TRUE)),
    stop("Unknown prior family code: ", fam)
  )
}


# ── The three likelihoods ────────────────────────────────────────────────────
#
# Common structure across families. For node j the linear predictor is
#
#     eta[, j] = [X B] %*% c(beta[, j], phi[, j])
#
# i.e. covariates (with the intercept as row 1 of beta for bernoulli/gaussian,
# absent for ordinal) followed by the lagged outcomes on their raw coding.
# Random effects are out of scope here: these tests use n_re = 0.

#' @noRd
oracle_eta <- function(beta, phi, sd_) {
  cbind(sd_$X, sd_$B) %*% rbind(beta, phi)      # n_obs x p
}

#' Adjacent-category log-likelihood for one node
#'
#' log P(Y = c) = (c - 1) * eta - sum_{h < c} kappa_h + const, i.e.
#' log P(Y = c + 1) / P(Y = c) = eta - kappa_c.
#' @noRd
oracle_ll_acat_node <- function(eta, kappa, y, C) {
  ck <- c(0, cumsum(kappa))
  L  <- outer(as.vector(eta), 0:(C - 1)) - rep(ck, each = length(eta))
  mx <- apply(L, 1, max)
  sum(L[cbind(seq_along(y), y)] - (mx + log(rowSums(exp(L - mx)))))
}

#' Unnormalised log posterior, bernoulli nodes
#' @noRd
oracle_lp_bernoulli <- function(pars, sd_) {
  eta <- oracle_eta(pars$beta, pars$phi, sd_)
  ll  <- sum(stats::dbinom(sd_$Y, 1, stats::plogis(eta), log = TRUE))
  ll +
    oracle_lprior(pars$beta[1, ], sd_$prior_intercept_fam, sd_$intercept_loc,
                  sd_$intercept_scale, sd_$intercept_df) +
    oracle_lprior(pars$beta[-1, ], sd_$prior_beta_fam, sd_$beta_loc,
                  sd_$beta_scale, sd_$beta_df) +
    oracle_lprior(pars$phi, sd_$prior_phi_fam, sd_$phi_loc,
                  sd_$phi_scale, sd_$phi_df)
}

#' Unnormalised log posterior, gaussian nodes
#' @noRd
oracle_lp_gaussian <- function(pars, sd_) {
  eta <- oracle_eta(pars$beta, pars$phi, sd_)
  ll  <- sum(stats::dnorm(sd_$Y, eta,
                          rep(pars$sigma, each = nrow(eta)), log = TRUE))
  ll +
    oracle_lprior(pars$beta[1, ], sd_$prior_intercept_fam, sd_$intercept_loc,
                  sd_$intercept_scale, sd_$intercept_df) +
    oracle_lprior(pars$beta[-1, ], sd_$prior_beta_fam, sd_$beta_loc,
                  sd_$beta_scale, sd_$beta_df) +
    oracle_lprior(pars$phi, sd_$prior_phi_fam, sd_$phi_loc,
                  sd_$phi_scale, sd_$phi_df) +
    # half-prior: the truncation constant does not move the mode
    oracle_lprior(pars$sigma, sd_$prior_sigma_fam, sd_$sigma_loc,
                  sd_$sigma_scale, sd_$sigma_df)
}

#' Unnormalised log posterior, ordinal (adjacent-category) nodes
#' @noRd
oracle_lp_ordinal <- function(pars, sd_) {
  eta <- oracle_eta(pars$beta, pars$phi, sd_)
  ll  <- sum(vapply(seq_len(sd_$p), function(j)
    oracle_ll_acat_node(eta[, j], pars$kappa[j, ], sd_$Y[, j], sd_$C),
    numeric(1)))
  ll +
    # no intercept row for ordinal: every beta row is a plain covariate
    oracle_lprior(pars$beta, sd_$prior_beta_fam, sd_$beta_loc,
                  sd_$beta_scale, sd_$beta_df) +
    oracle_lprior(pars$phi, sd_$prior_phi_fam, sd_$phi_loc,
                  sd_$phi_scale, sd_$phi_df) +
    oracle_lprior(pars$kappa, sd_$prior_kappa_fam, sd_$kappa_loc,
                  sd_$kappa_scale, sd_$kappa_df)
}


# ── Packing parameters into a flat vector for optim() ────────────────────────
#
# sigma is optimised on the log scale to keep it positive. Because no Jacobian
# is added, a monotone reparameterisation leaves the argmax unchanged — which
# is the same convention as `jacobian = FALSE` on the Stan side.

#' @noRd
oracle_par_spec <- function(family, sd_) {
  spec <- list(beta = c(sd_$n_fe, sd_$p), phi = c(sd_$p * sd_$K, sd_$p))
  if (family == "gaussian") spec$log_sigma <- c(sd_$p, 1L)
  if (family == "ordinal")  spec$kappa     <- c(sd_$p, sd_$C - 1L)
  spec
}

#' @noRd
oracle_unpack <- function(v, family, sd_) {
  spec <- oracle_par_spec(family, sd_)
  out <- list(); i <- 0L
  for (nm in names(spec)) {
    n <- prod(spec[[nm]])
    out[[nm]] <- matrix(v[i + seq_len(n)], spec[[nm]][1], spec[[nm]][2])
    i <- i + n
  }
  if (family == "gaussian") out$sigma <- as.vector(exp(out$log_sigma))
  out
}

#' Independent MAP fit, by direct optimisation of the oracle log posterior
#' @noRd
oracle_map <- function(family, sd_) {
  lp_fun <- switch(family,
    bernoulli = oracle_lp_bernoulli,
    gaussian  = oracle_lp_gaussian,
    ordinal   = oracle_lp_ordinal,
    stop("Unknown family: ", family)
  )
  spec  <- oracle_par_spec(family, sd_)
  n_par <- sum(vapply(spec, prod, numeric(1)))
  nll   <- function(v) -lp_fun(oracle_unpack(v, family, sd_), sd_)

  op <- stats::optim(rep(0, n_par), nll, method = "BFGS",
                     control = list(maxit = 10000, reltol = 1e-14))
  stopifnot(op$convergence == 0)
  oracle_unpack(op$par, family, sd_)
}


# ── Aligning Stan's flat draw names with the oracle's matrices ───────────────

#' Pull `name[i,j]`-indexed Stan estimates into a matrix of the given dim
#' @noRd
stan_matrix <- function(est, name, dim) {
  m <- matrix(NA_real_, dim[1], dim[2])
  for (i in seq_len(dim[1]))
    for (j in seq_len(dim[2]))
      m[i, j] <- est[[sprintf("%s[%d,%d]", name, i, j)]]
  m
}

#' Everything the comparison checks, as one flat named vector
#' @noRd
oracle_flatten <- function(pars, family, sd_) {
  v <- c(beta = as.vector(pars$beta), phi = as.vector(pars$phi))
  if (family == "gaussian") v <- c(v, sigma = as.vector(pars$sigma))
  if (family == "ordinal")  v <- c(v, kappa = as.vector(pars$kappa))
  v
}

#' @noRd
stan_flatten <- function(est, family, sd_) {
  beta <- stan_matrix(est, "beta", c(sd_$n_fe, sd_$p))
  phi  <- stan_matrix(est, "phi",  c(sd_$p * sd_$K, sd_$p))
  v <- c(beta = as.vector(beta), phi = as.vector(phi))
  if (family == "gaussian")
    v <- c(v, sigma = unname(est[paste0("sigma[", seq_len(sd_$p), "]")]))
  if (family == "ordinal")
    v <- c(v, kappa = as.vector(stan_matrix(est, "kappa",
                                            c(sd_$p, sd_$C - 1L))))
  v
}


# ── Test data ────────────────────────────────────────────────────────────────

#' Stan data for a likelihood-oracle test: p = 2 nodes, K = 1, no random effects
#'
#' p = 2 rather than 1 so the per-node loop and the `phi[, node]` column
#' indexing are exercised — a transposed phi would pass with a single node.
#'
#' For ordinal, `kappa` is deliberately NON-MONOTONE. That is the configuration
#' the `ordered` declaration could not represent, and it is the reason this
#' fixture is not just `get_sim_ordinal()`.
#' @noRd
make_oracle_standata <- function(family, C = 4L, N = 8L, T_obs = 40L,
                                 seed = 99L) {
  kappa <- if (family == "ordinal")
    list(c(1.6, 0.2, 1.1), c(0.9, -0.4, 1.2)) else NULL

  sim <- sim_var(N = N, T_obs = T_obs, p = 2L, K = 1L, q = 2L,
                 family = family, C = C, kappa = kappa, seed = seed)

  to_stan_data(
    data = sim$data, family = family,
    id_col = "id", time_col = "t",
    y_cols = paste0("y_", 1:2), x_cols = paste0("x_", 1:2),
    K = 1L, skip_lag = FALSE
  )
}
