# Extracted from test-simulate-predict.R:412

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "bvarnet", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
make_predictable_mock <- function(family = "gaussian", n_re = 0L, J = 5L,
                                   T_obs = 20L, p = 2L, q = 1L, K = 1L) {
  set.seed(99L)
  N <- J
  # Make long-format data
  df <- expand.grid(t = seq_len(T_obs), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p)) {
    col <- paste0("y_", j)
    if (family == "bernoulli")     df[[col]] <- rbinom(nrow(df), 1, 0.5)
    else if (family == "gaussian") df[[col]] <- rnorm(nrow(df))
    else                           df[[col]] <- sample(1:3, nrow(df), replace = TRUE)
  }
  for (j in seq_len(q)) df[[paste0("x_", j)]] <- rnorm(nrow(df))

  y_cols <- paste0("y_", seq_len(p))
  x_cols <- if (q > 0) paste0("x_", seq_len(q)) else character(0)

  # Use to_stan_data to get proper matrices
  sd <- to_stan_data(
    data     = df,
    family   = family,
    id_col   = "id",
    time_col = "t",
    y_cols   = y_cols,
    x_cols   = x_cols,
    K        = K,
    re_cols  = if (n_re > 0L) "Intercept" else character(0),
    re_temporal = (n_re > 1L)
  )

  n_fe   <- sd$n_fe
  actual_n_re <- sd$n_re
  PK     <- p * K
  n_iter <- 20L; n_chains <- 2L

  # Build draws array
  beta_nm <- character(0)
  for (node in seq_len(p)) for (fe in seq_len(n_fe))
    beta_nm <- c(beta_nm, sprintf("beta[%d,%d]", fe, node))
  phi_nm <- character(0)
  for (node in seq_len(p)) for (lag_idx in seq_len(PK))
    phi_nm <- c(phi_nm, sprintf("phi[%d,%d]", lag_idx, node))
  par_nms <- c(beta_nm, phi_nm)

  if (family == "gaussian")
    par_nms <- c(par_nms, paste0("sigma[", seq_len(p), "]"))
  if (family == "ordinal") {
    C <- sd$C
    for (node in seq_len(p)) for (k in seq_len(C - 1L))
      par_nms <- c(par_nms, sprintf("kappa[%d,%d]", node, k))
  }
  if (actual_n_re > 0L) {
    for (node in seq_len(p)) for (re in seq_len(actual_n_re))
      par_nms <- c(par_nms, sprintf("sd_u[%d,%d]", node, re))
    for (node in seq_len(p)) for (subj in seq_len(J)) for (re in seq_len(actual_n_re))
      par_nms <- c(par_nms, sprintf("u[%d,%d,%d]", node, subj, re))
  }

  n_par <- length(par_nms)
  draws <- array(rnorm(n_iter * n_chains * n_par, mean = 0.1, sd = 0.3),
                 dim = c(n_iter, n_chains, n_par),
                 dimnames = list(NULL, NULL, par_nms))

  # Make sigma positive, sd_u positive
  sigma_idx <- grep("^sigma\\[", par_nms)
  if (length(sigma_idx) > 0) draws[, , sigma_idx] <- abs(draws[, , sigma_idx]) + 0.1
  sd_u_idx <- grep("^sd_u\\[", par_nms)
  if (length(sd_u_idx) > 0) draws[, , sd_u_idx] <- abs(draws[, , sd_u_idx]) + 0.1

  # Make kappa ordered
  if (family == "ordinal") {
    for (node in seq_len(p)) {
      k_idx <- grep(sprintf("^kappa\\[%d,", node), par_nms)
      for (ch in 1:n_chains) {
        for (it in 1:n_iter) {
          vals <- draws[it, ch, k_idx]
          draws[it, ch, k_idx] <- sort(vals)
        }
      }
    }
  }

  smry <- data.frame(
    variable = par_nms, mean = 0, median = 0, sd = 0.1,
    mad = 0.1, q5 = -0.2, q95 = 0.2, rhat = 1.001,
    ess_bulk = 3000, ess_tail = 2800, stringsAsFactors = FALSE
  )

  structure(
    list(
      draws        = draws,
      summary      = smry,
      diagnostics  = data.frame(num_divergent = integer(n_chains),
                                num_max_treedepth = integer(n_chains),
                                ebfmi = rep(1.0, n_chains)),
      timing       = list(total = 5.0),
      metadata     = list(),
      return_codes = rep(0L, n_chains),
      family       = family,
      standata     = sd,
      priors       = set_priors()
    ),
    class = "bvarnet"
  )
}

# test -------------------------------------------------------------------------
skip_if_not(instantiate::stan_cmdstan_exists())
sim <- sim_var(N = 5, T_obs = 30, p = 2, K = 1,
                 family = "bernoulli", q = 0, seed = 2)
fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "bernoulli",
              iter = 200, warmup = 100, chains = 2, seed = 2)
