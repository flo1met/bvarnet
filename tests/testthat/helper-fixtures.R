# ──────────────────────────────────────────────────────────────────────────────
# helper-fixtures.R — session-cached fixtures & reusable test utilities
#
# This file is automatically sourced by testthat before any test file.
# Follows the bgms pattern of lazy-cached fixtures and shared helpers.
# ──────────────────────────────────────────────────────────────────────────────

# ── Load package functions ────────────────────────────────────────────────────
# When running via devtools::test() or test_check() the package is loaded
# automatically. When running via test_dir() or Rscript we source all R files
# so that the in-development versions are used instead of any installed copy.
# We always source to ensure local edits take precedence over the installed pkg.
local({
  here <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
  pkg_root <- normalizePath(file.path(here, "..", ".."))
  r_dir    <- file.path(pkg_root, "R")
  if (dir.exists(r_dir)) {
    for (f in list.files(r_dir, pattern = "\\.R$", full.names = TRUE)) source(f)
  }
})

# ── Session cache environment ────────────────────────────────────────────────
.test_cache <- new.env(parent = emptyenv())

# ── Fixture getters (each cached once per test session) ──────────────────────

#' Get a simulated bernoulli dataset (cached)
get_sim_bernoulli <- function() {

  if (!exists("sim_bernoulli", envir = .test_cache)) {
    .test_cache$sim_bernoulli <- sim_var(
      N = 10, T_obs = 50, p = 3, K = 1,
      family = "bernoulli", q = 2, seed = 42
    )
  }
  .test_cache$sim_bernoulli
}

#' Get a simulated gaussian dataset (cached)
get_sim_gaussian <- function() {
  if (!exists("sim_gaussian", envir = .test_cache)) {
    .test_cache$sim_gaussian <- sim_var(
      N = 10, T_obs = 50, p = 3, K = 1,
      family = "gaussian", q = 2, seed = 43
    )
  }
  .test_cache$sim_gaussian
}

#' Get a simulated ordinal dataset (cached)
get_sim_ordinal <- function() {
  if (!exists("sim_ordinal", envir = .test_cache)) {
    .test_cache$sim_ordinal <- sim_var(
      N = 10, T_obs = 50, p = 3, K = 1,
      family = "ordinal", q = 2, C = 5, seed = 44
    )
  }
  .test_cache$sim_ordinal
}

#' Get stan data for a bernoulli simulation (cached)
get_standata_bernoulli <- function() {
  if (!exists("standata_bernoulli", envir = .test_cache)) {
    sim <- get_sim_bernoulli()
    .test_cache$standata_bernoulli <- to_stan_data(
      data     = sim$data,
      family   = "bernoulli",
      id_col   = "id",
      time_col = "t",
      y_cols   = paste0("y_", 1:3),
      x_cols   = paste0("x_", 1:2),
      K        = 1
    )
  }
  .test_cache$standata_bernoulli
}

#' Get stan data for a gaussian simulation (cached)
get_standata_gaussian <- function() {
  if (!exists("standata_gaussian", envir = .test_cache)) {
    sim <- get_sim_gaussian()
    .test_cache$standata_gaussian <- to_stan_data(
      data     = sim$data,
      family   = "gaussian",
      id_col   = "id",
      time_col = "t",
      y_cols   = paste0("y_", 1:3),
      x_cols   = paste0("x_", 1:2),
      K        = 1
    )
  }
  .test_cache$standata_gaussian
}

#' Get stan data for an ordinal simulation (cached)
get_standata_ordinal <- function() {
  if (!exists("standata_ordinal", envir = .test_cache)) {
    sim <- get_sim_ordinal()
    .test_cache$standata_ordinal <- to_stan_data(
      data     = sim$data,
      family   = "ordinal",
      id_col   = "id",
      time_col = "t",
      y_cols   = paste0("y_", 1:3),
      x_cols   = paste0("x_", 1:2),
      K        = 1
    )
  }
  .test_cache$standata_ordinal
}


# ── Reusable test helpers ────────────────────────────────────────────────────

#' Check if all values lie within [lo, hi]
values_in_range <- function(x, lo, hi) {
  all(x >= lo & x <= hi, na.rm = TRUE)
}

#' Check if a matrix is symmetric
is_symmetric <- function(m, tol = .Machine$double.eps^0.5) {
  is.matrix(m) && nrow(m) == ncol(m) && all(abs(m - t(m)) < tol)
}

#' Build a minimal long-format data frame for testing to_stan_data
make_test_df <- function(N = 5, T_obs = 20, p = 2, q = 0,
                         family = "bernoulli", seed = 123) {
  set.seed(seed)
  df <- expand.grid(t = 1:T_obs, id = 1:N)
  df <- df[order(df$id, df$t), ]

  for (j in seq_len(p)) {
    col <- paste0("y_", j)
    if (family == "bernoulli") {
      df[[col]] <- rbinom(nrow(df), 1, 0.5)
    } else if (family == "gaussian") {
      df[[col]] <- rnorm(nrow(df))
    } else if (family == "ordinal") {
      df[[col]] <- sample(1:5, nrow(df), replace = TRUE)
    }
  }
  for (j in seq_len(q)) {
    df[[paste0("x_", j)]] <- rnorm(nrow(df))
  }
  df
}

#' Build a mock bvarnet object without running Stan
#'
#' Returns a minimal bvarnet list with a synthetic 3D draws array, a matching
#' summary data.frame, and stubbed diagnostics/timing/metadata fields.
#' p=2 outcomes, K=1 lag, n_fe=2 (Intercept + x_1), n_re=0.
#' For ordinal: C=3 (2 cutpoints).  For gaussian: sigma included.
make_mock_bvarnet <- function(family = "bernoulli",
                               n_iter  = 20L,
                               n_chains = 2L) {
  p <- 2L; K <- 1L; n_fe <- 2L; n_re <- 0L

  beta_nm <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]")
  phi_nm  <- c("phi[1,1]",  "phi[2,1]",  "phi[1,2]",  "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm)

  if (family == "gaussian")
    par_nms <- c(par_nms, "sigma[1]", "sigma[2]")
  if (family == "ordinal")
    par_nms <- c(par_nms, "kappa[1,1]", "kappa[2,1]", "kappa[1,2]", "kappa[2,2]")

  n_par <- length(par_nms)
  set.seed(42L)
  draws <- array(
    rnorm(n_iter * n_chains * n_par),
    dim      = c(n_iter, n_chains, n_par),
    dimnames = list(NULL, NULL, par_nms)
  )

  smry <- data.frame(
    variable  = par_nms,
    mean      = 0,
    median    = 0,
    sd        = 0.1,
    mad       = 0.1,
    q5        = -0.2,
    q95       = 0.2,
    rhat      = 1.001,
    ess_bulk  = 3000,
    ess_tail  = 2800,
    stringsAsFactors = FALSE
  )

  Y <- matrix(0L, 10L, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0,  10L, n_fe,
              dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0,  10L, p * K,
              dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))
  Z <- matrix(0,  10L, n_re)

  sd_list <- list(
    p = p, K = K, n_fe = n_fe, n_re = n_re, n = 10L,
    Y = Y, X = X, B = B, Z = Z
  )
  if (family == "ordinal") sd_list$C <- 3L

  structure(
    list(
      draws        = draws,
      summary      = smry,
      diagnostics  = data.frame(
        num_divergent     = integer(n_chains),
        num_max_treedepth = integer(n_chains),
        ebfmi             = rep(1.0, n_chains)
      ),
      timing       = list(total = 5.0),
      metadata     = list(),
      return_codes = rep(0L, n_chains),
      family       = family,
      standata     = sd_list
    ),
    class = "bvarnet"
  )
}
