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
#' p=2 outcomes, K=1 lag, n_fe=2 (Intercept + x_1).
#' For ordinal: C=3 (2 cutpoints).  For gaussian: sigma included.
#' When \code{n_re > 0}, \code{sd_u} and \code{u} draw columns are added
#' with \code{J} subjects and \code{n_re} random-effect columns.
#'
#' \code{family} can be a scalar (applied to all p nodes) or a character
#' vector of length p for mixed-family mocks. When mixed, sigma columns are
#' only generated for gaussian nodes and kappa columns only for ordinal nodes.
#' Ordinal nodes get \code{beta[1,j] = NA_real_} sentinel (D4).
make_mock_bvarnet <- function(family   = "bernoulli",
                               n_iter   = 20L,
                               n_chains = 2L,
                               n_re     = 0L,
                               J        = 5L) {
  p <- 2L; K <- 1L; n_fe <- 2L

  # Normalise family to a named vector of length p
  y_cols <- c("y_1", "y_2")
  if (length(family) == 1L) {
    family_vec <- setNames(rep(family, p), y_cols)
  } else {
    stopifnot(length(family) == p)
    family_vec <- setNames(family, y_cols)
  }

  beta_nm <- c("beta[1,1]", "beta[2,1]", "beta[1,2]", "beta[2,2]")
  phi_nm  <- c("phi[1,1]",  "phi[2,1]",  "phi[1,2]",  "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm)

  # sigma only for gaussian nodes
  gauss_idx <- which(family_vec == "gaussian")
  if (length(gauss_idx) > 0L) {
    sigma_nm <- paste0("sigma[", gauss_idx, "]")
    par_nms <- c(par_nms, sigma_nm)
  }

  # kappa only for ordinal nodes (C=3 → 2 cutpoints each)
  ord_idx <- which(family_vec == "ordinal")
  if (length(ord_idx) > 0L) {
    kappa_nm <- character(0)
    for (node in ord_idx)
      for (k in 1:2)
        kappa_nm <- c(kappa_nm, sprintf("kappa[%d,%d]", node, k))
    par_nms <- c(par_nms, kappa_nm)
  }

  # sd_u parameters: sd_u[node,re]
  if (n_re > 0L) {
    sd_u_nm <- character(0)
    for (node in seq_len(p))
      for (re in seq_len(n_re))
        sd_u_nm <- c(sd_u_nm, sprintf("sd_u[%d,%d]", node, re))
    par_nms <- c(par_nms, sd_u_nm)
  }

  # u parameters: array[p] matrix[J, n_re] → u[node, subject, re]
  u_nm <- character(0)
  if (n_re > 0L) {
    for (node in seq_len(p))
      for (subj in seq_len(J))
        for (re in seq_len(n_re))
          u_nm <- c(u_nm, sprintf("u[%d,%d,%d]", node, subj, re))
    par_nms <- c(par_nms, u_nm)
  }

  n_par <- length(par_nms)
  set.seed(42L)
  draws <- array(
    rnorm(n_iter * n_chains * n_par),
    dim      = c(n_iter, n_chains, n_par),
    dimnames = list(NULL, NULL, par_nms)
  )

  # Make sigma draws positive
  if (length(gauss_idx) > 0L) {
    sigma_idx_arr <- grep("^sigma\\[", par_nms)
    draws[, , sigma_idx_arr] <- abs(draws[, , sigma_idx_arr]) + 0.1
  }

  # Make kappa draws ordered (ascending per node)
  if (length(ord_idx) > 0L) {
    for (node in ord_idx) {
      k1 <- paste0("kappa[", node, ",1]")
      k2 <- paste0("kappa[", node, ",2]")
      draws[, , k1] <- -1 + runif(n_iter * n_chains, -0.2, 0.2)
      draws[, , k2] <-  1 + runif(n_iter * n_chains, -0.2, 0.2)
    }
  }

  # Set ordinal beta[1,j] to NA sentinel (D4)
  for (j in ord_idx) {
    nm <- paste0("beta[1,", j, "]")
    draws[, , nm] <- NA_real_
  }

  # Make sd_u draws positive (half-prior)
  if (n_re > 0L) {
    sd_u_idx <- grep("^sd_u\\[", par_nms)
    draws[, , sd_u_idx] <- abs(draws[, , sd_u_idx])
  }

  smry <- data.frame(
    variable  = par_nms,
    rhat      = 1.001,
    ess_bulk  = 3000,
    ess_tail  = 2800,
    stringsAsFactors = FALSE
  )

  n_obs <- 10L
  Y <- matrix(0L, n_obs, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0,  n_obs, n_fe,
              dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0,  n_obs, p * K,
              dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))

  if (n_re > 0L) {
    re_colnames <- paste0("z_", seq_len(n_re))
    Z <- matrix(rnorm(n_obs * n_re), n_obs, n_re,
                dimnames = list(NULL, re_colnames))
  } else {
    Z <- matrix(0, n_obs, 0L)
  }

  sd_list <- list(
    p = p, K = K, n_fe = n_fe, n_re = n_re, n = n_obs, n_obs = n_obs,
    J = J, Y = Y, X = X, B = B, Z = Z,
    id = rep(seq_len(J), length.out = n_obs),
    id_levels = as.character(seq_len(J)),
    x_center_means = NULL,
    fe_interaction_terms = list(),
    design_spec = list(
      id_col   = "id",
      time_col = "t",
      y_cols   = c("y_1", "y_2"),
      x_cols   = "x_1",
      center_x = FALSE,
      fe_interactions = NULL,
      re_interactions = NULL,
      re_cols    = character(0),
      re_temporal = FALSE,
      K          = K,
      skip_lag   = TRUE,
      na_action  = "listwise"
    )
  )
  if (length(ord_idx) > 0L) sd_list$C <- 3L
  if (length(ord_idx) > 0L) {
    C_per_node <- setNames(rep(NA_integer_, p), y_cols)
    C_per_node[ord_idx] <- 3L
    sd_list$C_per_node <- C_per_node
  }

  structure(
    list(
      draws        = draws,
      convergence  = smry,
      diagnostics  = data.frame(
        num_divergent     = integer(n_chains),
        num_max_treedepth = integer(n_chains),
        ebfmi             = rep(1.0, n_chains)
      ),
      timing       = list(total = 5.0),
      metadata     = list(),
      return_codes = rep(0L, n_chains),
      family       = family_vec,
      standata     = sd_list,
      priors       = set_priors()
    ),
    class = "bvarnet"
  )
}
