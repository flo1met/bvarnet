# ──────────────────────────────────────────────────────────────────────────────
# helper-fixtures.R — session-cached fixtures & reusable test utilities
#
# This file is automatically sourced by testthat before any test file.
# Follows the bgms pattern of lazy-cached fixtures and shared helpers.
# ──────────────────────────────────────────────────────────────────────────────

# ── Load package functions ────────────────────────────────────────────────────
# When running via devtools::test() or test_check(), the package is loaded
# automatically. For test_dir() or interactive use, we load all source files.
if (!isNamespace(tryCatch(asNamespace("bvarnet"), error = function(e) FALSE))) {
  pkg_root <- normalizePath(file.path(dirname(sys.frame(1)$ofile %||% "."), "..", ".."))
  r_files <- list.files(file.path(pkg_root, "R"), pattern = "\\.R$", full.names = TRUE)
  for (f in r_files) source(f)
}

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
