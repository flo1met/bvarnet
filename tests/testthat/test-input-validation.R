# ──────────────────────────────────────────────────────────────────────────────
# test-input-validation.R — error condition and boundary case tests
#
# Follows the bgms pattern: every user-facing function should give
# informative errors on invalid input, tested with expect_error(regexp = ...).
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 sim_var() validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("sim_var rejects invalid family", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 2, family = "poisson"),
    "arg"
  )
})


test_that("sim_var rejects N < 1", {
  expect_error(
    sim_var(N = 0, T_obs = 20, p = 2, family = "bernoulli", seed = 1),
    "N >= 1"
  )
})


test_that("sim_var rejects T_obs < K + 1", {
  expect_error(
    sim_var(N = 3, T_obs = 1, p = 2, K = 1, family = "bernoulli", seed = 1),
    "T_obs >= K"
  )
})


test_that("sim_var rejects p < 1", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 0, family = "bernoulli", seed = 1),
    "p >= 1"
  )
})


test_that("sim_var rejects wrong alpha length", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 3, family = "gaussian",
            alpha = c(1, 2), seed = 1, burnin = 0),
    "alpha"
  )
})


test_that("sim_var rejects wrong Phi dimensions", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 3, K = 1,
            family = "gaussian", Phi = matrix(0, 2, 2), seed = 1, burnin = 0),
    "Phi"
  )
})


test_that("sim_var rejects negative sigma", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 2, family = "gaussian",
            sigma = c(-1, 1), seed = 1, burnin = 0),
    "sigma > 0"
  )
})


test_that("sim_var rejects C < 2 for ordinal", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 2, family = "ordinal",
            C = 1, seed = 1, burnin = 0),
    "C >= 2"
  )
})


test_that("sim_var rejects negative burnin", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 2, family = "bernoulli",
            burnin = -1, seed = 1),
    "burnin >= 0"
  )
})


test_that("sim_var rejects unordered kappa", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 2, family = "ordinal", C = 3,
            kappa = list(c(1, 0), c(0, 1)),  # first one is unordered
            seed = 1, burnin = 0),
    "unsorted|kappa"
  )
})


test_that("sim_var rejects wrong kappa length", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 2, family = "ordinal", C = 4,
            kappa = list(c(-1, 0), c(-1, 0, 1)),  # first has 2 instead of 3
            seed = 1, burnin = 0),
    "kappa"
  )
})


test_that("sim_var rejects wrong gamma dimensions", {
  expect_error(
    sim_var(N = 3, T_obs = 20, p = 3, q = 2,
            family = "gaussian",
            gamma = matrix(0, 1, 3),  # should be 2x3
            seed = 1, burnin = 0),
    "gamma"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 to_stan_data() validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("to_stan_data rejects invalid family", {
  df <- make_test_df(N = 3, T_obs = 10, p = 2)
  expect_error(
    to_stan_data(df, "poisson", "id", "t", paste0("y_", 1:2), character(0), K = 1),
    "arg"
  )
})


test_that("to_stan_data where ordinal Y has values < 1 gives error", {
  df <- make_test_df(N = 3, T_obs = 10, p = 2, family = "ordinal")
  df$y_1 <- df$y_1 - 1  # now contains 0s

  expect_error(
    to_stan_data(df, "ordinal", "id", "t", paste0("y_", 1:2), character(0), K = 1),
    "< 1"
  )
})


test_that("to_stan_data stops when all observations are removed", {
  # A very small dataset where NA + K causes all data to be removed
  df <- data.frame(id = c(1, 1), t = c(1, 2), y_1 = c(NA, 0L), y_2 = c(1L, NA))

  expect_error(
    to_stan_data(df, "bernoulli", "id", "t", paste0("y_", 1:2), character(0), K = 1),
    "All observations removed"
  )
})


test_that("to_stan_data rejects re_cols not found in X", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 1, family = "bernoulli")

  expect_error(
    to_stan_data(df, "bernoulli", "id", "t",
                 paste0("y_", 1:2), "x_1",
                 re_cols = "nonexistent_col", K = 1),
    "re_cols not found"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 normalize_terms() validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("normalize_terms rejects non-list input", {
  expect_error(
    bvarnet:::normalize_terms("x_1"),
    "list"
  )
})


test_that("normalize_terms rejects single-factor terms", {
  expect_error(
    bvarnet:::normalize_terms(list(c("x_1"))),
    "at least 2"
  )
})


test_that("normalize_terms rejects duplicate factors", {
  expect_error(
    bvarnet:::normalize_terms(list(c("x_1", "x_1"))),
    "Duplicate"
  )
})


test_that("normalize_terms rejects empty factor names", {
  expect_error(
    bvarnet:::normalize_terms(list(c("x_1", ""))),
    "Empty"
  )
})


test_that("normalize_terms accepts valid terms", {
  result <- bvarnet:::normalize_terms(list(c("x_1", "x_2"), c("lag", "x_1")))
  expect_length(result, 2)
  expect_equal(result[[1]], c("x_1", "x_2"))
})


test_that("normalize_terms returns empty list for NULL", {
  expect_equal(bvarnet:::normalize_terms(NULL), list())
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 Interaction term validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("add_terms_to_X rejects unknown factor names", {
  X <- matrix(1, 10, 2, dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0, 10, 2, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))

  expect_error(
    bvarnet:::add_terms_to_X(X, B, list(c("x_1", "nonexistent"))),
    "Unknown factor"
  )
})


test_that("add_terms_to_X rejects duplicate columns", {
  X <- matrix(1, 10, 3, dimnames = list(NULL, c("Intercept", "x_1", "x_1:x_2")))
  B <- matrix(0, 10, 2, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))

  # x_1:x_2 already exists in X, so adding the interaction should fail
  # First we need x_2 in X for the term to be valid
  X2 <- cbind(X, x_2 = 1)
  expect_error(
    bvarnet:::add_terms_to_X(X2, B, list(c("x_1", "x_2"))),
    "already exist"
  )
})


test_that("add_re_interactions_from_X rejects missing FE interactions", {
  Z <- matrix(0, 10, 0)
  X <- matrix(1, 10, 2, dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0, 10, 2, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))

  # Ask for RE interaction on c("x_1","x_2") which has not been added to FE
  expect_error(
    bvarnet:::add_re_interactions_from_X(Z, X, B, list(c("x_1", "x_2"))),
    "not present in fixed effects"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §5 Edge cases
# ═══════════════════════════════════════════════════════════════════════════════

test_that("single subject (N=1) works", {
  result <- sim_var(N = 1, T_obs = 30, p = 2, family = "gaussian",
                    seed = 400, burnin = 10)
  expect_equal(nrow(result$data), 30)
  expect_equal(result$truth$N, 1L)
  expect_equal(nrow(result$truth$alpha_i), 1)
})


test_that("single node (p=1) works", {
  result <- sim_var(N = 3, T_obs = 20, p = 1, family = "bernoulli",
                    seed = 401, burnin = 10)
  expect_equal(ncol(result$truth$Phi), 1)
  expect_true("y_1" %in% names(result$data))
})


test_that("minimum ordinal C = 2 works", {
  result <- sim_var(N = 3, T_obs = 20, p = 2, family = "ordinal",
                    C = 2, seed = 402, burnin = 10)
  expect_true(all(result$data$y_1 %in% c(1L, 2L)))
})


test_that("minimum T_obs = K + 1 works", {
  result <- sim_var(N = 3, T_obs = 2, p = 2, K = 1,
                    family = "gaussian", seed = 403, burnin = 10)
  expect_equal(nrow(result$data), 3 * 2)
})


test_that("to_stan_data works with minimal data (1 subject, K+1 rows)", {
  df <- data.frame(id = 1L, t = 1:3, y_1 = rnorm(3), y_2 = rnorm(3))

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     c("y_1", "y_2"), character(0), K = 1)

  expect_equal(sd$n_obs, 2)  # T=3, K=1 → 2 modeled rows
  expect_equal(sd$J, 1)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §6 bvar() input validation (no Stan required)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("bvar rejects invalid family string", {
  df <- make_test_df(N = 3, T_obs = 10, p = 2)

  expect_error(
    bvar(id_col = "id", time_col = "t",
         y_cols = paste0("y_", 1:2), x_cols = character(0),
         data = df, family = "poisson", K = 1),
    "arg"
  )
})
