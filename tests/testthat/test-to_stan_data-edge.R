# ──────────────────────────────────────────────────────────────────────────────
# test-to_stan_data-edge.R — edge-case coverage for to_stan_data()
#
# Tests boundary dimensions, family-specific casting, design matrix mechanics,
# and the nodewise path (.to_stan_data_node).
# See dev/to_stan_data_coverage_plan.md for rationale.
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 Single subject (N=1, J=1)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("to_stan_data works with a single subject (J=1)", {
  df <- make_test_df(N = 1, T_obs = 20, p = 2, family = "gaussian")
  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_equal(sd$J, 1L)
  expect_equal(sd$n_obs, 19L)
  expect_true(all(sd$id == 1L))
  expect_equal(ncol(sd$Z), 0L)
})

# ═══════════════════════════════════════════════════════════════════════════════
# §2 Single node (p=1)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("to_stan_data works with a single node (p=1)", {
  df <- make_test_df(N = 3, T_obs = 15, p = 1, family = "bernoulli")
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     "y_1", character(0), K = 1)

  expect_equal(sd$p, 1L)
  expect_equal(ncol(sd$B), 1L)
  expect_equal(ncol(sd$Y), 1L)
  expect_equal(colnames(sd$B), "lag1_y_1")
})

# ═══════════════════════════════════════════════════════════════════════════════
# §3 Large K with short T (subject dropping)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("large K drops subjects with T_i <= K", {
  # 2 subjects: one with 10 rows (survives K=5), one with 4 rows (dropped)
  df_long  <- make_test_df(N = 1, T_obs = 10, p = 2, family = "gaussian", seed = 1)
  df_short <- make_test_df(N = 1, T_obs = 4,  p = 2, family = "gaussian", seed = 2)
  df_short$id <- 2L
  df <- rbind(df_long, df_short)

  sd <- suppressMessages(
    to_stan_data(df, "gaussian", "id", "t",
                 paste0("y_", 1:2), character(0), K = 5)
  )

  expect_equal(sd$J, 2L)
  # only subject 1 contributes rows: 10 - 5 = 5
  expect_equal(sd$n_obs, 5L)
  expect_true(all(sd$id == 1L))
  expect_equal(ncol(sd$B), 2 * 5)  # p * K = 10
})

# ═══════════════════════════════════════════════════════════════════════════════
# §4 K > 1 lag column ordering
# ═══════════════════════════════════════════════════════════════════════════════

test_that("K=2 lag columns are correctly ordered and filled", {
  df <- data.frame(
    id = rep(1L, 6),
    t  = 1:6,
    y_1 = c(10, 20, 30, 40, 50, 60),
    y_2 = c(11, 21, 31, 41, 51, 61)
  )
  sd <- to_stan_data(df, "gaussian", "id", "t",
                     c("y_1", "y_2"), character(0), K = 2)

  expect_equal(colnames(sd$B),
               c("lag1_y_1", "lag1_y_2", "lag2_y_1", "lag2_y_2"))

  # Row 1 of modeled data is t=3. lag1 = t=2 values, lag2 = t=1 values.
  # B is internally centered by grand-mean, so compare against centered values.
  b_cm <- unname(sd$b_center_means)
  expect_equal(unname(sd$B[1, "lag1_y_1"]), 20 - b_cm[1])
  expect_equal(unname(sd$B[1, "lag1_y_2"]), 21 - b_cm[2])
  expect_equal(unname(sd$B[1, "lag2_y_1"]), 10 - b_cm[3])
  expect_equal(unname(sd$B[1, "lag2_y_2"]), 11 - b_cm[4])
})

# ═══════════════════════════════════════════════════════════════════════════════
# §5 Ordinal Y must be integer-coded
# ═══════════════════════════════════════════════════════════════════════════════

test_that("ordinal Y with non-integer values is rejected", {
  df <- data.frame(
    id  = rep(1:2, each = 10),
    t   = rep(1:10, 2),
    y_1 = rep(c(1.5, 2.3, 3.0, 1.7, 2.0), 4)
  )
  expect_error(
    to_stan_data(df, "ordinal", "id", "t", "y_1", character(0), K = 1),
    "integer"
  )
})

test_that("ordinal Y with integer-valued doubles passes", {
  df <- data.frame(
    id  = rep(1:2, each = 10),
    t   = rep(1:10, 2),
    y_1 = rep(c(1.0, 2.0, 3.0, 1.0, 2.0), 4)
  )
  sd <- to_stan_data(df, "ordinal", "id", "t", "y_1", character(0), K = 1)
  expect_true(is.integer(sd$Y))
  expect_equal(sd$C, 3L)
})

# ═══════════════════════════════════════════════════════════════════════════════
# §6 Ordinal with C=2 (minimum categories)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("ordinal C=2 returns correct metadata", {
  df <- data.frame(
    id  = rep(1:3, each = 10),
    t   = rep(1:10, 3),
    y_1 = rep(c(1L, 2L, 1L, 2L, 1L), 6)
  )
  sd <- to_stan_data(df, "ordinal", "id", "t", "y_1", character(0), K = 1)

  expect_equal(sd$C, 2L)
  expect_true("prior_kappa_fam" %in% names(sd))
  expect_true("kappa_scale" %in% names(sd))
})

# ═══════════════════════════════════════════════════════════════════════════════
# §7 Gaussian prior scaling with near-zero SD
# ═══════════════════════════════════════════════════════════════════════════════

test_that("gaussian constant Y keeps priors at unscaled defaults", {
  df <- data.frame(
    id  = rep(1:3, each = 10),
    t   = rep(1:10, 3),
    y_1 = rep(5.0, 30),
    y_2 = rep(3.0, 30)
  )
  priors <- set_priors()
  default_beta_scale  <- priors$beta$scale
  default_sigma_scale <- priors$sigma$scale

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     c("y_1", "y_2"), character(0), K = 1,
                     priors = priors)

  # s_y = 0, so scaling should be skipped; priors stay at defaults
  expect_equal(sd$beta_scale, default_beta_scale)
  expect_equal(sd$sigma_scale, default_sigma_scale)
})

# ═══════════════════════════════════════════════════════════════════════════════
# §8 Empty fe_interactions = list()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("fe_interactions = list() is same as NULL", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 1, family = "gaussian")
  sd_null <- to_stan_data(df, "gaussian", "id", "t",
                          paste0("y_", 1:2), "x_1", K = 1,
                          fe_interactions = NULL)
  sd_empty <- to_stan_data(df, "gaussian", "id", "t",
                           paste0("y_", 1:2), "x_1", K = 1,
                           fe_interactions = list())

  expect_equal(ncol(sd_null$X), ncol(sd_empty$X))
  expect_equal(sd_null$n_fe, sd_empty$n_fe)
})

# ═══════════════════════════════════════════════════════════════════════════════
# §9 Ordinal + re_cols includes "Intercept"
# ═══════════════════════════════════════════════════════════════════════════════

test_that("ordinal with re_cols = 'Intercept' has random intercept in Z", {
  df <- make_test_df(N = 5, T_obs = 20, p = 2, family = "ordinal")
  sd <- to_stan_data(df, "ordinal", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     re_cols = "Intercept")

  # X should NOT contain Intercept (ordinal strips it)
  expect_false("Intercept" %in% colnames(sd$X))
  # Z should contain the random intercept (build_Z reintroduces it)
  expect_equal(sd$n_re, 1L)
  expect_true(all(sd$Z[, 1] == 1))
})

# ═══════════════════════════════════════════════════════════════════════════════
# §10 Multiple time gaps within one subject
# ═══════════════════════════════════════════════════════════════════════════════

test_that("multiple time gaps within one subject are handled correctly", {
  # Times: 1,2,5,6,7,10,11 — gaps at t=5 (from 2) and t=10 (from 7)
  df <- data.frame(
    id  = rep(1L, 7),
    t   = c(1, 2, 5, 6, 7, 10, 11),
    y_1 = c(0, 1, 1, 0, 1, 0, 1)
  )

  # skip_lag = TRUE: keep rows, zero B at gaps

  sd_skip <- to_stan_data(df, "bernoulli", "id", "t",
                          "y_1", character(0), K = 1, skip_lag = TRUE)

  # Modeled rows: t=2,5,6,7,10,11 (6 rows, indices start at K+1=2)
  expect_equal(sd_skip$n_obs, 6L)

  # B is internally centered. Gap rows had raw B = 0, now 0 - b_cm.
  # Non-gap rows had raw B = y_lag, now y_lag - b_cm.
  b_cm <- unname(sd_skip$b_center_means)

  # B for t=5 (gap from t=2, diff=3) should be zeroed (raw 0 → -b_cm)
  # t=5 is the 2nd modeled row
  expect_equal(unname(sd_skip$B[2, 1]), 0 - b_cm[1])

  # B for t=10 (gap from t=7, diff=3) should be zeroed (raw 0 → -b_cm)
  # t=10 is the 5th modeled row
  expect_equal(unname(sd_skip$B[5, 1]), 0 - b_cm[1])

  # B for t=6 (lag from t=5, diff=1) should be non-zero
  # t=6 is the 3rd modeled row; raw B = y_1 at t=5 = 1
  expect_equal(unname(sd_skip$B[3, 1]), 1 - b_cm[1])

  # skip_lag = FALSE: drop gap rows
  sd_drop <- suppressMessages(
    to_stan_data(df, "bernoulli", "id", "t",
                 "y_1", character(0), K = 1, skip_lag = FALSE)
  )
  # Rows with valid lags: t=2, t=6, t=7, t=11 → 4 rows
  expect_equal(sd_drop$n_obs, 4L)
})

# ═══════════════════════════════════════════════════════════════════════════════
# §11 Node-level extraction per family (.to_stan_data_node)
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".to_stan_data_node extracts correct column and applies family rules", {
  # Build shared data with 2 nodes (gaussian-like Y, will cast per family)
  set.seed(99)
  df <- data.frame(
    id  = rep(1:3, each = 15),
    t   = rep(1:15, 3),
    y_1 = rnorm(45, mean = 0, sd = 5),
    y_2 = rnorm(45, mean = 0, sd = 0.1)
  )
  # Make y_1 ordinal-safe (integer 1..3) for ordinal tests below
  df$y_1_ord <- sample(1:3, 45, replace = TRUE)
  # Make y_1 binary for bernoulli tests
  df$y_1_bin <- rbinom(45, 1, 0.5)

  priors <- set_priors()

  # --- Gaussian: scales priors by per-node sd ---
  shared_g <- .to_stan_data_shared(df, "id", "t",
                                    c("y_1", "y_2"), character(0), K = 1)
  node1 <- .to_stan_data_node(shared_g, 1L, "gaussian", set_priors())
  node2 <- .to_stan_data_node(shared_g, 2L, "gaussian", set_priors())

  # Y should be the correct column
  expect_equal(as.numeric(node1$Y), shared_g$Y[, 1])
  expect_equal(as.numeric(node2$Y), shared_g$Y[, 2])

  # Priors should be scaled by per-node sd, not mean sd
  sd_y1 <- sd(shared_g$Y[, 1])
  sd_y2 <- sd(shared_g$Y[, 2])
  expect_equal(node1$beta_scale, 1.0 * sd_y1)  # default scale=1
  expect_equal(node2$beta_scale, 1.0 * sd_y2)
  # They should differ because y_1 has sd≈5, y_2 has sd≈0.1
  expect_true(abs(node1$beta_scale - node2$beta_scale) > 1)

  # X should include Intercept for gaussian
  expect_true("Intercept" %in% colnames(node1$X))

  # --- Bernoulli: keeps intercept, integer Y ---
  shared_b <- .to_stan_data_shared(df, "id", "t",
                                    "y_1_bin", character(0), K = 1)
  node_b <- .to_stan_data_node(shared_b, 1L, "bernoulli", set_priors())
  expect_true(is.integer(node_b$Y))
  expect_true("Intercept" %in% colnames(node_b$X))

  # --- Ordinal: strips intercept from X ---
  shared_o <- .to_stan_data_shared(df, "id", "t",
                                    "y_1_ord", character(0), K = 1)
  node_o <- .to_stan_data_node(shared_o, 1L, "ordinal", set_priors())
  expect_true(is.integer(node_o$Y))
  expect_false("Intercept" %in% colnames(node_o$X))
  expect_true("C" %in% names(node_o))
})
