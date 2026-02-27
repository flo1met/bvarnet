# ──────────────────────────────────────────────────────────────────────────────
# test-to_stan_data.R — data preparation / design matrix tests
#
# Tests the full pipeline from long-format data to Stan-ready data list.
# Covers: dimensions, intercept handling, B matrix (lags), Z matrix (RE),
#         interaction terms, gap handling, centering, and edge cases.
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 Basic dimensions
# ═══════════════════════════════════════════════════════════════════════════════

test_that("stan data has correct p, J, K, n_obs for all families", {
  for (fam in c("bernoulli", "gaussian", "ordinal")) {
    N <- 5; T_obs <- 20; p <- 3; q <- 2; K <- 1

    df <- make_test_df(N = N, T_obs = T_obs, p = p, q = q, family = fam)

    sd <- to_stan_data(
      data = df, family = fam,
      id_col = "id", time_col = "t",
      y_cols = paste0("y_", 1:p), x_cols = paste0("x_", 1:q),
      K = K
    )

    expect_equal(sd$p, p, info = paste("p for", fam))
    expect_equal(sd$J, N, info = paste("J for", fam))
    expect_equal(sd$K, K, info = paste("K for", fam))
    expect_equal(sd$n_obs, N * (T_obs - K), info = paste("n_obs for", fam))
  }
})


test_that("Y matrix has correct dimensions", {
  df <- make_test_df(N = 4, T_obs = 25, p = 3, family = "gaussian")

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:3), character(0), K = 1)

  expect_equal(dim(sd$Y), c(4 * (25 - 1), 3))
  expect_true(all(colnames(sd$Y) == paste0("y_", 1:3)))
})


test_that("B matrix has p*K columns", {
  df <- make_test_df(N = 3, T_obs = 30, p = 4, family = "bernoulli")
  K <- 2

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:4), character(0), K = K)

  expect_equal(ncol(sd$B), 4 * K)  # p * K
  expect_equal(nrow(sd$B), sd$n_obs)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 Intercept handling
# ═══════════════════════════════════════════════════════════════════════════════

test_that("bernoulli and gaussian X matrix includes intercept column", {
  for (fam in c("bernoulli", "gaussian")) {
    df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 1, family = fam)

    sd <- to_stan_data(df, fam, "id", "t",
                       paste0("y_", 1:2), "x_1", K = 1)

    expect_true("Intercept" %in% colnames(sd$X),
                info = paste("Intercept column for", fam))
    expect_equal(sd$n_fe, 2)  # Intercept + x_1

    # Intercept column should be all 1s
    expect_true(all(sd$X[, "Intercept"] == 1),
                info = paste("Intercept values for", fam))
  }
})


test_that("ordinal X matrix does NOT include intercept column", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 1, family = "ordinal")

  sd <- to_stan_data(df, "ordinal", "id", "t",
                     paste0("y_", 1:2), "x_1", K = 1)

  expect_false("Intercept" %in% colnames(sd$X))
  expect_equal(sd$n_fe, 1)  # just x_1
})


test_that("no-covariate model: X is intercept-only (bernoulli/gaussian)", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 0, family = "bernoulli")

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_equal(sd$n_fe, 1)
  expect_equal(ncol(sd$X), 1)
  expect_equal(colnames(sd$X), "Intercept")
})


test_that("no-covariate ordinal: X has zero columns", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 0, family = "ordinal")

  sd <- to_stan_data(df, "ordinal", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_equal(sd$n_fe, 0)
  expect_equal(ncol(sd$X), 0)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 Lag matrix B
# ═══════════════════════════════════════════════════════════════════════════════

test_that("B column names follow lag_y pattern", {
  df <- make_test_df(N = 3, T_obs = 20, p = 3, family = "bernoulli")
  K <- 2

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:3), character(0), K = K)

  expected_names <- c(
    paste0("lag1_y_", 1:3),
    paste0("lag2_y_", 1:3)
  )
  expect_equal(colnames(sd$B), expected_names)
})


test_that("B values match actual lagged Y values (no gaps)", {
  N <- 2; T_obs <- 10; p <- 2; K <- 1
  df <- make_test_df(N = N, T_obs = T_obs, p = p, family = "gaussian")

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:p), character(0), K = K)

  # For the first subject, the first modeled row (t=2) should have
  # B = [y1[t=1], y2[t=1]]
  sub1 <- df[df$id == 1, ]
  expect_equal(unname(sd$B[1, 1]), sub1$y_1[1])
  expect_equal(unname(sd$B[1, 2]), sub1$y_2[1])
})


test_that("K = 1 and K = 2 produce different n_obs per subject", {
  df <- make_test_df(N = 3, T_obs = 20, p = 2, family = "bernoulli")

  sd1 <- to_stan_data(df, "bernoulli", "id", "t",
                      paste0("y_", 1:2), character(0), K = 1)
  sd2 <- to_stan_data(df, "bernoulli", "id", "t",
                      paste0("y_", 1:2), character(0), K = 2)

  expect_equal(sd1$n_obs, 3 * (20 - 1))
  expect_equal(sd2$n_obs, 3 * (20 - 2))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 Random effects matrix Z
# ═══════════════════════════════════════════════════════════════════════════════

test_that("Z is empty (0 columns) with no RE specification", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 1, family = "bernoulli")

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), "x_1", K = 1)

  expect_equal(sd$n_re, 0)
  expect_equal(ncol(sd$Z), 0)
})


test_that("re_cols creates Z with correct columns", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 2, family = "bernoulli")

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), paste0("x_", 1:2),
                     re_cols = "x_1", K = 1)

  expect_equal(sd$n_re, 1)
  expect_equal(ncol(sd$Z), 1)
  # Z values should match the corresponding x_1 values in X
  x1_idx <- which(colnames(sd$X) == "x_1")
  expect_equal(sd$Z[, 1], sd$X[, x1_idx])
})


test_that("re_temporal = TRUE adds p*K lag columns to Z", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 0, family = "bernoulli")
  K <- 1

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0),
                     re_temporal = TRUE, K = K)

  expect_equal(sd$n_re, 2 * K)  # p * K
  expect_equal(ncol(sd$Z), 2 * K)
})


test_that("re_cols + re_temporal combine correctly", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, q = 1, family = "gaussian")
  K <- 1

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), "x_1",
                     re_cols = "x_1", re_temporal = TRUE, K = K)

  # 1 from re_cols + p*K from re_temporal = 1 + 2 = 3
  expect_equal(sd$n_re, 3)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §5 Centering
# ═══════════════════════════════════════════════════════════════════════════════

test_that("center_x = TRUE centers covariate columns but not intercept", {
  df <- make_test_df(N = 5, T_obs = 30, p = 2, q = 2, family = "gaussian")

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), paste0("x_", 1:2),
                     center_x = TRUE, K = 1)

  # Intercept should still be all 1s
  expect_true(all(sd$X[, "Intercept"] == 1))

  # Covariate columns should be approximately centered (mean ≈ 0)
  for (j in 2:ncol(sd$X)) {
    col_mean <- abs(mean(sd$X[, j]))
    expect_lt(col_mean, 0.01)
  }
})


# ═══════════════════════════════════════════════════════════════════════════════
# §6 Time gap handling (skip_lag)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("skip_lag = TRUE zeros out B rows with time gaps", {
  df <- make_test_df(N = 1, T_obs = 10, p = 2, family = "bernoulli")

  # Introduce a gap: remove time point 5
  df <- df[df$t != 5, ]

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     skip_lag = TRUE)

  # Time 6 → lag should be zeroed because t[6] - t[5] != 1 (gap)
  # Row for t=6 is the 5th modeled row (t=2,3,4,6,7,8,9,10 → modeled from t=2)
  # After removing t=5, the sub has 9 rows; first modeled row is t=2 → row 1
  # t=6 comes after removing t=5: raw rows are t=1,2,3,4,6,7,8,9,10
  # modeled rows start from t=2: t=2(1),t=3(2),t=4(3),t=6(4),...
  # Row 4 in the modeled data is t=6 whose lag is t=4 (gap of 2)
  gap_row <- 4  # t=6 in modeled data
  expect_true(all(sd$B[gap_row, ] == 0),
              info = "B row with time gap should be zeroed")
})


test_that("skip_lag = FALSE removes rows with time gaps", {
  df <- make_test_df(N = 1, T_obs = 10, p = 2, family = "bernoulli")

  # Introduce a gap
  df <- df[df$t != 5, ]

  sd_skip <- to_stan_data(df, "bernoulli", "id", "t",
                          paste0("y_", 1:2), character(0), K = 1,
                          skip_lag = TRUE)

  sd_noskip <- to_stan_data(df, "bernoulli", "id", "t",
                            paste0("y_", 1:2), character(0), K = 1,
                            skip_lag = FALSE)

  # skip_lag = FALSE should have fewer rows (gaps removed)
  expect_lt(sd_noskip$n_obs, sd_skip$n_obs)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §7 Subject ID mapping
# ═══════════════════════════════════════════════════════════════════════════════

test_that("id vector maps to sequential integers 1..J", {
  df <- make_test_df(N = 5, T_obs = 15, p = 2, family = "bernoulli")

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_true(all(sd$id %in% 1:5))
  expect_equal(length(unique(sd$id)), 5)
})


test_that("subjects with <= K rows are dropped without error", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "bernoulli")

  # Add a subject with only 1 row (< K+1 = 2)
  short_sub <- data.frame(id = 99, t = 1, y_1 = 0L, y_2 = 1L)
  df <- rbind(df, short_sub)

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  # The short subject should be silently dropped
  expect_equal(sd$J, 4)  # 3 original + 1 short (present in unique IDs)
  # But n_obs should only reflect the 3 contributing subjects
  expect_equal(sd$n_obs, 3 * (15 - 1))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §8 Ordinal-specific fields
# ═══════════════════════════════════════════════════════════════════════════════

test_that("ordinal model adds C and kappa priors to stan data", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "ordinal")

  sd <- to_stan_data(df, "ordinal", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_true("C" %in% names(sd))
  expect_true("prior_kappa_fam" %in% names(sd))
  expect_true(sd$C >= 2)
})


test_that("gaussian model adds sigma priors to stan data", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "gaussian")

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_true("prior_sigma_fam" %in% names(sd))
})


test_that("bernoulli model does NOT add C, sigma, or kappa priors", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "bernoulli")

  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)

  expect_false("C" %in% names(sd))
  expect_false("prior_sigma_fam" %in% names(sd))
  expect_false("prior_kappa_fam" %in% names(sd))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §9 Prior fields always present
# ═══════════════════════════════════════════════════════════════════════════════

test_that("all families include beta, phi, sd prior fields", {
  for (fam in c("bernoulli", "gaussian", "ordinal")) {
    df <- make_test_df(N = 3, T_obs = 15, p = 2, family = fam)

    sd <- to_stan_data(df, fam, "id", "t",
                       paste0("y_", 1:2), character(0), K = 1)

    for (par in c("beta", "phi", "sd")) {
      expect_true(paste0("prior_", par, "_fam") %in% names(sd),
                  info = paste(par, "fam for", fam))
      expect_true(paste0(par, "_loc") %in% names(sd),
                  info = paste(par, "loc for", fam))
      expect_true(paste0(par, "_scale") %in% names(sd),
                  info = paste(par, "scale for", fam))
    }
  }
})


# ═══════════════════════════════════════════════════════════════════════════════
# §10 Interaction terms
# ═══════════════════════════════════════════════════════════════════════════════

test_that("fe_interactions adds columns to X", {
  df <- make_test_df(N = 3, T_obs = 20, p = 2, q = 2, family = "bernoulli")

  sd_no_int <- to_stan_data(df, "bernoulli", "id", "t",
                            paste0("y_", 1:2), paste0("x_", 1:2), K = 1)
  sd_int <- to_stan_data(df, "bernoulli", "id", "t",
                         paste0("y_", 1:2), paste0("x_", 1:2),
                         fe_interactions = list(c("x_1", "x_2")), K = 1)

  expect_gt(sd_int$n_fe, sd_no_int$n_fe)
  expect_true("x_1:x_2" %in% colnames(sd_int$X))
})


test_that("lag interaction adds p*K columns to X", {
  df <- make_test_df(N = 3, T_obs = 20, p = 2, q = 1, family = "bernoulli")
  K <- 1

  sd_int <- to_stan_data(df, "bernoulli", "id", "t",
                         paste0("y_", 1:2), "x_1",
                         fe_interactions = list(c("lag", "x_1")), K = K)

  # lag interaction should add p*K = 2 columns
  expected_cols <- paste0("lag1_y_", 1:2, ":x_1")
  expect_true(all(expected_cols %in% colnames(sd_int$X)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §11 Consistency: sim_var → to_stan_data round-trip
# ═══════════════════════════════════════════════════════════════════════════════

test_that("to_stan_data on sim_var output produces valid stan data (bernoulli)", {
  sd <- get_standata_bernoulli()

  expect_true(all(sd$Y %in% c(0L, 1L)))
  expect_equal(sd$p, 3)
  expect_equal(sd$K, 1)
  expect_true(sd$n_obs > 0)
  expect_true(sd$n_fe >= 1)  # at least intercept
})


test_that("to_stan_data on sim_var output produces valid stan data (gaussian)", {
  sd <- get_standata_gaussian()

  expect_true(all(is.finite(sd$Y)))
  expect_equal(sd$p, 3)
})


test_that("to_stan_data on sim_var output produces valid stan data (ordinal)", {
  sd <- get_standata_ordinal()

  expect_true(all(sd$Y >= 1))
  expect_true(all(sd$Y <= sd$C))
  expect_equal(sd$p, 3)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §12 Listwise deletion
# ═══════════════════════════════════════════════════════════════════════════════

test_that("listwise deletion removes rows with NA in y or x columns", {
  df <- make_test_df(N = 3, T_obs = 20, p = 2, q = 1, family = "gaussian")

  # Inject some NAs
  df$y_1[5]  <- NA
  df$x_1[10] <- NA

  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), "x_1", K = 1)

  # Should have fewer obs than clean data (some rows deleted)
  sd_clean <- to_stan_data(make_test_df(N = 3, T_obs = 20, p = 2, q = 1, family = "gaussian"),
                           "gaussian", "id", "t",
                           paste0("y_", 1:2), "x_1", K = 1)

  expect_lt(sd$n_obs, sd_clean$n_obs)
  # No NA in Y or X
  expect_false(any(is.na(sd$Y)))
  expect_false(any(is.na(sd$X)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §13 Prior injection
# ═══════════════════════════════════════════════════════════════════════════════

test_that("to_stan_data runs without error using default set_priors()", {
  for (fam in c("bernoulli", "gaussian", "ordinal")) {
    df <- make_test_df(N = 3, T_obs = 15, p = 2, family = fam)
    expect_no_error(
      to_stan_data(df, fam, "id", "t",
                   paste0("y_", 1:2), character(0), K = 1,
                   priors = set_priors())
    )
  }
})

test_that("default priors produce prior_beta_fam == 1L", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2)
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     priors = set_priors())
  expect_equal(sd$prior_beta_fam, 1L)
})

test_that("default priors produce beta_scale == 1 for bernoulli", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2)
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     priors = set_priors())
  expect_equal(sd$beta_scale, 1)
})

test_that("custom cauchy beta prior is injected correctly", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2)
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     priors = set_priors(beta = prior("cauchy", 0, 0.5)))
  expect_equal(sd$prior_beta_fam, 3L)
  expect_equal(sd$beta_scale, 0.5)
})

test_that("custom student_t phi prior is injected with correct df", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2)
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     priors = set_priors(phi = prior("student_t", 0, 0.3, df = 5)))
  expect_equal(sd$prior_phi_fam, 2L)
  expect_equal(sd$phi_df, 5)
})

test_that("gaussian model output includes prior_sigma_fam and sigma_scale", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "gaussian")
  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)
  expect_true("prior_sigma_fam" %in% names(sd))
  expect_true("sigma_scale"     %in% names(sd))
})

test_that("gaussian model with default priors scales beta_scale by s_y", {
  set.seed(1)
  df <- make_test_df(N = 5, T_obs = 20, p = 2, family = "gaussian", seed = 1)
  # multiply Y values so sd >> 1
  df$y_1 <- df$y_1 * 10
  df$y_2 <- df$y_2 * 10
  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     priors = set_priors())
  expect_gt(sd$beta_scale, 1)   # scaled up because sd(Y) >> 1
})

test_that("gaussian model with user-specified beta prior is NOT scaled", {
  set.seed(2)
  df <- make_test_df(N = 5, T_obs = 20, p = 2, family = "gaussian")
  df$y_1 <- df$y_1 * 10
  df$y_2 <- df$y_2 * 10
  sd <- to_stan_data(df, "gaussian", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1,
                     priors = set_priors(beta = prior("normal", 0, 0.5)))
  # User set scale = 0.5; must NOT be multiplied by s_y
  expect_equal(sd$beta_scale, 0.5)
})

test_that("ordinal model includes prior_kappa_fam", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "ordinal")
  sd <- to_stan_data(df, "ordinal", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)
  expect_true("prior_kappa_fam" %in% names(sd))
})

test_that("bernoulli model does NOT include prior_kappa_fam", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2, family = "bernoulli")
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)
  expect_false("prior_kappa_fam" %in% names(sd))
})

test_that("to_stan_data without explicit priors argument uses set_priors() defaults", {
  df <- make_test_df(N = 3, T_obs = 15, p = 2)
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), character(0), K = 1)
  expect_equal(sd$prior_beta_fam, 1L)
  expect_equal(sd$beta_scale, 1)
})

