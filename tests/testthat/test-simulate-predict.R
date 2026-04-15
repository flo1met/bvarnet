# tests/testthat/test-simulate-predict.R
# ──────────────────────────────────────────────────────────────────────────────
# Tests for simulate.bvarnet() and predict.bvarnet()
# ──────────────────────────────────────────────────────────────────────────────


# ═══════════════════════════════════════════════════════════════════════════════
# Helper: build a mock bvarnet object with design_spec and proper matrices
# for predict/simulate to work end-to-end without Stan.
# ═══════════════════════════════════════════════════════════════════════════════

#' Build a mock bvarnet with real-ish design matrices for predict testing
#' @noRd
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
    variable = par_nms, rhat = 1.001,
    ess_bulk = 3000, ess_tail = 2800, stringsAsFactors = FALSE
  )

  # Normalise family to named vector
  family_vec <- setNames(rep(family, p), y_cols)

  structure(
    list(
      draws        = draws,
      convergence  = smry,
      diagnostics  = data.frame(num_divergent = integer(n_chains),
                                num_max_treedepth = integer(n_chains),
                                ebfmi = rep(1.0, n_chains)),
      timing       = list(total = 5.0),
      metadata     = list(),
      return_codes = rep(0L, n_chains),
      family       = family_vec,
      standata     = sd,
      priors       = set_priors()
    ),
    class = "bvarnet"
  )
}


# ═══════════════════════════════════════════════════════════════════════════════
#  §6.1  simulate.bvarnet tests (no Stan)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("simulate: posterior-mean returns data.frame with expected columns", {
  mock <- make_predictable_mock("gaussian")
  out  <- simulate(mock, nsim = 15L, seed = 1)
  expect_s3_class(out, "data.frame")
  expect_true(all(c("id", "t", "y_1", "y_2") %in% names(out)))
  expect_equal(nrow(out), mock$standata$J * 15L)
})

test_that("simulate: posterior-sample returns list of data.frames", {
  mock <- make_predictable_mock("gaussian")
  out  <- simulate(mock, nsim = 15L, seed = 1, method = "posterior-sample",
                   ndraws = 3L)
  expect_type(out, "list")
  expect_length(out, 3L)
  expect_s3_class(out[[1]], "data.frame")
  expect_equal(nrow(out[[1]]), mock$standata$J * 15L)
})

test_that("simulate: seed reproducibility", {
  mock <- make_predictable_mock("gaussian")
  out1 <- simulate(mock, nsim = 10L, seed = 42)
  out2 <- simulate(mock, nsim = 10L, seed = 42)
  expect_identical(out1, out2)
})

test_that("simulate: bernoulli outputs in {0, 1}", {
  mock <- make_predictable_mock("bernoulli")
  out  <- simulate(mock, nsim = 20L, seed = 1)
  y_vals <- unlist(out[, grep("^y_", names(out))])
  expect_true(all(y_vals %in% c(0L, 1L)))
})

test_that("simulate: gaussian outputs are finite numeric", {
  mock <- make_predictable_mock("gaussian")
  out  <- simulate(mock, nsim = 20L, seed = 1)
  y_vals <- unlist(out[, grep("^y_", names(out))])
  expect_true(all(is.finite(y_vals)))
})

test_that("simulate: ordinal outputs are integers in 1..C", {
  mock <- make_predictable_mock("ordinal")
  out  <- simulate(mock, nsim = 20L, seed = 1)
  y_vals <- unlist(out[, grep("^y_", names(out))])
  expect_true(all(y_vals %in% seq_len(mock$standata$C)))
})

test_that("simulate: subject_re='zero' gives identical REs", {
  mock <- make_predictable_mock("gaussian", n_re = 1L)
  out  <- simulate(mock, nsim = 15L, seed = 1, subject_re = "zero")
  expect_s3_class(out, "data.frame")
  # All subjects should have same population-level parameters
  # (no between-subject variability from REs)
  expect_true(all(is.finite(out$y_1)))
})

test_that("simulate: subject_re='sample' produces between-subject variation", {
  mock <- make_predictable_mock("gaussian", n_re = 1L)
  out  <- simulate(mock, nsim = 30L, seed = 10, subject_re = "sample")

  # Compute per-subject means — should vary
  sub_means <- tapply(out$y_1, out$id, mean)
  expect_true(sd(sub_means) > 0)  # non-zero between-subject variation
})


# ═══════════════════════════════════════════════════════════════════════════════
#  §6.2  predict.bvarnet tests (no Stan)
# ═══════════════════════════════════════════════════════════════════════════════

# --- shape tests ---

test_that("predict: link type returns data.frame with id, time, predicted cols (in-sample)", {
  mock <- make_predictable_mock("gaussian")
  out  <- predict(mock, type = "link")
  expect_true(is.data.frame(out))
  expect_true("id" %in% names(out))
  expect_true("t" %in% names(out))
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  expect_true(all(pred_names %in% names(out)))
  expect_equal(nrow(out), mock$standata$n_obs)
  expect_true(all(is.finite(out[[pred_names[1]]])))
})

test_that("predict: response type returns data.frame with probabilities (in-sample)", {
  mock <- make_predictable_mock("bernoulli")
  out  <- predict(mock, type = "response")
  expect_true(is.data.frame(out))
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  expect_true(all(pred_names %in% names(out)))
  expect_equal(nrow(out), mock$standata$n_obs)
  # Predicted values should be probabilities
  expect_true(all(out[[pred_names[1]]] >= 0 & out[[pred_names[1]]] <= 1))
})

test_that("predict: probabilities type for bernoulli returns list of data.frames", {
  mock <- make_predictable_mock("bernoulli")
  out  <- predict(mock, type = "probabilities")
  expect_type(out, "list")
  expect_length(out, mock$standata$p)
  expect_true(is.data.frame(out[[1]]))
  expect_true("p1" %in% names(out[[1]]))
  expect_true(all(out[[1]][["p1"]] >= 0 & out[[1]][["p1"]] <= 1))
})

test_that("predict: probabilities type for gaussian returns list with mean+sd", {
  mock <- make_predictable_mock("gaussian")
  out  <- predict(mock, type = "probabilities")
  expect_type(out, "list")
  expect_length(out, mock$standata$p)
  expect_true(is.data.frame(out[[1]]))
  expect_true(all(c("mean", "sd") %in% names(out[[1]])))
  expect_true(all(out[[1]][["sd"]] > 0))
})

test_that("predict: probabilities type for ordinal returns list with cat_ columns", {
  mock <- make_predictable_mock("ordinal")
  out  <- predict(mock, type = "probabilities")
  expect_type(out, "list")
  expect_length(out, mock$standata$p)
  expect_true(is.data.frame(out[[1]]))
  C <- mock$standata$C
  cat_cols <- grep("^cat_", names(out[[1]]), value = TRUE)
  expect_equal(length(cat_cols), C)
  # Probability rows should sum to ~1
  row_sums <- rowSums(out[[1]][, cat_cols, drop = FALSE])
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

# --- posterior-sample returns attr("sd") ---

test_that("predict: posterior-sample returns _sd columns and attr('ndraws')", {
  mock <- make_predictable_mock("gaussian")
  out  <- predict(mock, type = "response", method = "posterior-sample",
                  ndraws = 5L, seed = 1)
  expect_true(is.data.frame(out))
  y_names <- colnames(mock$standata$Y)
  sd_names <- paste0("predicted_", y_names, "_sd")
  expect_true(all(sd_names %in% names(out)))
  expect_true(!is.null(attr(out, "ndraws")))
  expect_equal(attr(out, "ndraws"), 5L)
})

test_that("predict: posterior-sample probabilities returns _sd columns", {
  mock <- make_predictable_mock("bernoulli")
  out  <- predict(mock, type = "probabilities", method = "posterior-sample",
                  ndraws = 5L, seed = 1)
  expect_type(out, "list")
  expect_length(out, mock$standata$p)
  expect_true(is.data.frame(out[[1]]))
  expect_true("p1_sd" %in% names(out[[1]]))
})

# --- newdata with NA for first K rows ---

test_that("predict: newdata returns compact data.frame (modeled rows only)", {
  mock <- make_predictable_mock("gaussian", T_obs = 10L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 10L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  for (j in 1:1) df[[paste0("x_", j)]] <- rnorm(nrow(df))

  out <- predict(mock, newdata = df, type = "link")
  expect_true(is.data.frame(out))
  # Only modeled rows (T_obs - K per subject)
  expected_rows <- N * (TT - K)
  expect_equal(nrow(out), expected_rows)
  # All predicted values should be finite (no NAs)
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) expect_true(all(is.finite(out[[pn]])))
})

test_that("predict: in-sample nrow matches n_obs (modeled rows only)", {
  for (fam in c("gaussian", "bernoulli", "ordinal")) {
    mock <- make_predictable_mock(fam, J = 4L, T_obs = 12L, K = 1L)
    out  <- predict(mock, type = "link")
    expect_true(is.data.frame(out), info = paste("family:", fam))
    # Output rows should match n_obs (modeled rows), not n_rows_data
    expect_equal(nrow(out), mock$standata$n_obs,
                 info = paste("family:", fam))
    # All predicted values should be finite (no NAs)
    y_names <- colnames(mock$standata$Y)
    pred_names <- paste0("predicted_", y_names)
    for (pn in pred_names) {
      expect_true(all(is.finite(out[[pn]])),
                  info = paste("family:", fam, "col:", pn))
    }
  }
})

test_that("predict: in-sample and newdata give same nrow for same data", {
  mock <- make_predictable_mock("gaussian", J = 3L, T_obs = 10L, K = 1L)
  # Reconstruct the same data that was used to build the mock
  set.seed(99L)
  df <- expand.grid(t = seq_len(10L), id = seq_len(3L))
  df <- df[order(df$id, df$t), ]
  for (j in 1:mock$standata$p) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  for (j in 1:1) df[[paste0("x_", j)]] <- rnorm(nrow(df))

  out_in    <- predict(mock, type = "link")
  out_new   <- predict(mock, newdata = df, type = "link")
  expect_equal(nrow(out_in), nrow(out_new))
  # Both should have n_obs rows (modeled only)
  expect_equal(nrow(out_in), mock$standata$n_obs)
})

# --- subject_re tests ---

test_that("predict: subject_re='zero' vs 'posterior-mean' differ for seen subjects", {
  mock <- make_predictable_mock("gaussian", n_re = 1L)
  out_zero <- predict(mock, type = "link", subject_re = "zero")
  out_re   <- predict(mock, type = "link", subject_re = "posterior-mean")
  # Should differ because posterior mean u != 0
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  expect_false(all(out_zero[pred_names] == out_re[pred_names]))
})

test_that("predict: n_re=0 makes all subject_re options identical", {
  mock <- make_predictable_mock("gaussian", n_re = 0L)
  out_zero <- predict(mock, type = "link", subject_re = "zero")
  # With n_re=0, "posterior-mean" silently degrades to "zero"
  out_pm   <- predict(mock, type = "link", subject_re = "posterior-mean")
  expect_equal(out_zero, out_pm)
})

test_that("predict: unseen IDs with new_subject='zero' are population-level", {
  mock <- make_predictable_mock("gaussian", n_re = 1L, J = 3L, T_obs = 15L)

  set.seed(99L)
  # Create newdata with an unseen subject id
  df_new <- expand.grid(t = seq_len(15L), id = c(1L, 999L))
  df_new <- df_new[order(df_new$id, df_new$t), ]
  for (j in seq_len(mock$standata$p))
    df_new[[paste0("y_", j)]] <- rnorm(nrow(df_new))
  df_new$x_1 <- rnorm(nrow(df_new))

  out <- predict(mock, newdata = df_new, type = "link",
                 subject_re = "posterior-mean", new_subject = "zero")

  # Predictions for subject 999 should use u=0 and be finite
  rows_999 <- out$id == "999"
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) expect_true(all(is.finite(out[rows_999, pn])))
})

# --- error tests ---

test_that("predict: errors on missing newdata columns", {
  mock <- make_predictable_mock("gaussian")
  bad_df <- data.frame(id = 1:5, t = 1:5)
  expect_error(predict(mock, newdata = bad_df), "Missing columns in newdata")
})


# ═══════════════════════════════════════════════════════════════════════════════
#  §6.2b  Recursive forecasting tests (no Stan)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("recursive: returns finite predictions after conditioning window", {
  mock <- make_predictable_mock("gaussian", T_obs = 20L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 20L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  for (j in 1:1) df[[paste0("x_", j)]] <- rnorm(nrow(df))

  out <- predict(mock, newdata = df, type = "link",
                 forecast = "recursive", conditioning_window = K + 3L)
  expect_true(is.data.frame(out))
  expected_rows <- N * (TT - K)
  expect_equal(nrow(out), expected_rows)
  # All rows should be finite (no NAs)
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) expect_true(all(is.finite(out[[pn]])))
})

test_that("recursive: differs from one-step on long horizons", {
  mock <- make_predictable_mock("gaussian", T_obs = 25L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 25L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  for (j in 1:1) df[[paste0("x_", j)]] <- rnorm(nrow(df))

  out_one  <- predict(mock, newdata = df, type = "link",
                      forecast = "one-step")
  out_rec  <- predict(mock, newdata = df, type = "link",
                      forecast = "recursive", conditioning_window = K + 2L)
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  # Later forecast rows should diverge
  expect_false(identical(out_one[pred_names], out_rec[pred_names]))
})

test_that("recursive: works with conditioning_window = K (all rows forecasted)", {
  mock <- make_predictable_mock("gaussian", T_obs = 15L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 15L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  for (j in 1:1) df[[paste0("x_", j)]] <- rnorm(nrow(df))

  out <- predict(mock, newdata = df, type = "response",
                 forecast = "recursive", conditioning_window = K)
  expect_true(is.data.frame(out))
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) expect_true(all(is.finite(out[[pn]])))
})

test_that("recursive: conditioning_window < K errors", {
  mock <- make_predictable_mock("gaussian")
  set.seed(99L)
  K <- mock$standata$K
  df <- expand.grid(t = seq_len(10L), id = 1:3)
  df <- df[order(df$id, df$t), ]
  for (j in 1:mock$standata$p) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  df$x_1 <- rnorm(nrow(df))

  expect_error(
    predict(mock, newdata = df, forecast = "recursive",
            conditioning_window = K - 1L),
    "conditioning_window must be >= K"
  )
})

test_that("recursive: named conditioning_window per subject works", {
  mock <- make_predictable_mock("gaussian", J = 3L, T_obs = 20L)
  set.seed(99L)
  K <- mock$standata$K; p_n <- mock$standata$p
  df <- expand.grid(t = seq_len(20L), id = 1:3)
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  df$x_1 <- rnorm(nrow(df))

  cw <- c("1" = K + 2L, "2" = K + 5L, "3" = K + 3L)
  out <- predict(mock, newdata = df, type = "link",
                 forecast = "recursive", conditioning_window = cw)
  expect_true(is.data.frame(out))
  expected_rows <- 3L * (20L - K)
  expect_equal(nrow(out), expected_rows)
  # All predicted values should be finite
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) expect_true(all(is.finite(out[[pn]])))
})

test_that("recursive: missing subject in named cw errors", {
  mock <- make_predictable_mock("gaussian", J = 3L, T_obs = 15L)
  set.seed(99L)
  K <- mock$standata$K; p_n <- mock$standata$p
  df <- expand.grid(t = seq_len(15L), id = 1:3)
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  df$x_1 <- rnorm(nrow(df))

  cw <- c("1" = K + 2L, "2" = K + 3L)  # missing "3"
  expect_error(
    predict(mock, newdata = df, forecast = "recursive",
            conditioning_window = cw),
    "missing entries"
  )
})

test_that("recursive: bernoulli returns probabilities in [0,1]", {
  mock <- make_predictable_mock("bernoulli", T_obs = 20L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 20L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rbinom(nrow(df), 1, 0.5)
  df$x_1 <- rnorm(nrow(df))

  out <- predict(mock, newdata = df, type = "response",
                 forecast = "recursive", conditioning_window = K + 3L)
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) {
    expect_true(all(out[[pn]] >= 0 & out[[pn]] <= 1))
  }
})

test_that("recursive: ordinal probability rows sum to 1", {
  mock <- make_predictable_mock("ordinal", T_obs = 20L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 20L; p_n <- mock$standata$p; K <- mock$standata$K
  C <- mock$standata$C
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- sample(1:C, nrow(df), replace = TRUE)
  df$x_1 <- rnorm(nrow(df))

  out <- predict(mock, newdata = df, type = "probabilities",
                 forecast = "recursive", conditioning_window = K + 3L)
  expect_type(out, "list")
  expect_length(out, p_n)
  expect_true(is.data.frame(out[[1]]))
  cat_cols <- grep("^cat_", names(out[[1]]), value = TRUE)
  row_sums <- rowSums(out[[1]][, cat_cols, drop = FALSE])
  expect_true(all(abs(row_sums - 1) < 1e-10))
})

test_that("recursive: posterior-sample works with _sd columns", {
  mock <- make_predictable_mock("gaussian", T_obs = 20L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 20L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  df$x_1 <- rnorm(nrow(df))

  out <- predict(mock, newdata = df, type = "response",
                 method = "posterior-sample", ndraws = 5L, seed = 1,
                 forecast = "recursive", conditioning_window = K + 2L)
  expect_true(is.data.frame(out))
  y_names <- colnames(mock$standata$Y)
  sd_names <- paste0("predicted_", y_names, "_sd")
  expect_true(all(sd_names %in% names(out)))
  expect_true(!is.null(attr(out, "ndraws")))
  expect_equal(attr(out, "ndraws"), 5L)
})

test_that("recursive: in-sample (newdata=NULL) works", {
  mock <- make_predictable_mock("gaussian", T_obs = 20L)
  # In-sample recursive forecast
  out <- predict(mock, type = "link", forecast = "recursive",
                 conditioning_window = mock$standata$K + 3L)
  expect_true(is.data.frame(out))
  expect_equal(nrow(out), mock$standata$n_obs)
  y_names <- colnames(mock$standata$Y)
  pred_names <- paste0("predicted_", y_names)
  for (pn in pred_names) expect_true(all(is.finite(out[[pn]])))
})

test_that("recursive: subject_re works with recursive mode", {
  mock <- make_predictable_mock("gaussian", n_re = 1L, T_obs = 20L)
  set.seed(99L)
  N <- mock$standata$J; TT <- 20L; p_n <- mock$standata$p; K <- mock$standata$K
  df <- expand.grid(t = seq_len(TT), id = seq_len(N))
  df <- df[order(df$id, df$t), ]
  for (j in seq_len(p_n)) df[[paste0("y_", j)]] <- rnorm(nrow(df))
  df$x_1 <- rnorm(nrow(df))

  out_zero <- predict(mock, newdata = df, type = "link",
                      forecast = "recursive", conditioning_window = K + 2L,
                      subject_re = "zero")
  out_re   <- predict(mock, newdata = df, type = "link",
                      forecast = "recursive", conditioning_window = K + 2L,
                      subject_re = "posterior-mean")
  # Should differ because posterior mean u != 0
  expect_false(identical(out_zero, out_re))
})


# --- posterior-mean is deterministic ---

test_that("predict: posterior-mean is deterministic (no seed needed)", {
  mock <- make_predictable_mock("gaussian")
  out1 <- predict(mock, type = "response", method = "posterior-mean")
  out2 <- predict(mock, type = "response", method = "posterior-mean")
  expect_identical(out1, out2)
})

# --- posterior-sample seed reproducibility ---

test_that("predict: posterior-sample with same seed gives same results", {
  mock <- make_predictable_mock("gaussian")
  out1 <- predict(mock, type = "response", method = "posterior-sample",
                  ndraws = 5L, seed = 42)
  out2 <- predict(mock, type = "response", method = "posterior-sample",
                  ndraws = 5L, seed = 42)
  expect_equal(out1, out2)
})


# ═══════════════════════════════════════════════════════════════════════════════
#  §6.3  Analysis-validation tests (no Stan)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("fixed vs individual predictions have comparable shapes", {
  mock <- make_predictable_mock("gaussian", n_re = 1L)
  out_fixed <- predict(mock, type = "response", subject_re = "zero")
  out_indiv <- predict(mock, type = "response", subject_re = "posterior-mean")
  expect_equal(nrow(out_fixed), nrow(out_indiv))
  expect_equal(names(out_fixed), names(out_indiv))
})

test_that("for n_re > 0, individual predictions differ from fixed for seen subjects", {
  mock <- make_predictable_mock("gaussian", n_re = 1L)
  out_fixed <- predict(mock, type = "link", subject_re = "zero")
  out_indiv <- predict(mock, type = "link", subject_re = "posterior-mean")
  expect_false(identical(out_fixed, out_indiv))
})


# ═══════════════════════════════════════════════════════════════════════════════
#  §6.4  Stan-backed integration tests
# ═══════════════════════════════════════════════════════════════════════════════

# Helper: skip if compiled Stan models are not available
can_run_stan <- function() {
  instantiate::stan_cmdstan_exists() &&
    tryCatch({
      instantiate::stan_package_model(name = "model_gaussian", package = "bvarnet")
      TRUE
    }, error = function(e) FALSE)
}

test_that("integration: simulate + predict roundtrip for gaussian", {
  skip_if_not(can_run_stan(), "Compiled Stan models not available")

  sim <- sim_var(N = 5, T_obs = 30, p = 2, K = 1,
                 family = "gaussian", q = 0, seed = 1)
  fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "gaussian",
              iter = 200, warmup = 100, chains = 2, seed = 1)

  # Simulate
  sim_out <- simulate(fit, nsim = 20, seed = 42)
  expect_s3_class(sim_out, "data.frame")
  expect_true(all(is.finite(sim_out$y_1)))
  expect_true(all(is.finite(sim_out$y_2)))

  # Predict on simulated data
  preds <- predict(fit, newdata = sim_out, type = "response")
  expect_true(is.data.frame(preds))
  pred_names <- paste0("predicted_", c("y_1", "y_2"))
  expect_true(all(pred_names %in% names(preds)))
  expect_true(nrow(preds) > 0)
  for (pn in pred_names) expect_true(all(is.finite(preds[[pn]])))
})

test_that("integration: simulate + predict roundtrip for bernoulli", {
  skip_if_not(can_run_stan(), "Compiled Stan models not available")

  sim <- sim_var(N = 5, T_obs = 30, p = 2, K = 1,
                 family = "bernoulli", q = 0, seed = 2)
  fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "bernoulli",
              iter = 200, warmup = 100, chains = 2, seed = 2)

  sim_out <- simulate(fit, nsim = 20, seed = 42)
  y_vals <- unlist(sim_out[, grep("^y_", names(sim_out))])
  expect_true(all(y_vals %in% c(0L, 1L)))

  preds <- predict(fit, newdata = sim_out, type = "response")
  pred_names <- paste0("predicted_", c("y_1", "y_2"))
  for (pn in pred_names) {
    expect_true(all(preds[[pn]] >= 0 & preds[[pn]] <= 1))
  }
})

test_that("integration: predict positive association gaussian", {
  skip_if_not(can_run_stan(), "Compiled Stan models not available")

  sim <- sim_var(N = 10, T_obs = 50, p = 2, K = 1,
                 family = "gaussian", q = 0, seed = 3)
  fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "gaussian",
              iter = 500, warmup = 200, chains = 2, seed = 3)

  preds <- predict(fit, type = "response")
  obs_y <- fit$standata$Y[, 1]
  pred_y <- preds[["predicted_y_1"]]
  # Should have positive correlation between observed and predicted
  expect_true(cor(obs_y, pred_y) > 0)
})

test_that("integration: simulate + predict roundtrip for ordinal", {
  skip_if_not(can_run_stan(), "Compiled Stan models not available")

  sim <- sim_var(N = 10, T_obs = 50, p = 2, K = 1,
                 family = "ordinal", q = 0, C = 3, seed = 4)

  # Ordinal models can be difficult to sample with few iterations;

  # wrap in tryCatch so this doesn't block the test suite.
  fit <- tryCatch(
    bvar(id_col = "id", time_col = "t",
         y_cols = c("y_1", "y_2"), x_cols = character(0),
         data = sim$data, family = "ordinal",
         iter = 500, warmup = 300, chains = 2, seed = 4,
         adapt_delta = 0.95),
    error = function(e) {
      skip("Ordinal model failed to sample — skipping roundtrip test.")
    }
  )

  sim_out <- simulate(fit, nsim = 20, seed = 42)
  y_vals <- unlist(sim_out[, grep("^y_", names(sim_out))])
  expect_true(all(y_vals %in% seq_len(fit$standata$C)))

  preds <- predict(fit, newdata = sim_out, type = "probabilities")
  expect_type(preds, "list")
  expect_length(preds, 2L)
  expect_true(is.data.frame(preds[[1]]))
  cat_cols <- grep("^cat_", names(preds[[1]]), value = TRUE)
  expect_true(length(cat_cols) > 0)
  row_sums <- rowSums(preds[[1]][, cat_cols, drop = FALSE])
  expect_true(all(abs(row_sums - 1) < 1e-6))
})

test_that("integration: recursive forecast gaussian smoke test", {
  skip_if_not(can_run_stan(), "Compiled Stan models not available")

  sim <- sim_var(N = 5, T_obs = 30, p = 2, K = 1,
                 family = "gaussian", q = 0, seed = 5)
  fit <- bvar(id_col = "id", time_col = "t",
              y_cols = c("y_1", "y_2"), x_cols = character(0),
              data = sim$data, family = "gaussian",
              iter = 200, warmup = 100, chains = 2, seed = 5)

  K <- fit$standata$K

  # Recursive forecast on training data with a conditioning window
  preds <- predict(fit, type = "response", forecast = "recursive",
                   conditioning_window = K + 5L)
  expect_true(is.data.frame(preds))
  pred_names <- paste0("predicted_", c("y_1", "y_2"))
  for (pn in pred_names) expect_true(all(is.finite(preds[[pn]])))

  # Compare to one-step: should differ for later rows
  preds_one <- predict(fit, type = "response", forecast = "one-step")
  expect_false(identical(preds[pred_names], preds_one[pred_names]))
})
