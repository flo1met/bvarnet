# ──────────────────────────────────────────────────────────────────────────────
# test-summary-extractors.R — tests for summary.bvarnet(), extract_temporal(),
# extract_random_effects(), extract_network_matrix(), and exported extract_draws()
# ──────────────────────────────────────────────────────────────────────────────


# ═══════════════════════════════════════════════════════════════════════════════
# §1 — summary.bvarnet()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("summary.bvarnet returns summary.bvarnet class", {
  obj <- make_mock_bvarnet("gaussian")
  s   <- summary(obj)
  expect_s3_class(s, "summary.bvarnet")
})

test_that("summary.bvarnet contains expected elements", {
  obj <- make_mock_bvarnet("bernoulli")
  s   <- summary(obj)
  expect_true(all(c("table", "family", "p", "K", "n",
                     "rhat_max", "n_divergences") %in% names(s)))
})

test_that("summary.bvarnet table is a data.frame with correct columns", {
  obj <- make_mock_bvarnet("gaussian")
  s   <- summary(obj)
  expect_true(is.data.frame(s$table))
  expected_cols <- c("type", "predictor", "outcome",
                     "mean", "median", "q5", "q95",
                     "rhat", "ess_bulk", "ess_tail")
  expect_true(all(expected_cols %in% names(s$table)))
})

test_that("summary.bvarnet prints without error", {
  obj <- make_mock_bvarnet("gaussian")
  s   <- summary(obj)
  expect_output(print(s), "BVAR Network Summary")
  expect_output(print(s), "gaussian")
  expect_output(print(s), "Autoregressive")
})

test_that("summary.bvarnet warns for high Rhat", {
  obj <- make_mock_bvarnet()
  obj$convergence$rhat <- rep(1.05, nrow(obj$convergence))
  s <- summary(obj)
  expect_output(print(s), "WARNING")
})

test_that("summary.bvarnet with bayes_factor = TRUE includes BF10", {
  obj <- make_mock_bvarnet("gaussian")
  s   <- summary(obj, bayes_factor = TRUE)
  expect_true("BF10" %in% names(s$table))
})

test_that("print.summary.bvarnet returns object invisibly", {
  obj <- make_mock_bvarnet()
  s   <- summary(obj)
  ret <- withVisible(print(s))
  expect_false(ret$visible)
  expect_identical(ret$value, s)
})

test_that("print.summary.bvarnet truncates groups and shows footer", {
  obj <- make_mock_bvarnet("gaussian")
  s   <- summary(obj)
  # max_rows=1 should truncate almost every group
  out <- capture.output(print(s, max_rows = 1))
  txt <- paste(out, collapse = "\n")
  expect_true(grepl("more rows", txt))
  expect_true(grepl("extract_param\\(\\)", txt))
  expect_true(grepl("extract_temporal\\(\\)", txt))
  expect_true(grepl("extract_network_matrix\\(\\)", txt))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 — extract_temporal()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_temporal returns only temporal types", {
  obj <- make_mock_bvarnet("bernoulli")
  out <- extract_temporal(obj)
  expect_true(is.data.frame(out))
  expect_true(all(out$type %in% c("Autoregressive", "Cross-lagged")))
})

test_that("extract_temporal effect='ar' returns only AR rows", {
  obj <- make_mock_bvarnet("bernoulli")
  out <- extract_temporal(obj, effect = "ar")
  expect_true(all(out$type == "Autoregressive"))
  expect_true(nrow(out) > 0)
})

test_that("extract_temporal effect='cl' returns only CL rows", {
  obj <- make_mock_bvarnet("bernoulli")
  out <- extract_temporal(obj, effect = "cl")
  expect_true(all(out$type == "Cross-lagged"))
  expect_true(nrow(out) > 0)
})

test_that("extract_temporal lag filter works", {
  obj <- make_mock_bvarnet("bernoulli")  # K=1
  out <- extract_temporal(obj, lag = 1)
  expect_true(nrow(out) > 0)
  expect_true(all(grepl("^lag1_", out$predictor)))
})

test_that("extract_temporal has convergence columns", {
  obj <- make_mock_bvarnet("bernoulli")
  out <- extract_temporal(obj)
  expect_true(all(c("rhat", "ess_bulk", "ess_tail") %in% names(out)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 — extract_random_effects()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_random_effects errors when n_re = 0", {
  obj <- make_mock_bvarnet("bernoulli", n_re = 0L)
  expect_error(extract_random_effects(obj), "n_re = 0")
})

test_that("extract_random_effects what='sd' returns data.frame", {
  obj <- make_mock_bvarnet("bernoulli", n_re = 2L, J = 5L)
  out <- extract_random_effects(obj, what = "sd")
  expect_true(is.data.frame(out))
  expect_true(all(out$type == "Random Effect SD"))
  expect_equal(nrow(out), obj$standata$p * obj$standata$n_re)
})

test_that("extract_random_effects what='mean_u' returns 3D array", {
  obj <- make_mock_bvarnet("bernoulli", n_re = 2L, J = 5L)
  out <- extract_random_effects(obj, what = "mean_u")
  expect_true(is.array(out))
  expect_equal(length(dim(out)), 3L)
  expect_equal(dim(out), c(obj$standata$p, obj$standata$J, obj$standata$n_re))
})

test_that("extract_random_effects what='draws_u' returns 4D array", {
  obj <- make_mock_bvarnet("bernoulli", n_re = 2L, J = 5L)
  out <- extract_random_effects(obj, what = "draws_u")
  expect_true(is.array(out))
  expect_equal(length(dim(out)), 4L)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 — extract_network_matrix()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_network_matrix returns named p x p matrix", {
  obj <- make_mock_bvarnet("bernoulli")
  mat <- extract_network_matrix(obj)
  expect_true(is.matrix(mat))
  expect_equal(dim(mat), c(obj$standata$p, obj$standata$p))
  expect_false(is.null(rownames(mat)))
  expect_false(is.null(colnames(mat)))
  expect_equal(rownames(mat), colnames(mat))
})

test_that("extract_network_matrix respects stat argument", {
  obj  <- make_mock_bvarnet("gaussian")
  m1 <- extract_network_matrix(obj, stat = "mean")
  m2 <- extract_network_matrix(obj, stat = "median")
  # They should both be p x p but may have different values
  expect_equal(dim(m1), dim(m2))
})

test_that("extract_network_matrix errors for invalid lag", {
  obj <- make_mock_bvarnet("bernoulli")  # K=1
  expect_error(extract_network_matrix(obj, lag = 2))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §5 — extract_draws() (now exported)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_draws is accessible as an exported function", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_draws(obj, "beta")
  expect_true(is.matrix(res))
  expect_true(ncol(res) > 0)
})
