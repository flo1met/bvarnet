# ──────────────────────────────────────────────────────────────────────────────
# test-extract_param.R — tests for extract_param() and extract_draws()
#
# These functions require a fitted bvarnet object (Stan fit), which is
# heavy and slow to produce. This file uses mock objects to test the
# structure and logic without running MCMC. Full integration tests with
# real Stan fits should be gated behind skip_on_cran().
# make_mock_bvarnet() is defined in helper-fixtures.R and auto-sourced.
# ──────────────────────────────────────────────────────────────────────────────


# Also note: extract_draws is internal. When sourced (test_dir), it lives in the
# global env; when the package is loaded (devtools::test), use bvarnet:::.
# The helper always sources local R files so plain extract_draws() works here.

# ═══════════════════════════════════════════════════════════════════════════════
# §1 extract_draws() — class and family validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_draws rejects non-bvarnet objects", {
  expect_error(
    extract_draws(list(a = 1), "beta"),
    "inherits"
  )
})


test_that("extract_draws validates parameter argument", {
  mock_obj <- structure(list(), class = "bvarnet")

  expect_error(
    extract_draws(mock_obj, "invalid_param"),
    "arg"
  )
})


test_that("extract_draws rejects sigma for non-gaussian family", {
  mock_obj <- make_mock_bvarnet("bernoulli")

  expect_error(
    extract_draws(mock_obj, "sigma"),
    "gaussian"
  )
})


test_that("extract_draws rejects kappa for non-ordinal family", {
  mock_obj <- make_mock_bvarnet("bernoulli")

  expect_error(
    extract_draws(mock_obj, "kappa"),
    "ordinal"
  )
})


test_that("extract_draws rejects sd_u when n_re = 0", {
  mock_obj <- make_mock_bvarnet("gaussian")   # n_re = 0 by default

  expect_error(
    extract_draws(mock_obj, "sd_u"),
    "no random effects"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 extract_param() — class validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_param rejects non-bvarnet objects", {
  expect_error(
    extract_param(list(a = 1)),
    "inherits"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 compare_to_truth() — class validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("compare_to_truth rejects non-bvarnet fit object", {
  expect_error(
    compare_to_truth(list(a = 1), list()),
    "inherits"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 extract_param() — return structure
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_param returns a plain data.frame for bernoulli", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  expect_true(is.data.frame(res))
  expect_identical(class(res), "data.frame")   # not classed
})


test_that("extract_param result has all required columns", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  expected_cols <- c("type", "predictor", "outcome",
                     "mean", "median", "q5", "q95",
                     "rhat", "ess_bulk", "ess_tail")
  expect_true(all(expected_cols %in% names(res)))
})


test_that("extract_param includes beta and phi rows for bernoulli", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  expect_true(any(res$type %in% c("Intercept", "Fixed Effect")))
  expect_true(any(res$type == "Autoregressive"))
  expect_true(any(res$type == "Cross-lagged"))
})


test_that("extract_param does not include sigma/kappa for bernoulli", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  expect_false(any(res$type == "Residual SD"))
  expect_false(any(res$type == "Threshold"))
})


test_that("extract_param includes sigma rows for gaussian", {
  obj <- make_mock_bvarnet("gaussian")
  res <- extract_param(obj)

  expect_true(any(res$type == "Residual SD"))
  expect_equal(sum(res$type == "Residual SD"), obj$standata$p)
})


test_that("extract_param includes kappa rows for ordinal", {
  obj <- make_mock_bvarnet("ordinal")
  res <- extract_param(obj)

  expect_true(any(res$type == "Threshold"))
  # p=2, C-1=2 cutpoints: expect 4 threshold rows
  expect_equal(sum(res$type == "Threshold"),
               obj$standata$p * (obj$standata$C - 1L))
})


test_that("extract_param filtering by type works", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  ar_rows <- subset(res, type == "Autoregressive")
  cl_rows <- subset(res, type == "Cross-lagged")
  expect_true(nrow(ar_rows) > 0)
  expect_true(nrow(cl_rows) > 0)
  expect_true(all(ar_rows$type == "Autoregressive"))
  expect_true(all(cl_rows$type == "Cross-lagged"))
})


test_that("extract_param rhat and ess columns are numeric", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  expect_true(is.numeric(res$rhat))
  expect_true(is.numeric(res$ess_bulk))
  expect_true(is.numeric(res$ess_tail))
})


test_that("extract_param predictor labels use variable names from standata", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj)

  intercept_rows <- subset(res, type == "Intercept")
  expect_true(all(intercept_rows$predictor == "Intercept"))

  phi_rows <- subset(res, type %in% c("Autoregressive", "Cross-lagged"))
  expect_true(all(grepl("lag1_y_[0-9]+", phi_rows$predictor)))
})
