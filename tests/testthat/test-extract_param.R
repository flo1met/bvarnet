# ──────────────────────────────────────────────────────────────────────────────
# test-extract_param.R — tests for extract_param() and extract_draws()
#
# These functions require a fitted bvarnet object (Stan fit), which is
# heavy and slow to produce. This file uses mock objects to test the
# structure and logic without running MCMC. Full integration tests with
# real Stan fits should be gated behind skip_on_cran().
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 extract_draws() — class validation
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_draws rejects non-bvarnet objects", {
  expect_error(
    extract_draws(list(a = 1), "beta"),
    "inherits"
  )
})


test_that("extract_draws validates parameter argument", {
  # Even with a mock object, match.arg should reject invalid param
  mock_obj <- structure(list(), class = "bvarnet")

  expect_error(
    extract_draws(mock_obj, "invalid_param"),
    "arg"
  )
})


test_that("extract_draws rejects sigma for non-gaussian family", {
  mock_obj <- structure(
    list(family = "bernoulli", fit = NULL, standata = list()),
    class = "bvarnet"
  )

  expect_error(
    extract_draws(mock_obj, "sigma"),
    "gaussian"
  )
})


test_that("extract_draws rejects kappa for non-ordinal family", {
  mock_obj <- structure(
    list(family = "bernoulli", fit = NULL, standata = list()),
    class = "bvarnet"
  )

  expect_error(
    extract_draws(mock_obj, "kappa"),
    "ordinal"
  )
})


test_that("extract_draws rejects sd_u when n_re = 0", {
  mock_obj <- structure(
    list(
      family   = "gaussian",
      fit      = NULL,
      standata = list(n_re = 0, p = 2, K = 1, n_fe = 1,
                      Y = matrix(0, 1, 2), X = matrix(0, 1, 1),
                      B = matrix(0, 1, 2), Z = matrix(0, 1, 0))
    ),
    class = "bvarnet"
  )

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
# §4 print.bvarnet_params (structural test)
# ═══════════════════════════════════════════════════════════════════════════════

test_that("print.bvarnet_params runs without error on valid structure", {
  mock_params <- structure(
    list(
      beta = data.frame(
        type = c("Intercept", "Fixed Effect"),
        predictor = c("Intercept", "x_1"),
        outcome = c("y_1", "y_1"),
        mean = c(0.5, -0.3),
        median = c(0.48, -0.31),
        q5 = c(0.1, -0.7),
        q95 = c(0.9, 0.1),
        stringsAsFactors = FALSE
      ),
      phi = data.frame(
        type = "Temporal",
        predictor = "lag1_y_1",
        outcome = "y_1",
        mean = 0.3,
        median = 0.29,
        q5 = 0.1,
        q95 = 0.5,
        stringsAsFactors = FALSE
      ),
      re_sd = data.frame(
        type = "Random Effect SD",
        outcome = "y_1",
        random_effect = "x_1",
        mean = 0.5,
        median = 0.49,
        q5 = 0.2,
        q95 = 0.8,
        stringsAsFactors = FALSE
      ),
      sigma = data.frame(
        type = "Residual SD",
        predictor = "y_1",
        outcome = "sigma",
        mean = 1.0,
        median = 0.99,
        q5 = 0.8,
        q95 = 1.2,
        stringsAsFactors = FALSE
      ),
      standata = list(),
      fit = NULL
    ),
    class = "bvarnet_params"
  )

  expect_output(print(mock_params), "BVAR Network")
})


test_that("print.bvarnet_params handles NULL re_sd and sigma", {
  mock_params <- structure(
    list(
      beta = data.frame(
        type = "Intercept", predictor = "Int", outcome = "y1",
        mean = 0, median = 0, q5 = -1, q95 = 1,
        stringsAsFactors = FALSE
      ),
      phi = data.frame(
        type = "Temporal", predictor = "l", outcome = "y1",
        mean = 0, median = 0, q5 = -1, q95 = 1,
        stringsAsFactors = FALSE
      ),
      re_sd = NULL,
      sigma = NULL,
      standata = list(),
      fit = NULL
    ),
    class = "bvarnet_params"
  )

  expect_output(print(mock_params), "BVAR Network")
})
