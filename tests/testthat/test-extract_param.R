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
# §1b extract_draws() — multiple parameter blocks
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_draws returns several blocks side by side in one matrix", {
  obj <- make_mock_bvarnet("gaussian")
  res <- extract_draws(obj, c("beta", "phi", "sigma"))

  expect_true(is.matrix(res))
  expect_equal(nrow(res), 40L)   # n_iter (20) * n_chains (2)
  expect_equal(
    ncol(res),
    ncol(extract_draws(obj, "beta")) +
      ncol(extract_draws(obj, "phi")) +
      ncol(extract_draws(obj, "sigma"))
  )
  # Columns keep their Stan names and appear in the order requested
  expect_equal(colnames(res)[1L], "beta[1,1]")
  expect_true(all(grepl("^(beta|phi|sigma)\\[", colnames(res))))
  expect_equal(res[, "phi[2,1]"], extract_draws(obj, "phi")[, "phi[2,1]"])
})


test_that("extract_draws deduplicates repeated parameter names", {
  obj <- make_mock_bvarnet("bernoulli")

  expect_equal(
    extract_draws(obj, c("phi", "phi")),
    extract_draws(obj, "phi")
  )
})


test_that("extract_draws supports partial matching of parameter names", {
  obj <- make_mock_bvarnet("gaussian")

  expect_equal(extract_draws(obj, "sig"), extract_draws(obj, "sigma"))
})


test_that("extract_draws rejects unknown parameter names in a vector", {
  obj <- make_mock_bvarnet("bernoulli")

  # A bad name must not be silently dropped just because a good one matched
  expect_error(extract_draws(obj, c("phi", "bogus")), "bogus")
})


test_that("extract_draws rejects 'all' combined with other names", {
  obj <- make_mock_bvarnet("bernoulli")

  expect_error(extract_draws(obj, c("all", "phi")), "cannot be combined")
})


test_that("extract_draws 'all' returns every block the model has", {
  obj <- make_mock_bvarnet("gaussian", n_re = 2L)
  res <- extract_draws(obj, "all")

  expect_true(is.matrix(res))
  expect_true(all(c("beta[1,1]", "phi[1,1]", "sd_u[1,1]", "sigma[1]") %in%
                    colnames(res)))
  # kappa belongs to ordinal models only, and u is not one of the blocks
  expect_false(any(grepl("^(kappa|u)\\[", colnames(res))))
})


test_that("extract_draws defaults to 'all'", {
  obj <- make_mock_bvarnet("ordinal")

  expect_equal(extract_draws(obj), extract_draws(obj, "all"))
})


test_that("extract_draws does not offer u — that is extract_random_effects' job", {
  obj <- make_mock_bvarnet("gaussian", n_re = 2L, J = 5L)

  expect_error(extract_draws(obj, "u"), "Unknown `parameter` value")
  expect_false("u" %in% .draw_param_choices)
  expect_false(any(grepl("^u\\[", colnames(extract_draws(obj, "all")))))
})


test_that(".extract_draws_block still reaches u internally", {
  # u is off the public menu but must stay available to package internals,
  # which is why .extract_draws_block() does not consult .draw_param_choices.
  obj <- make_mock_bvarnet("gaussian", n_re = 2L, J = 5L)
  res <- .extract_draws_block(obj, "u")

  expect_true(is.matrix(res))
  sd <- obj$standata
  expect_equal(ncol(res), sd$p * sd$J * sd$n_re)
  expect_true(all(grepl("^u\\[", colnames(res))))
  # Same draws the shaped extractor reports, just flattened into columns
  expect_equal(unname(res[, "u[2,3,1]"]),
               unname(.extract_u_draws(obj)[, 2L, 3L, 1L]))
})


test_that("extract_draws returns lp__ as a one-column matrix", {
  obj <- add_mock_lp(make_mock_bvarnet("bernoulli"))
  res <- extract_draws(obj, "lp__")

  expect_true(is.matrix(res))
  expect_equal(dim(res), c(40L, 1L))
  expect_equal(colnames(res), "lp__")
  # lp__ is the one name matched whole rather than by an "[index]" suffix
  expect_equal(as.numeric(res), as.numeric(obj$draws[, , "lp__"]))
})


test_that("extract_draws errors on lp__ when the draws do not carry it", {
  obj <- make_mock_bvarnet("bernoulli")   # mocks have no lp__ unless added

  expect_error(extract_draws(obj, "lp__"), "lp__")
})


test_that("extract_draws 'all' includes lp__ and orders blocks", {
  obj <- add_mock_lp(make_mock_bvarnet("ordinal", n_re = 2L))
  res <- extract_draws(obj, "all")

  block <- sub("\\[.*$", "", colnames(res))
  expect_equal(unique(block), c("beta", "phi", "sd_u", "kappa", "lp__"))
  expect_equal(ncol(res), sum(vapply(
    c("beta", "phi", "sd_u", "kappa", "lp__"),
    function(p) ncol(extract_draws(obj, p)), integer(1L)
  )))
})


test_that("extract_draws 'all' skips blocks this model does not have", {
  # Pure ordinal with no covariates: no sigma (not gaussian), no sd_u
  # (n_re = 0), and beta is declared matrix[0, p] so it has no draws.
  obj <- make_mock_ordinal_no_fe()
  res <- extract_draws(obj, "all")

  expect_true(all(grepl("^(phi|kappa)\\[", colnames(res))))
  expect_false(any(grepl("^beta\\[", colnames(res))))
  # Naming an absent block explicitly is still an error
  expect_error(extract_draws(obj, c("phi", "sigma")), "gaussian")
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
                     "mean", "median", "ci_lower", "ci_upper",
                     "rhat", "ess_bulk", "ess_tail")
  expect_true(all(expected_cols %in% names(res)))
})


test_that("extract_param ci_level widens/narrows the interval", {
  obj <- make_mock_bvarnet("bernoulli")

  res95 <- extract_param(obj)
  res50 <- extract_param(obj, ci_level = 0.50)
  res99 <- extract_param(obj, ci_level = 0.99)

  expect_true(all(res50$ci_lower >= res95$ci_lower))
  expect_true(all(res50$ci_upper <= res95$ci_upper))
  expect_true(all(res99$ci_lower <= res95$ci_lower))
  expect_true(all(res99$ci_upper >= res95$ci_upper))
})


test_that("extract_param rejects an invalid ci_level", {
  obj <- make_mock_bvarnet("bernoulli")
  expect_error(extract_param(obj, ci_level = 95), "strictly between 0 and 1")
  expect_error(extract_param(obj, ci_level = 0), "strictly between 0 and 1")
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


test_that("extract_param retains FE rows for pure ordinal (no sentinel filtering)", {
  obj <- make_mock_bvarnet("ordinal")
  res <- extract_param(obj)

  # Pure ordinal: beta[1,j] is x_1, not an intercept sentinel
  fe_rows <- res[res$type == "Fixed Effect", ]
  expect_equal(nrow(fe_rows), obj$standata$p)
  expect_true(all(fe_rows$predictor == "x_1"))
  expect_false(any(is.na(fe_rows$mean)))
  # No Intercept rows expected

  expect_equal(sum(res$type == "Intercept"), 0L)
})


test_that("extract_draws returns a zero-column matrix when n_fe == 0", {
  obj <- make_mock_ordinal_no_fe()
  d <- extract_draws(obj, "beta")

  expect_true(is.matrix(d))
  expect_equal(ncol(d), 0L)
  expect_equal(nrow(d), 40L)   # n_iter (20) * n_chains (2)
})


test_that("extract_param works for pure ordinal with no fixed effects", {
  obj <- make_mock_ordinal_no_fe()
  res <- extract_param(obj)

  expect_s3_class(res, "data.frame")
  expect_gt(nrow(res), 0L)
  expect_equal(sum(res$type %in% c("Intercept", "Fixed Effect")), 0L)
  expect_true(is.character(res$type))
  expect_true(all(c("Autoregressive", "Cross-lagged", "Threshold") %in% res$type))
  expect_equal(sum(res$type == "Threshold"),
               obj$standata$p * (obj$standata$C - 1L))
})


test_that("summary() works for pure ordinal with no fixed effects", {
  obj <- make_mock_ordinal_no_fe()

  expect_no_error(s <- summary(obj))
  expect_s3_class(s, "summary.bvarnet")
  expect_output(print(s), "Threshold")
})


test_that(".summarize_draws handles a zero-column draws matrix", {
  s <- .summarize_draws(matrix(numeric(0), nrow = 10L, ncol = 0L), c(0.025, 0.975))

  expect_equal(s$mean,     numeric(0))
  expect_equal(s$median,   numeric(0))
  expect_equal(s$ci_lower, numeric(0))
  expect_equal(s$ci_upper, numeric(0))
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


# ═══════════════════════════════════════════════════════════════════════════════
# §N — extract_param() type= argument
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_param type=NULL returns all rows (backward compat)", {
  obj  <- make_mock_bvarnet("gaussian")
  full <- extract_param(obj)
  with_null <- extract_param(obj, type = NULL)
  expect_equal(full, with_null)
})

test_that("extract_param type='Threshold' returns only threshold rows for ordinal", {
  obj <- make_mock_bvarnet("ordinal")
  res <- extract_param(obj, type = "Threshold")
  expect_true(is.data.frame(res))
  expect_true(nrow(res) > 0L)
  expect_true(all(res$type == "Threshold"))
})

test_that("extract_param type='Residual SD' returns only residual SD rows for gaussian", {
  obj <- make_mock_bvarnet("gaussian")
  res <- extract_param(obj, type = "Residual SD")
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), obj$standata$p)
  expect_true(all(res$type == "Residual SD"))
})

test_that("extract_param type= returns empty data.frame when type absent in model", {
  obj <- make_mock_bvarnet("bernoulli")
  res <- extract_param(obj, type = "Threshold")
  expect_true(is.data.frame(res))
  expect_equal(nrow(res), 0L)
})

test_that("extract_param type= accepts multiple types", {
  obj <- make_mock_bvarnet("gaussian")
  res <- extract_param(obj, type = c("Autoregressive", "Cross-lagged"))
  expect_true(all(res$type %in% c("Autoregressive", "Cross-lagged")))
  expect_true(nrow(res) > 0L)
})

test_that("extract_param type= errors on unknown type value", {
  obj <- make_mock_bvarnet("bernoulli")
  expect_error(
    extract_param(obj, type = "NotAType"),
    "Unknown type value"
  )
})
