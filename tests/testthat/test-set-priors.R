# ─────────────────────────────────────────────────────────────────────────────
# test-set-priors.R — tests for prior(), set_priors(), format(), print(),
#                     and get_default_priors()
# ─────────────────────────────────────────────────────────────────────────────

# ═════════════════════════════════════════════════════════════════════════════
# §1 — prior() construction
# ═════════════════════════════════════════════════════════════════════════════

test_that("prior('normal') returns a bvarnet_prior object", {
  p <- prior("normal", 0, 1)
  expect_s3_class(p, "bvarnet_prior")
})

test_that("prior('normal') has correct family_int", {
  expect_equal(prior("normal", 0, 1)$family_int, 1L)
})

test_that("prior('student_t') has correct family_int and preserves df", {
  p <- prior("student_t", 0, 1, df = 7)
  expect_equal(p$family_int, 2L)
  expect_equal(p$df, 7)
})

test_that("prior('cauchy') has correct family_int", {
  expect_equal(prior("cauchy", 0, 1)$family_int, 3L)
})

test_that("all Phase-1 families produce objects with expected family_int values", {
  supported <- list(
    list("normal",    1L),
    list("student_t", 2L),
    list("cauchy",    3L)
  )
  for (s in supported) {
    p <- prior(s[[1]], 0, 1)
    expect_equal(p$family_int, s[[2]], info = s[[1]])
  }
})

test_that("prior() built directly by user has is_default == FALSE", {
  expect_false(prior("normal", 0, 1)$is_default)
})

test_that("prior() stores loc and scale correctly", {
  p <- prior("normal", loc = 0.5, scale = 2)
  expect_equal(p$loc,   0.5)
  expect_equal(p$scale, 2)
})

test_that("unrecognised family stops with informative error", {
  expect_error(prior("horseshoe"), "Unrecognised prior family")
  expect_error(prior("flat"),      "Unrecognised prior family")
  expect_error(prior("laplace"),   "Unrecognised prior family")
})

test_that("scale <= 0 stops", {
  expect_error(prior("normal", 0, 0),  "scale")
  expect_error(prior("normal", 0, -1), "scale")
})

test_that("df <= 0 for student_t stops", {
  expect_error(prior("student_t", 0, 1, df = 0),  "df")
  expect_error(prior("student_t", 0, 1, df = -1), "df")
})

test_that("very large scale produces a warning", {
  expect_warning(prior("normal", 0, 100), "very large")
})

test_that("very small scale produces a warning", {
  expect_warning(prior("normal", 0, 0.001), "strongly informative")
})

test_that("non-student_t family ignores df (sets it to 0)", {
  expect_equal(prior("normal", 0, 1)$df, 0)
  expect_equal(prior("cauchy", 0, 1)$df, 0)
})


# ═════════════════════════════════════════════════════════════════════════════
# §2 — set_priors() defaults
# ═════════════════════════════════════════════════════════════════════════════

test_that("set_priors() returns a bvarnet_priors object", {
  expect_s3_class(set_priors(), "bvarnet_priors")
})

test_that("set_priors() has the correct names", {
  expect_named(set_priors(), c("beta", "phi", "sd_u", "kappa", "sigma"))
})

test_that("default beta is Normal(0, 1)", {
  sp <- set_priors()
  expect_equal(sp$beta$family_int, 1L)
  expect_equal(sp$beta$loc,   0)
  expect_equal(sp$beta$scale, 1)
})

test_that("default phi scale is 0.5", {
  expect_equal(set_priors()$phi$scale, 0.5)
})

test_that("default sd_u scale is 1", {
  expect_equal(set_priors()$sd_u$scale, 1)
})

test_that("default sigma scale is 2.5", {
  expect_equal(set_priors()$sigma$scale, 2.5)
})

test_that("default kappa scale is 2", {
  expect_equal(set_priors()$kappa$scale, 2)
})

test_that("all defaults have is_default == TRUE", {
  sp <- set_priors()
  for (nm in names(sp)) {
    expect_true(sp[[nm]]$is_default, info = nm)
  }
})

test_that("overriding beta sets is_default to FALSE for beta only", {
  sp <- set_priors(beta = prior("cauchy", 0, 0.5))
  expect_equal(sp$beta$family_int, 3L)
  expect_false(sp$beta$is_default)
  # Other slots unchanged
  expect_true(sp$phi$is_default)
  expect_true(sp$sd_u$is_default)
})

test_that("overriding phi does not affect other slots", {
  sp <- set_priors(phi = prior("student_t", 0, 0.3, df = 5))
  expect_equal(sp$phi$family_int, 2L)
  expect_equal(sp$phi$df, 5)
  expect_true(sp$beta$is_default)
})

test_that("passing non-bvarnet_prior to set_priors() stops", {
  expect_error(set_priors(beta = "normal"), "bvarnet_prior")
  expect_error(set_priors(phi  = list(family = "normal")), "bvarnet_prior")
})


# ═════════════════════════════════════════════════════════════════════════════
# §3 — format.bvarnet_prior() and print.bvarnet_priors()
# ═════════════════════════════════════════════════════════════════════════════

test_that("format(prior('normal', 0, 1)) equals 'Normal(0, 1)'", {
  expect_equal(format(prior("normal", 0, 1)), "Normal(0, 1)")
})

test_that("format with half = TRUE prepends 'Half-'", {
  expect_equal(format(prior("normal", 0, 1), half = TRUE), "Half-Normal(0, 1)")
})

test_that("format(prior('student_t', ...)) contains 'Student-t'", {
  expect_match(format(prior("student_t", 0, 1, df = 7)), "Student-t")
})

test_that("format(prior('cauchy', ...)) contains 'Cauchy'", {
  expect_match(format(prior("cauchy", 0, 1)), "Cauchy")
})

test_that("print(set_priors()) outputs 'beta', 'phi', 'sd_u' without error", {
  out <- capture.output(print(set_priors()))
  expect_true(any(grepl("beta",  out)))
  expect_true(any(grepl("phi",   out)))
  expect_true(any(grepl("sd_u",  out)))
})

test_that("print(set_priors()) shows 'Normal' in output", {
  out <- capture.output(print(set_priors()))
  expect_true(any(grepl("Normal", out)))
})

test_that("print(set_priors()) shows 'Half-' for sd_u and sigma", {
  out <- capture.output(print(set_priors()))
  expect_true(any(grepl("Half-", out)))
})


# ═════════════════════════════════════════════════════════════════════════════
# §4 — get_default_priors()
# ═════════════════════════════════════════════════════════════════════════════

test_that("get_default_priors('bernoulli') returns a bvarnet_priors object", {
  expect_s3_class(get_default_priors("bernoulli"), "bvarnet_priors")
})

test_that("get_default_priors('ordinal') has a kappa slot", {
  sp <- get_default_priors("ordinal")
  expect_true("kappa" %in% names(sp))
  expect_s3_class(sp$kappa, "bvarnet_prior")
})

test_that("get_default_priors('gaussian') has a sigma slot", {
  sp <- get_default_priors("gaussian")
  expect_true("sigma" %in% names(sp))
  expect_s3_class(sp$sigma, "bvarnet_prior")
})

test_that("get_default_priors() with invalid family stops", {
  expect_error(get_default_priors("poisson"), "arg")
})
