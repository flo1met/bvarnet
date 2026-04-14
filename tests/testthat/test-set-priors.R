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

test_that("print(set_priors()) shows 'all defaults' message", {
  out <- capture.output(print(set_priors()))
  expect_true(any(grepl("all defaults", out)))
})

test_that("print(set_priors()) does not show individual prior lines when all defaults", {
  out <- capture.output(print(set_priors()))
  # Should NOT have individual "beta ~" lines
  expect_false(any(grepl("^\\s*beta\\s+~", out)))
})

test_that("print with one user prior shows only that prior", {
  sp <- set_priors(beta = prior("cauchy", 0, 2))
  out <- capture.output(print(sp))
  # beta should appear on its own line
  expect_true(any(grepl("^\\s*beta\\s+~", out)))
  # defaults should NOT appear
  expect_false(any(grepl("phi", out)))
  expect_false(any(grepl("sd_u", out)))
  expect_false(any(grepl("Defaults:", out)))
})

test_that("print with multiple user priors shows only those", {
  sp <- set_priors(
    beta = prior("cauchy", 0, 2),
    phi  = prior("student_t", 0, 0.3, df = 5)
  )
  out <- capture.output(print(sp))
  expect_true(any(grepl("beta", out)))
  expect_true(any(grepl("phi", out)))
  expect_false(any(grepl("sd_u", out)))
  expect_false(any(grepl("sigma", out)))
})

test_that("print with all user-set priors shows all of them", {
  sp <- set_priors(
    beta  = prior("normal", 0, 2),
    phi   = prior("normal", 0, 1),
    sd_u  = prior("cauchy", 0, 1),
    kappa = prior("normal", 0, 3),
    sigma = prior("cauchy", 0, 2)
  )
  out <- capture.output(print(sp))
  expect_false(any(grepl("all defaults", out)))
  expect_true(any(grepl("beta",  out)))
  expect_true(any(grepl("phi",   out)))
  expect_true(any(grepl("sd_u",  out)))
  expect_true(any(grepl("sigma", out)))
  expect_true(any(grepl("kappa", out)))
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

test_that("get_default_priors('bernoulli') omits sigma and kappa", {
  sp <- get_default_priors("bernoulli")
  expect_null(sp$sigma)
  expect_null(sp$kappa)
  expect_true("beta" %in% names(sp))
  expect_true("phi"  %in% names(sp))
  expect_true("sd_u" %in% names(sp))
})

test_that("get_default_priors('gaussian') omits kappa but keeps sigma", {
  sp <- get_default_priors("gaussian")
  expect_null(sp$kappa)
  expect_true("sigma" %in% names(sp))
})

test_that("get_default_priors('ordinal') omits sigma but keeps kappa", {
  sp <- get_default_priors("ordinal")
  expect_null(sp$sigma)
  expect_true("kappa" %in% names(sp))
})

test_that("get_default_priors(has_re = FALSE) omits sd_u", {
  sp <- get_default_priors("bernoulli", has_re = FALSE)
  expect_null(sp$sd_u)
  expect_true("beta" %in% names(sp))
  expect_true("phi"  %in% names(sp))
})

test_that("get_default_priors() with no args returns all five", {
  sp <- get_default_priors()
  expect_named(sp, c("beta", "phi", "sd_u", "kappa", "sigma"))
})


# ═════════════════════════════════════════════════════════════════════════════
# §5 — .ensure_prior_slots()
# ═════════════════════════════════════════════════════════════════════════════

test_that(".ensure_prior_slots fills missing sigma for gaussian family", {
  p <- get_default_priors("bernoulli")  # no sigma
  family_vec <- c(y_1 = "gaussian")
  expect_warning(
    filled <- bvarnet:::.ensure_prior_slots(p, family_vec),
    "sigma"
  )
  expect_s3_class(filled$sigma, "bvarnet_prior")
  expect_true(filled$sigma$is_default)
})

test_that(".ensure_prior_slots fills missing kappa for ordinal family", {
  p <- get_default_priors("bernoulli")  # no kappa
  family_vec <- c(y_1 = "ordinal")
  expect_warning(
    filled <- bvarnet:::.ensure_prior_slots(p, family_vec),
    "kappa"
  )
  expect_s3_class(filled$kappa, "bvarnet_prior")
})

test_that(".ensure_prior_slots fills missing sd_u from filtered object", {
  p <- get_default_priors("bernoulli", has_re = FALSE)  # no sd_u
  family_vec <- c(y_1 = "bernoulli")
  expect_warning(
    filled <- bvarnet:::.ensure_prior_slots(p, family_vec),
    "sd_u"
  )
  expect_s3_class(filled$sd_u, "bvarnet_prior")
})

test_that(".ensure_prior_slots does nothing for complete priors", {
  p <- set_priors()
  family_vec <- c(y_1 = "gaussian", y_2 = "ordinal")
  # No warning expected
  expect_silent(bvarnet:::.ensure_prior_slots(p, family_vec))
})


# ═════════════════════════════════════════════════════════════════════════════
# §6 — .prior_warnings()
# ═════════════════════════════════════════════════════════════════════════════

test_that(".prior_warnings returns correct needed set for bernoulli without REs", {
  p <- set_priors()
  needed <- bvarnet:::.prior_warnings(p, c(y_1 = "bernoulli"), n_re = 0L)
  expect_equal(needed, c("beta", "phi"))
})

test_that(".prior_warnings returns correct needed set for gaussian with REs", {
  p <- set_priors()
  needed <- bvarnet:::.prior_warnings(p, c(y_1 = "gaussian"), n_re = 2L)
  expect_equal(needed, c("beta", "phi", "sd_u", "sigma"))
})

test_that(".prior_warnings returns correct needed set for ordinal", {
  p <- set_priors()
  needed <- bvarnet:::.prior_warnings(p, c(y_1 = "ordinal"), n_re = 0L)
  expect_equal(needed, c("beta", "phi", "kappa"))
})

test_that(".prior_warnings warns about unused user-set priors", {
  p <- set_priors(sigma = prior("cauchy", 0, 2))
  expect_warning(
    bvarnet:::.prior_warnings(p, c(y_1 = "bernoulli"), n_re = 0L),
    "not used.*sigma"
  )
})

test_that(".prior_warnings warns about unused sd_u when no REs", {
  p <- set_priors(sd_u = prior("cauchy", 0, 1))
  expect_warning(
    bvarnet:::.prior_warnings(p, c(y_1 = "bernoulli"), n_re = 0L),
    "not used.*sd_u"
  )
})

test_that(".prior_warnings messages about auto-filled defaults", {
  p <- set_priors(beta = prior("cauchy", 0, 2))
  expect_message(
    bvarnet:::.prior_warnings(p, c(y_1 = "gaussian"), n_re = 0L),
    "Using default priors for: phi, sigma"
  )
})

test_that(".prior_warnings does NOT message when all priors are defaults", {
  p <- set_priors()
  expect_silent(
    bvarnet:::.prior_warnings(p, c(y_1 = "gaussian"), n_re = 0L)
  )
})

test_that(".prior_warnings does NOT message when all needed priors are user-set", {
  p <- set_priors(
    beta = prior("cauchy", 0, 2),
    phi  = prior("normal", 0, 0.3)
  )
  expect_silent(
    bvarnet:::.prior_warnings(p, c(y_1 = "bernoulli"), n_re = 0L)
  )
})
