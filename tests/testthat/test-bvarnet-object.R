# ════──────────────────────────────────────────────────────────────────────────
# test-bvarnet-object.R — tests for the bvarnet object structure
#
# Covers: object slots, draws array properties, extract_draws() subsetting,
# and print.bvarnet().
# make_mock_bvarnet() is defined in helper-fixtures.R and auto-sourced.
# extract_draws() is internal; the helper always sources local R files so
# it is available in the global env as plain extract_draws().
# ════──────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 — Object structure
# ═══════════════════════════════════════════════════════════════════════════════

test_that("bvarnet object has class 'bvarnet'", {
  obj <- make_mock_bvarnet()
  expect_s3_class(obj, "bvarnet")
})


test_that("bvarnet object contains expected top-level names", {
  obj <- make_mock_bvarnet()
  expected <- c("draws", "convergence", "diagnostics", "timing",
                "metadata", "return_codes", "family", "standata", "priors")
  expect_true(all(expected %in% names(obj)))
})


test_that("bvarnet object does not contain a 'fit' field", {
  obj <- make_mock_bvarnet()
  expect_null(obj$fit)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 — draws slot
# ═══════════════════════════════════════════════════════════════════════════════

test_that("draws slot is a plain base-R array", {
  obj <- make_mock_bvarnet()
  expect_true(is.array(obj$draws))
  expect_false(inherits(obj$draws, "draws_array"))
  expect_identical(class(obj$draws), "array")
})


test_that("draws slot has 3 dimensions", {
  obj <- make_mock_bvarnet()
  expect_equal(length(dim(obj$draws)), 3L)
})


test_that("draws dimnames[[3]] is a non-null character vector", {
  obj <- make_mock_bvarnet()
  dn3 <- dimnames(obj$draws)[[3]]
  expect_false(is.null(dn3))
  expect_true(is.character(dn3))
  expect_true(length(dn3) > 0)
})


test_that("draws parameter names are in Stan bracket format", {
  obj <- make_mock_bvarnet()
  dn3 <- dimnames(obj$draws)[[3]]
  expect_true(all(grepl("^[a-z_]+\\[[0-9,]+\\]$", dn3)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 — convergence slot
# ═══════════════════════════════════════════════════════════════════════════════

test_that("convergence slot is a plain data.frame (not tibble)", {
  obj <- make_mock_bvarnet()
  expect_identical(class(obj$convergence), "data.frame")
})


test_that("convergence slot has required convergence columns", {
  obj  <- make_mock_bvarnet()
  cols <- names(obj$convergence)
  expect_true("variable"  %in% cols)
  expect_true("rhat"      %in% cols)
  expect_true("ess_bulk"  %in% cols)
  expect_true("ess_tail"  %in% cols)
})


test_that("convergence slot has one row per parameter", {
  obj <- make_mock_bvarnet("bernoulli")
  expect_equal(nrow(obj$convergence), dim(obj$draws)[3])
})


test_that("convergence rhat values are numeric", {
  obj <- make_mock_bvarnet()
  expect_true(is.numeric(obj$convergence$rhat))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 — diagnostics slot
# ═══════════════════════════════════════════════════════════════════════════════

test_that("diagnostics slot is a data.frame", {
  obj <- make_mock_bvarnet()
  expect_true(is.data.frame(obj$diagnostics))
})


test_that("diagnostics slot has required columns", {
  obj  <- make_mock_bvarnet()
  cols <- names(obj$diagnostics)
  expect_true("num_divergent"     %in% cols)
  expect_true("num_max_treedepth" %in% cols)
  expect_true("ebfmi"             %in% cols)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §5 — extract_draws() base-R subsetting
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_draws returns a matrix with Stan column names", {
  obj    <- make_mock_bvarnet("bernoulli")
  result <- extract_draws(obj, "beta")

  expect_true(is.matrix(result))
  expect_true(all(grepl("^beta\\[", colnames(result))))
})


test_that("extract_draws returns correct row count (iter * chains)", {
  n_iter <- 20L; n_chains <- 2L
  obj    <- make_mock_bvarnet("bernoulli", n_iter = n_iter, n_chains = n_chains)
  result <- extract_draws(obj, "beta")

  expect_equal(nrow(result), n_iter * n_chains)
})


test_that("extract_draws returns correct column count for beta", {
  obj    <- make_mock_bvarnet("bernoulli")  # p=2, n_fe=2 → 2*2=4 beta params
  result <- extract_draws(obj, "beta")

  expect_equal(ncol(result), 4L)
})


test_that("extract_draws for 'beta' does not include phi, sigma, or kappa columns", {
  obj    <- make_mock_bvarnet("gaussian")
  result <- extract_draws(obj, "beta")

  expect_false(any(grepl("^phi\\[",   colnames(result))))
  expect_false(any(grepl("^sigma\\[", colnames(result))))
})


test_that("extract_draws for sigma returns correct column count", {
  obj    <- make_mock_bvarnet("gaussian")   # p=2 → 2 sigma params
  result <- extract_draws(obj, "sigma")

  expect_equal(ncol(result), 2L)
  expect_true(all(grepl("^sigma\\[", colnames(result))))
})


test_that("extract_draws for kappa returns correct column count", {
  obj    <- make_mock_bvarnet("ordinal")    # p=2, C-1=2 → 4 kappa params
  result <- extract_draws(obj, "kappa")

  expect_equal(ncol(result), 4L)
  expect_true(all(grepl("^kappa\\[", colnames(result))))
})


test_that("extract_draws column values match expected draws slice", {
  obj  <- make_mock_bvarnet("bernoulli")
  res  <- extract_draws(obj, "beta")
  drws <- obj$draws

  # Manually replicate what extract_draws should do
  idx   <- grep("^beta\\[", dimnames(drws)[[3]])
  chunk <- drws[, , idx, drop = FALSE]
  dim(chunk) <- c(prod(dim(chunk)[1:2]), dim(chunk)[3])
  colnames(chunk) <- dimnames(drws)[[3]][idx]

  expect_equal(res, chunk)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §6 — print.bvarnet()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("print.bvarnet outputs model family", {
  obj <- make_mock_bvarnet("gaussian")
  expect_output(print(obj), "gaussian")
})


test_that("print.bvarnet outputs dimension information", {
  obj <- make_mock_bvarnet("bernoulli")
  expect_output(print(obj), "2")    # p = 2 appears
})


test_that("print.bvarnet outputs convergence information", {
  obj <- make_mock_bvarnet()
  expect_output(print(obj), "Rhat")
})


test_that("print.bvarnet outputs divergence information", {
  obj <- make_mock_bvarnet()
  expect_output(print(obj), "Divergences")
})


test_that("print.bvarnet outputs timing", {
  obj <- make_mock_bvarnet()
  expect_output(print(obj), "5")   # timing$total = 5.0
})


test_that("print.bvarnet warns when divergences > 0", {
  obj <- make_mock_bvarnet()
  obj$diagnostics$num_divergent <- c(3L, 0L)
  expect_output(print(obj), "WARNING")
})


test_that("print.bvarnet warns when Rhat > 1.01", {
  obj <- make_mock_bvarnet()
  obj$convergence$rhat <- rep(1.05, nrow(obj$convergence))
  expect_output(print(obj), "WARNING")
})


test_that("print.bvarnet returns object invisibly", {
  obj <- make_mock_bvarnet()
  ret <- withVisible(print(obj))
  expect_false(ret$visible)
  expect_identical(ret$value, obj)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §7 — priors slot
# ═══════════════════════════════════════════════════════════════════════════════

test_that("priors slot is present on mock bvarnet object", {
  obj <- make_mock_bvarnet()
  expect_true("priors" %in% names(obj))
})

test_that("priors slot inherits bvarnet_priors", {
  obj <- make_mock_bvarnet()
  expect_s3_class(obj$priors, "bvarnet_priors")
})

test_that("priors slot has correct default beta", {
  obj <- make_mock_bvarnet()
  expect_equal(obj$priors$beta$family_int, 1L)
  expect_equal(obj$priors$beta$loc,   0)
  expect_equal(obj$priors$beta$scale, 1)
})

test_that("priors slot has all expected named entries", {
  obj <- make_mock_bvarnet()
  expect_named(obj$priors, c("intercept", "beta", "phi", "sd_u", "kappa", "sigma"))
})
