# ──────────────────────────────────────────────────────────────────────────────
# test-helpers.R — tests for helper functions (no Stan fit required)
#
# Tests get_param_names(), build_summary_table(), and the internal
# data manipulation utilities. extract_draws() requires a real Stan fit
# so is tested separately in test-extract_param.R (skip_on_cran).
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 get_param_names()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("get_param_names returns correct structure with named columns", {
  # Create a mock standata-like list
  sd <- list(
    p = 3, K = 1, n_fe = 2, n_re = 1,
    Y = matrix(0, 10, 3, dimnames = list(NULL, c("y1", "y2", "y3"))),
    X = matrix(0, 10, 2, dimnames = list(NULL, c("Intercept", "x1"))),
    B = matrix(0, 10, 3, dimnames = list(NULL, c("lag1_y1", "lag1_y2", "lag1_y3"))),
    Z = matrix(0, 10, 1, dimnames = list(NULL, c("x1")))
  )

  nm <- bvarnet:::get_param_names(sd)

  expect_type(nm, "list")
  expect_named(nm, c("y", "fe", "b", "re"))
  expect_equal(nm$y, c("y1", "y2", "y3"))
  expect_equal(nm$fe, c("Intercept", "x1"))
  expect_equal(nm$b, c("lag1_y1", "lag1_y2", "lag1_y3"))
  expect_equal(nm$re, "x1")
})


test_that("get_param_names uses fallback names when colnames are NULL", {
  sd <- list(
    p = 2, K = 1, n_fe = 1, n_re = 0,
    Y = matrix(0, 10, 2),  # no colnames
    X = matrix(0, 10, 1),
    B = matrix(0, 10, 2),
    Z = matrix(0, 10, 0)
  )

  nm <- bvarnet:::get_param_names(sd)

  expect_equal(nm$y, c("y1", "y2"))
  expect_equal(nm$fe, "fe1")
  expect_equal(length(nm$re), 0)
})


test_that("get_param_names with n_re = 0 returns empty re vector", {
  sd <- list(p = 2, K = 1, n_fe = 1, n_re = 0,
             Y = matrix(0, 5, 2, dimnames = list(NULL, c("a", "b"))),
             X = matrix(0, 5, 1, dimnames = list(NULL, "Int")),
             B = matrix(0, 5, 2, dimnames = list(NULL, c("l1_a", "l1_b"))),
             Z = matrix(0, 5, 0))

  nm <- bvarnet:::get_param_names(sd)
  expect_equal(nm$re, character(0))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 build_summary_table()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("build_summary_table returns correct data frame structure", {
  set.seed(500)
  draws <- matrix(rnorm(200 * 6), nrow = 200, ncol = 6)
  row_names <- c("pred1", "pred2")
  col_names <- c("out1", "out2", "out3")

  tab <- bvarnet:::build_summary_table(draws, row_names, col_names, "TestType")

  expect_s3_class(tab, "data.frame")
  expect_equal(nrow(tab), 2 * 3)  # nr * nc
  expect_true(all(c("type", "predictor", "outcome",
                     "mean", "median", "ci_lower", "ci_upper") %in% names(tab)))
})


test_that("build_summary_table type column is filled correctly", {
  draws <- matrix(rnorm(100 * 4), nrow = 100, ncol = 4)
  tab <- bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "MyType")

  expect_true(all(tab$type == "MyType"))
})


test_that("build_summary_table predictor × outcome layout is correct", {
  draws <- matrix(rnorm(100 * 6), nrow = 100, ncol = 6)
  row_names <- c("r1", "r2", "r3")
  col_names <- c("c1", "c2")

  tab <- bvarnet:::build_summary_table(draws, row_names, col_names, "T")

  # Column order: rows cycle within columns
  # (r1,c1), (r2,c1), (r3,c1), (r1,c2), (r2,c2), (r3,c2)
  expect_equal(tab$predictor, rep(row_names, times = 2))
  expect_equal(tab$outcome, rep(col_names, each = 3))
})


test_that("build_summary_table ci_lower <= median <= ci_upper", {
  set.seed(501)
  draws <- matrix(rnorm(500 * 4), nrow = 500, ncol = 4)
  tab <- bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "T")

  expect_true(all(tab$ci_lower <= tab$median))
  expect_true(all(tab$median <= tab$ci_upper))
})


test_that("build_summary_table defaults to a 95% equal-tailed interval", {
  set.seed(502)
  draws <- matrix(rnorm(2000 * 4), nrow = 2000, ncol = 4)
  tab <- bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "T")

  expect_equal(tab$ci_lower,
               unname(apply(draws, 2, stats::quantile, probs = 0.025)))
  expect_equal(tab$ci_upper,
               unname(apply(draws, 2, stats::quantile, probs = 0.975)))
})


test_that("build_summary_table honours ci_level", {
  set.seed(503)
  draws <- matrix(rnorm(2000 * 4), nrow = 2000, ncol = 4)
  tab90 <- bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "T",
                                         ci_level = 0.90)

  expect_equal(tab90$ci_lower,
               unname(apply(draws, 2, stats::quantile, probs = 0.05)))
  expect_equal(tab90$ci_upper,
               unname(apply(draws, 2, stats::quantile, probs = 0.95)))

  # Narrower level nests inside the default 95% interval
  tab95 <- bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "T")
  expect_true(all(tab90$ci_lower >= tab95$ci_lower))
  expect_true(all(tab90$ci_upper <= tab95$ci_upper))
})


test_that(".ci_probs rejects invalid credible-interval levels", {
  expect_equal(bvarnet:::.ci_probs(0.95), c(0.025, 0.975))
  expect_equal(bvarnet:::.ci_probs(0.90), c(0.05, 0.95))

  for (bad in list(0, 1, -0.1, 1.5, NA_real_, NaN, Inf, "0.95", c(0.9, 0.95)))
    expect_error(bvarnet:::.ci_probs(bad), "strictly between 0 and 1")

  # Error message names the calling argument
  expect_error(bvarnet:::.ci_probs(2, "ci_width"), "`ci_width`")
})


test_that("build_summary_table rejects wrong number of columns", {
  draws <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)  # 5 columns

  expect_error(
    bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "T"),
    # should fail: nr*nc = 2*2 = 4 != 5
    "ncol_draws == nr"
  )
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 build_Z()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("build_Z returns empty matrix when no RE specified", {
  X <- matrix(1, 10, 2, dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0, 10, 3, dimnames = list(NULL, paste0("lag1_y_", 1:3)))

  Z <- bvarnet:::build_Z(X, B)

  expect_equal(ncol(Z), 0)
  expect_equal(nrow(Z), 10)
})


test_that("build_Z with re_cols selects correct X columns", {
  X <- matrix(rnorm(30), 10, 3,
              dimnames = list(NULL, c("Intercept", "x_1", "x_2")))
  B <- matrix(0, 10, 2)

  Z <- bvarnet:::build_Z(X, B, re_cols = c("x_1", "x_2"))

  expect_equal(ncol(Z), 2)
  expect_equal(Z[, 1], X[, "x_1"])
  expect_equal(Z[, 2], X[, "x_2"])
})


test_that("build_Z with re_temporal includes all B columns", {
  X <- matrix(1, 10, 1, dimnames = list(NULL, "Intercept"))
  B <- matrix(rnorm(20), 10, 2,
              dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))

  Z <- bvarnet:::build_Z(X, B, re_temporal = TRUE)

  expect_equal(ncol(Z), 2)
  expect_equal(Z[, 1], B[, 1])
  expect_equal(Z[, 2], B[, 2])
})


test_that("build_Z rejects invalid re_cols", {
  X <- matrix(1, 10, 2, dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0, 10, 2)

  expect_error(
    bvarnet:::build_Z(X, B, re_cols = "missing_col"),
    "re_cols not found"
  )
})


test_that("build_Z creates Intercept column when absent from X (ordinal case)", {
  # Ordinal X has no Intercept column (absorbed by kappa cutpoints)
  X <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("x_1", "x_2")))
  B <- matrix(0, 10, 3)

  Z <- bvarnet:::build_Z(X, B, re_cols = "Intercept")

  expect_equal(ncol(Z), 1)
  expect_equal(nrow(Z), 10)
  expect_true(all(Z[, 1] == 1))
  expect_equal(colnames(Z), "Intercept")
})


test_that("build_Z creates Intercept + selects other re_cols when Intercept absent", {
  X <- matrix(rnorm(20), 10, 2, dimnames = list(NULL, c("x_1", "x_2")))
  B <- matrix(0, 10, 3)

  Z <- bvarnet:::build_Z(X, B, re_cols = c("Intercept", "x_1"))

  expect_equal(ncol(Z), 2)
  expect_true(all(Z[, "Intercept"] == 1))
  expect_equal(Z[, "x_1"], X[, "x_1"])
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 make_term_matrix()
# ═══════════════════════════════════════════════════════════════════════════════

test_that("make_term_matrix creates simple interaction from X columns", {
  X <- matrix(c(rep(1, 5), 1:5, 6:10), 5, 3,
              dimnames = list(NULL, c("Intercept", "x_1", "x_2")))
  B <- matrix(0, 5, 2)

  M <- bvarnet:::make_term_matrix(X, B, c("x_1", "x_2"))

  expect_equal(ncol(M), 1)
  expect_equal(colnames(M), "x_1:x_2")
  expect_equal(M[, 1], X[, "x_1"] * X[, "x_2"])
})


test_that("make_term_matrix with 'lag' expands across B columns", {
  X <- matrix(c(rep(1, 5), 1:5), 5, 2,
              dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(rnorm(10), 5, 2,
              dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))

  M <- bvarnet:::make_term_matrix(X, B, c("lag", "x_1"))

  expect_equal(ncol(M), 2)
  expect_equal(colnames(M), c("lag1_y_1:x_1", "lag1_y_2:x_1"))
  expect_equal(M[, 1], B[, 1] * X[, "x_1"])
})


# ═══════════════════════════════════════════════════════════════════════════════
# §N — print.bvarnet() priors display
# ═══════════════════════════════════════════════════════════════════════════════

test_that("print.bvarnet() shows 'all defaults' when no priors are user-set", {
  obj <- make_mock_bvarnet()
  out <- capture.output(print(obj))
  expect_true(any(grepl("all defaults", out)))
})

test_that("print.bvarnet() includes 'phi' in priors summary", {
  obj <- make_mock_bvarnet()
  expect_output(print(obj), "phi")
})

test_that("print.bvarnet() includes 'Half-' for sd_u in priors when model has REs", {
  obj <- make_mock_bvarnet(n_re = 1L)
  expect_output(print(obj), "Half-")
})

test_that("print.bvarnet() omits sd_u from priors when model has no REs", {
  obj <- make_mock_bvarnet()
  out <- capture.output(print(obj))
  priors_line <- out[grepl("Priors", out)]
  expect_false(any(grepl("sd_u", priors_line)))
})

test_that("print.bvarnet() shows 'sigma' in priors for gaussian family", {
  obj <- make_mock_bvarnet(family = "gaussian")
  expect_output(print(obj), "sigma")
})

test_that("print.bvarnet() shows 'kappa' in priors for ordinal family", {
  obj <- make_mock_bvarnet(family = "ordinal")
  expect_output(print(obj), "kappa")
})

test_that("print.bvarnet() omits sigma for bernoulli family", {
  obj <- make_mock_bvarnet(family = "bernoulli")
  out <- capture.output(print(obj))
  expect_false(any(grepl("sigma", out)))
})

test_that("print.bvarnet() shows (default) tag for user-set priors", {
  obj <- make_mock_bvarnet(family = "gaussian")
  obj$priors$beta <- prior("cauchy", 0, 2)
  out <- capture.output(print(obj))
  # beta should NOT have (default) tag
  beta_line <- out[grepl("beta", out)]
  expect_false(any(grepl("\\(default\\)", beta_line)))
  # phi and sigma should have (default) tag
  phi_line <- out[grepl("phi", out)]
  expect_true(any(grepl("\\(default\\)", phi_line)))
  sigma_line <- out[grepl("sigma", out)]
  expect_true(any(grepl("\\(default\\)", sigma_line)))
})

test_that("print.bvarnet() backward compat: old objects without priors_needed still work", {
  obj <- make_mock_bvarnet(family = "gaussian")
  obj$priors_needed <- NULL
  out <- capture.output(print(obj))
  expect_true(any(grepl("beta", out)))
  expect_true(any(grepl("sigma", out)))
})


# ═══════════════════════════════════════════════════════════════════════════════
# .check_sampler_dots() — bvar(...) passthrough to CmdStanR's $sample()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".check_sampler_dots() returns an empty list when no dots are given", {
  expect_identical(.check_sampler_dots(), list())
})

test_that(".check_sampler_dots() passes through genuine sampler arguments", {
  expect_identical(
    .check_sampler_dots(refresh = 0, init = 2, thin = 5),
    list(refresh = 0, init = 2, thin = 5)
  )
})

test_that(".check_sampler_dots() preserves NULL values without dropping names", {
  dots <- .check_sampler_dots(init = NULL)
  expect_named(dots, "init")
  expect_null(dots$init)
})

test_that(".check_sampler_dots() rejects sampler args bvar() sets itself", {
  expect_error(.check_sampler_dots(iter_warmup = 100), "use bvar\\(warmup = ")
  expect_error(.check_sampler_dots(iter_sampling = 100), "use bvar\\(iter = ")
  expect_error(.check_sampler_dots(parallel_chains = 4), "use bvar\\(cores = ")
})

test_that(".check_sampler_dots() rejects CmdStanR's deprecated aliases", {
  # $sample() resolves these by overwriting the modern argument, which would
  # silently override the value bvar() supplied.
  expect_error(.check_sampler_dots(num_warmup = 100), "use bvar\\(warmup = ")
  expect_error(.check_sampler_dots(num_samples = 100), "use bvar\\(iter = ")
  expect_error(.check_sampler_dots(num_chains = 2), "use bvar\\(chains = ")
  expect_error(.check_sampler_dots(num_cores = 2), "use bvar\\(cores = ")
  expect_error(.check_sampler_dots(cores = 2), "use bvar\\(cores = ")
  expect_error(.check_sampler_dots(max_depth = 12), "use bvar\\(max_treedepth = ")
})

test_that(".check_sampler_dots() reports every reserved argument it found", {
  err <- expect_error(
    .check_sampler_dots(iter_warmup = 1, parallel_chains = 2, refresh = 0)
  )
  expect_match(conditionMessage(err), "iter_warmup")
  expect_match(conditionMessage(err), "parallel_chains")
  expect_no_match(conditionMessage(err), "refresh")
})

test_that(".check_sampler_dots() requires all dots to be named", {
  expect_error(.check_sampler_dots(5), "must be named")
  expect_error(.check_sampler_dots(refresh = 0, 5), "must be named")
})

test_that(".check_sampler_dots() rejects duplicated argument names", {
  expect_error(.check_sampler_dots(refresh = 0, refresh = 1), "Duplicated")
})

test_that("every modern reserved name is a real CmdStanR $sample() argument", {
  skip_if_not_installed("cmdstanr")
  sample_args <- names(formals(
    getFromNamespace("CmdStanModel", "cmdstanr")$public_methods$sample
  ))
  expect_length(setdiff(.bvarnet_fixed_sample_args, sample_args), 0L)
  expect_length(
    setdiff(.bvarnet_fixed_sample_args, names(.bvarnet_reserved_sampler_args)),
    0L
  )
})

test_that(".check_sampler_dots() follows CmdStanR's available alias formals", {
  skip_if_not_installed("cmdstanr")
  sample_args <- names(formals(
    getFromNamespace("CmdStanModel", "cmdstanr")$public_methods$sample
  ))

  # $sample() has no '...', so do.call() partial-matches any dot name that
  # isn't an exact formal name against the formals bvar() hasn't already
  # claimed by exact name (e.g. 'num_warm' -> 'num_warmup'). Left unchecked,
  # this silently overrides the value bvar() set for 'warmup'. CmdStanR 1.0
  # removes these deprecated formals, in which case there is no partial match
  # for this helper to intercept and the argument is left for $sample() to
  # reject.
  check_alias <- function(alias, formal, value, replacement) {
    args <- stats::setNames(list(value), alias)
    if (formal %in% sample_args) {
      expect_error(do.call(.check_sampler_dots, args), replacement)
    } else {
      expect_identical(do.call(.check_sampler_dots, args), args)
    }
  }

  check_alias("num_warm", "num_warmup", 50, "use bvar\\(warmup = ")
  check_alias("num_sam", "num_samples", 50, "use bvar\\(iter = ")
  check_alias("num_ch", "num_chains", 2, "use bvar\\(chains = ")
  check_alias("core", "cores", 2, "use bvar\\(cores = ")
  check_alias("max_de", "max_depth", 8, "use bvar\\(max_treedepth = ")
})

test_that(".check_sampler_dots() still allows genuine abbreviations of non-reserved args", {
  skip_if_not_installed("cmdstanr")
  # 'ref' uniquely abbreviates the passthrough-only formal 'refresh', not any
  # reserved name, so it must not be rejected.
  expect_identical(.check_sampler_dots(ref = 0), list(ref = 0))
})

test_that(".check_sampler_dots() degrades to exact-name matching without cmdstanr", {
  testthat::local_mocked_bindings(
    .bvarnet_sample_partial_pool = function() character(0)
  )
  # No partial-match pool available: exact reserved/alias names are still
  # rejected, but an abbreviation can no longer be resolved and silently
  # passes the guard (it will fail later at do.call(), same as before this
  # fix existed).
  expect_error(.check_sampler_dots(num_warmup = 50), "use bvar\\(warmup = ")
  expect_identical(.check_sampler_dots(num_warm = 50), list(num_warm = 50))
})


# ═══════════════════════════════════════════════════════════════════════════════
# .bvarnet_node_dots() — keeps per-node $sample() output files distinct
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_node_dots() is a no-op when output_basename isn't set", {
  dots <- list(refresh = 0, thin = 2)
  expect_identical(.bvarnet_node_dots(dots, 1), dots)
  expect_identical(.bvarnet_node_dots(list(), 3), list())
})

test_that(".bvarnet_node_dots() suffixes output_basename with the node index", {
  dots <- list(output_basename = "myfit", refresh = 0)
  expect_identical(
    .bvarnet_node_dots(dots, 2),
    list(output_basename = "myfit_node2", refresh = 0)
  )
})

test_that(".bvarnet_node_dots() gives every node a distinct output_basename", {
  dots <- list(output_basename = "myfit")
  suffixed <- lapply(1:3, function(j) .bvarnet_node_dots(dots, j)$output_basename)
  expect_length(unique(unlist(suffixed)), 3L)
})
