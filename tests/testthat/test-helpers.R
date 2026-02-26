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
                     "mean", "median", "q5", "q95") %in% names(tab)))
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


test_that("build_summary_table q5 <= median <= q95", {
  set.seed(501)
  draws <- matrix(rnorm(500 * 4), nrow = 500, ncol = 4)
  tab <- bvarnet:::build_summary_table(draws, c("a", "b"), c("c", "d"), "T")

  expect_true(all(tab$q5 <= tab$median))
  expect_true(all(tab$median <= tab$q95))
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
