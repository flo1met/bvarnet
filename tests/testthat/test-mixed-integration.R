# test-mixed-integration.R — mock-based tests for mixed-family support
# ──────────────────────────────────────────────────────────────────────────────

# ── Mock construction ────────────────────────────────────────────────────────

test_that("make_mock_bvarnet() creates valid mixed GB object", {
  obj <- make_mock_bvarnet(family = c("gaussian", "bernoulli"))
  expect_s3_class(obj, "bvarnet")
  expect_equal(as.character(obj$family), c("gaussian", "bernoulli"))
  expect_equal(names(obj$family), c("y_1", "y_2"))

  par_nms <- dimnames(obj$draws)[[3]]
  # sigma only for node 1 (gaussian)
  expect_true("sigma[1]" %in% par_nms)
  expect_false("sigma[2]" %in% par_nms)
  # no kappa
  expect_false(any(grepl("^kappa\\[", par_nms)))
})

test_that("make_mock_bvarnet() creates valid mixed GO object", {
  obj <- make_mock_bvarnet(family = c("gaussian", "ordinal"))
  par_nms <- dimnames(obj$draws)[[3]]
  # sigma only for node 1

  expect_true("sigma[1]" %in% par_nms)
  expect_false("sigma[2]" %in% par_nms)
  # kappa only for node 2
  expect_true("kappa[2,1]" %in% par_nms)
  expect_true("kappa[2,2]" %in% par_nms)
  expect_false("kappa[1,1]" %in% par_nms)
  # ordinal beta[1,2] is NA sentinel
  expect_true(all(is.na(obj$draws[, , "beta[1,2]"])))
  # gaussian beta[1,1] is real
  expect_false(any(is.na(obj$draws[, , "beta[1,1]"])))
})

test_that("make_mock_bvarnet() creates valid mixed BO object", {
  obj <- make_mock_bvarnet(family = c("bernoulli", "ordinal"))
  par_nms <- dimnames(obj$draws)[[3]]
  expect_false(any(grepl("^sigma\\[", par_nms)))
  expect_true("kappa[2,1]" %in% par_nms)
  expect_true(all(is.na(obj$draws[, , "beta[1,2]"])))
})

# ── extract_param on mixed objects ───────────────────────────────────────────

test_that("extract_param() works on mixed GB mock", {
  obj <- make_mock_bvarnet(family = c("gaussian", "bernoulli"))
  result <- extract_param(obj)
  expect_s3_class(result, "data.frame")

  # sigma rows only for gaussian nodes
  sigma_rows <- result[result$type == "Residual SD", ]
  expect_equal(nrow(sigma_rows), 1L)
  expect_equal(sigma_rows$predictor, "y_1")

  # no kappa rows
  kappa_rows <- result[result$type == "Threshold", ]
  expect_equal(nrow(kappa_rows), 0L)
})

test_that("extract_param() filters ordinal intercept rows in mixed GO mock", {
  obj <- make_mock_bvarnet(family = c("gaussian", "ordinal"))
  result <- extract_param(obj)

  # beta section: ordinal node should NOT have an Intercept row
  intercept_rows <- result[result$type == "Intercept", ]
  ord_intercept <- intercept_rows[intercept_rows$outcome == "y_2", ]
  expect_equal(nrow(ord_intercept), 0L)

  # gaussian node SHOULD have an Intercept row
  gauss_intercept <- intercept_rows[intercept_rows$outcome == "y_1", ]
  expect_equal(nrow(gauss_intercept), 1L)
})

# ── print method on mixed objects ────────────────────────────────────────────

test_that("print.bvarnet() works on mixed objects", {
  obj <- make_mock_bvarnet(family = c("gaussian", "bernoulli"))
  out <- capture.output(print(obj))
  out_str <- paste(out, collapse = "\n")
  expect_true(grepl("gaussian", out_str))
  expect_true(grepl("bernoulli", out_str))
})

# ── sim_var with mixed family ───────────────────────────────────────────────

test_that("sim_var() generates data with mixed GB family", {
  sim <- sim_var(N = 5, T_obs = 20, p = 2, family = c("gaussian", "bernoulli"),
                 seed = 100)
  expect_true(is.data.frame(sim$data))
  expect_equal(as.character(sim$truth$family), c("gaussian", "bernoulli"))
  # sigma: non-NA for gaussian node
  expect_false(is.na(sim$truth$sigma[1]))
  expect_true(is.na(sim$truth$sigma[2]))
})

test_that("sim_var() generates data with mixed GO family", {
  sim <- sim_var(N = 5, T_obs = 20, p = 2, family = c("gaussian", "ordinal"),
                 C = 4, seed = 101)
  expect_equal(sim$truth$C, 4L)
  expect_true(is.numeric(sim$truth$kappa[[2]]))
  expect_null(sim$truth$kappa[[1]])
})

test_that("sim_var() generates data with mixed GBO family (p=3)", {
  sim <- sim_var(N = 5, T_obs = 20, p = 3,
                 family = c("gaussian", "bernoulli", "ordinal"),
                 C = 3, seed = 102)
  expect_equal(length(sim$truth$family), 3L)
  # y columns present
  expect_true(all(c("y_1", "y_2", "y_3") %in% names(sim$data)))
})

test_that("generate_response_node() dispatches correctly", {
  set.seed(99)
  # bernoulli
  y_b <- generate_response_node(0, "bernoulli")
  expect_true(y_b %in% c(0L, 1L))
  # gaussian
  y_g <- generate_response_node(0, "gaussian", sigma = 1)
  expect_true(is.numeric(y_g))
  # ordinal
  y_o <- generate_response_node(0, "ordinal", kappa = c(-1, 0, 1))
  expect_true(y_o %in% 1:4)
})
