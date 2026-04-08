# ──────────────────────────────────────────────────────────────────────────────
# test-bayes_factor.R — Tests for savage_dickey(), bf_table(), and helpers
# ──────────────────────────────────────────────────────────────────────────────

# ── Layer 1: Analytical ground truth (no MCMC) ───────────────────────────────

test_that("eval_prior_density: normal prior", {
  pr <- prior("normal", loc = 0, scale = 1)
  expect_equal(eval_prior_density(pr, 0), dnorm(0, 0, 1))
  expect_equal(eval_prior_density(pr, 1.5), dnorm(1.5, 0, 1))
})

test_that("eval_prior_density: student_t prior", {
  pr <- prior("student_t", loc = 0, scale = 2, df = 5)
  expected <- dt(0 / 2, df = 5) / 2
  expect_equal(eval_prior_density(pr, 0), expected)
  expected2 <- dt(1.3 / 2, df = 5) / 2
  expect_equal(eval_prior_density(pr, 1.3), expected2)
})

test_that("eval_prior_density: cauchy prior", {
  pr <- prior("cauchy", loc = 0, scale = 1)
  expect_equal(eval_prior_density(pr, 0), dcauchy(0, 0, 1))
})

test_that("eval_joint_prior_density: product of independent marginals", {
  priors <- set_priors()
  types <- c("phi", "phi", "beta")
  null  <- c(0, 0, 0)
  expected <- eval_prior_density(priors$phi, 0) *
              eval_prior_density(priors$phi, 0) *
              eval_prior_density(priors$beta, 0)
  expect_equal(
    eval_joint_prior_density(priors, types, null),
    expected
  )
})

test_that("logspline SDDR matches analytical normal-normal BF", {
  set.seed(1)
  S <- 50000

  # Conjugate normal-normal: prior N(0,1), likelihood N(theta, 1), n=30, ybar=0.15
  # (moderate effect so BF is not too extreme for stable estimation)
  n_obs <- 30; ybar <- 0.15; sigma <- 1; prior_sd <- 1
  post_var  <- 1 / (n_obs / sigma^2 + 1 / prior_sd^2)
  post_mean <- post_var * (n_obs * ybar / sigma^2)
  draws <- rnorm(S, post_mean, sqrt(post_var))

  BF_true <- dnorm(0, post_mean, sqrt(post_var)) / dnorm(0, 0, 1)

  pr <- prior("normal", loc = 0, scale = 1)
  res <- .compute_sddr_logspline(draws, prior = pr, null = 0)
  expect_equal(res$BF01, BF_true, tolerance = 0.10)
})

test_that("MVN SDDR matches analytical bivariate normal-normal BF", {
  set.seed(2)
  S <- 30000

  # Two independent conjugate normals
  prior_sd <- c(1, 0.5)
  n_obs <- 100; ybar <- c(0.3, -0.2); sigma <- c(1, 1)

  post_var  <- 1 / (n_obs / sigma^2 + 1 / prior_sd^2)
  post_mean <- post_var * (n_obs * ybar / sigma^2)

  draws_mat <- cbind(
    rnorm(S, post_mean[1], sqrt(post_var[1])),
    rnorm(S, post_mean[2], sqrt(post_var[2]))
  )
  colnames(draws_mat) <- c("phi[1,1]", "phi[2,1]")

  # Analytical joint posterior density at (0, 0) — independent, so product
  post_den_true <- dnorm(0, post_mean[1], sqrt(post_var[1])) *
                   dnorm(0, post_mean[2], sqrt(post_var[2]))
  prior_den_true <- dnorm(0, 0, prior_sd[1]) * dnorm(0, 0, prior_sd[2])
  BF_true <- post_den_true / prior_den_true

  # Use two phi priors with different scales — need to pass them via prior_list
  priors <- list(phi = prior("normal", 0, 1))
  # Since both use "phi" type, the function uses the same prior for both — but
  # our analytical uses different scales.  So let's use same scale for both.
  draws_mat2 <- cbind(
    rnorm(S, post_mean[1], sqrt(post_var[1])),
    rnorm(S, post_mean[1], sqrt(post_var[1]))  # same marginal
  )
  colnames(draws_mat2) <- c("phi[1,1]", "phi[2,1]")

  # Equal-prior version
  post_var_eq  <- 1 / (n_obs / 1 + 1 / 1)
  post_mean_eq <- post_var_eq * (n_obs * 0.3)
  post_den_eq  <- dnorm(0, post_mean_eq, sqrt(post_var_eq))^2
  prior_den_eq <- dnorm(0, 0, 1)^2
  BF_eq <- post_den_eq / prior_den_eq

  res <- .compute_sddr_mvn(
    draws_mat2,
    prior_list  = list(phi = prior("normal", 0, 1)),
    param_types = c("phi", "phi"),
    null_vec    = c(0, 0)
  )
  expect_equal(res$BF01, BF_eq, tolerance = 0.10)
})


# ── Layer 2: Consistency across methods ──────────────────────────────────────

test_that("logspline and MVN give consistent single-parameter BFs", {
  set.seed(3)
  S <- 50000
  n_obs <- 30; ybar <- 0.15; prior_sd <- 1
  post_var  <- 1 / (n_obs + 1)
  post_mean <- post_var * (n_obs * ybar)
  draws <- rnorm(S, post_mean, sqrt(post_var))

  pr <- prior("normal", 0, 1)

  res_ls <- .compute_sddr_logspline(draws, prior = pr, null = 0)

  draws_mat <- matrix(draws, ncol = 1)
  colnames(draws_mat) <- "phi[1,1]"
  res_mvn <- .compute_sddr_mvn(draws_mat,
                                prior_list = list(phi = pr),
                                param_types = "phi",
                                null_vec = 0)

  expect_equal(res_ls$BF01, res_mvn$BF01, tolerance = 0.15)
})


# ── get_phi_indices tests ────────────────────────────────────────────────────

test_that("get_phi_indices returns correct AR indices for p=3, K=1", {
  sd <- list(p = 3L, K = 1L)
  ar <- get_phi_indices(sd, lag = 1, effect = "ar")
  expect_equal(ar, c("phi[1,1]", "phi[2,2]", "phi[3,3]"))
})

test_that("get_phi_indices returns correct CL indices for p=3, K=1", {
  sd <- list(p = 3L, K = 1L)
  cl <- get_phi_indices(sd, lag = 1, effect = "cl")
  expect_length(cl, 6)  # p*(p-1) = 6
  expect_false("phi[1,1]" %in% cl)
  expect_false("phi[2,2]" %in% cl)
  expect_false("phi[3,3]" %in% cl)
  expect_true("phi[2,1]" %in% cl)
  expect_true("phi[1,2]" %in% cl)
})

test_that("get_phi_indices returns all for p=2, K=1", {
  sd <- list(p = 2L, K = 1L)
  all <- get_phi_indices(sd, lag = 1, effect = "all")
  expect_length(all, 4)  # p*p = 4
})

test_that("get_phi_indices lag=2 uses correct offset", {
  sd <- list(p = 2L, K = 2L)
  ar2 <- get_phi_indices(sd, lag = 2, effect = "ar")
  expect_equal(ar2, c("phi[3,1]", "phi[4,2]"))
})

test_that("get_phi_indices errors on invalid lag", {
  sd <- list(p = 2L, K = 1L)
  expect_error(get_phi_indices(sd, lag = 2, effect = "ar"))
})


# ── get_beta_indices tests ───────────────────────────────────────────────────

test_that("get_beta_indices returns intercepts correctly", {
  X <- matrix(0, 10, 2, dimnames = list(NULL, c("Intercept", "x_1")))
  sd <- list(p = 3L, n_fe = 2L, X = X)
  ic <- get_beta_indices(sd, type = "intercepts")
  expect_equal(ic, c("beta[1,1]", "beta[1,2]", "beta[1,3]"))
})

test_that("get_beta_indices returns FE (non-intercept) correctly", {
  sd <- list(p = 2L, n_fe = 3L)
  fe <- get_beta_indices(sd, type = "fe")
  expect_length(fe, 4)  # (3-1) * 2 = 4
  expect_true("beta[2,1]" %in% fe)
  expect_true("beta[3,2]" %in% fe)
  expect_false("beta[1,1]" %in% fe)
})


# ── savage_dickey() with mock bvarnet object ─────────────────────────────────

test_that("savage_dickey errors on missing param names", {
  mock <- make_mock_bvarnet("gaussian")
  expect_error(
    savage_dickey(mock, params = "phi[99,99]"),
    "not found in draws"
  )
})

test_that("savage_dickey errors on any half-prior param", {
  mock <- make_mock_bvarnet("gaussian")
  expect_error(
    savage_dickey(mock, params = "sigma[1]", null_value = 0),
    "Cannot compute SDDR"
  )
  # Also blocked at positive interior null values
  expect_error(
    savage_dickey(mock, params = "sigma[1]", null_value = 0.5),
    "Cannot compute SDDR"
  )
})

test_that("savage_dickey returns correct structure for univariate", {
  mock <- make_mock_bvarnet("gaussian")
  res <- savage_dickey(mock, params = "phi[1,1]", null_value = 0)
  expect_true(is.list(res))
  expect_named(res, c("BF01", "BF10", "log_BF01", "post_density",
                       "prior_density", "method", "params", "null_value"))
  expect_equal(res$method, "logspline")
  expect_equal(res$BF01 * res$BF10, 1)
  expect_equal(res$log_BF01, log(res$BF01))
})

test_that("savage_dickey returns correct structure for joint MVN", {
  mock <- make_mock_bvarnet("gaussian")
  res <- savage_dickey(mock, params = c("phi[1,1]", "phi[2,1]"),
                       null_value = 0, method = "mvn")
  expect_true(is.list(res))
  expect_equal(res$method, "mvn")
  expect_true(is.finite(res$BF01))
  expect_equal(res$BF01 * res$BF10, 1)
})



# ── bf_table() output shape tests ────────────────────────────────────────────

test_that("bf_table returns correct data frame shape for ar + cl", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = c("ar", "cl"))

  expect_s3_class(res, "data.frame")
  expect_true(all(c("type", "predictor", "outcome", "BF01", "BF10",
                     "log_BF01", "post_density", "prior_density", "method")
                  %in% names(res)))

  # p=2, K=1 → ar: 2 + 1 joint + cl: 2 + 1 joint = 6
  expect_equal(nrow(res), 6)
  expect_true(all(res$BF01 > 0))
  expect_true(any(grepl("joint", res$type)))
})

test_that("bf_table returns correct rows for ar type (p=2, K=1)", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = "ar")

  # p=2: 2 AR params + 1 joint row = 3
  expect_equal(nrow(res), 3)
  expect_true(any(grepl("joint", res$type)))
})

test_that("bf_table returns correct rows for cl type (p=2, K=1)", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = "cl")

  # p=2: 2 CL params + 1 joint row = 3
  expect_equal(nrow(res), 3)
})

test_that("bf_table returns correct rows for intercepts", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = "intercepts")

  # p=2: 2 intercept params + 1 joint row = 3
  expect_equal(nrow(res), 3)
})

test_that("bf_table handles multiple types", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = c("ar", "cl"))

  # ar: 2 + 1 + cl: 2 + 1 = 6
  expect_equal(nrow(res), 6)
})

test_that("bf_table joint row BF01 is positive and finite", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = c("ar", "cl"))
  joint_rows <- res[grepl("joint", res$type), ]
  expect_equal(nrow(joint_rows), 2)  # one for AR, one for CL
  expect_true(all(is.finite(joint_rows$BF01)))
  expect_true(all(joint_rows$BF01 > 0))
})


# ── extract_param integration ────────────────────────────────────────────────

test_that("extract_param with bayes_factor = TRUE adds BF columns", {
  mock <- make_mock_bvarnet("gaussian")
  res <- extract_param(mock, bayes_factor = TRUE)
  expect_true("BF01" %in% names(res))
  expect_true("BF10" %in% names(res))

  # BFs should be computed for Intercept, Fixed Effect, Autoregressive, Cross-lagged rows
  bf_rows <- res[res$type %in% c("Intercept", "Fixed Effect", "Autoregressive", "Cross-lagged"), ]
  expect_true(all(!is.na(bf_rows$BF01)))
  expect_true(all(bf_rows$BF01 > 0))
})

test_that("extract_param with bayes_factor = FALSE has no BF columns", {
  mock <- make_mock_bvarnet("gaussian")
  res <- extract_param(mock, bayes_factor = FALSE)
  expect_false("BF01" %in% names(res))
})


# ── Edge cases ───────────────────────────────────────────────────────────────

test_that("savage_dickey warns on very few draws", {
  mock <- make_mock_bvarnet("gaussian", n_iter = 50L, n_chains = 2L)
  # 50 draws * 2 chains = 100 < 1000
  expect_warning(
    savage_dickey(mock, params = "phi[1,1]", null_value = 0),
    "Very few posterior draws"
  )
})

test_that("p = 1 case: AR and phi produce same params", {
  # Create a p=1 mock
  sd <- list(p = 1L, K = 1L)
  ar <- get_phi_indices(sd, lag = 1, effect = "ar")
  all <- get_phi_indices(sd, lag = 1, effect = "all")
  cl <- get_phi_indices(sd, lag = 1, effect = "cl")
  expect_equal(ar, all)
  expect_length(cl, 0)
})


# ── get_beta_indices_by_predictor tests ──────────────────────────────────────

test_that("get_beta_indices_by_predictor groups by predictor row", {
  X <- matrix(0, 10, 4,
              dimnames = list(NULL, c("Intercept", "x_1", "x_2", "x_1:x_2")))
  sd <- list(p = 3L, n_fe = 4L, X = X)
  by_pred <- get_beta_indices_by_predictor(sd, type = "fe")

  # 3 predictor rows (rows 2, 3, 4); each should have p = 3 beta names
  expect_length(by_pred, 3)
  expect_named(by_pred, c("x_1", "x_2", "x_1:x_2"))
  expect_equal(by_pred[["x_1"]], c("beta[2,1]", "beta[2,2]", "beta[2,3]"))
  expect_equal(by_pred[["x_2"]], c("beta[3,1]", "beta[3,2]", "beta[3,3]"))
})

test_that("get_beta_indices_by_predictor handles intercepts", {
  X <- matrix(0, 10, 2, dimnames = list(NULL, c("Intercept", "x_1")))
  sd <- list(p = 2L, n_fe = 2L, X = X)
  by_pred <- get_beta_indices_by_predictor(sd, type = "intercepts")

  expect_length(by_pred, 1)
  expect_named(by_pred, "Intercept")
  expect_equal(by_pred[["Intercept"]], c("beta[1,1]", "beta[1,2]"))
})

test_that("get_beta_indices_by_predictor errors for ordinal intercepts", {
  X <- matrix(0, 10, 2, dimnames = list(NULL, c("x_1", "x_2")))
  sd <- list(p = 2L, n_fe = 2L, X = X, C = 5L)
  expect_error(
    get_beta_indices_by_predictor(sd, type = "intercepts"),
    "not valid for ordinal"
  )
})

test_that("get_beta_indices errors for ordinal intercepts", {
  X <- matrix(0, 10, 2, dimnames = list(NULL, c("x_1", "x_2")))
  sd <- list(p = 2L, n_fe = 2L, X = X, C = 5L)
  expect_error(
    get_beta_indices(sd, type = "intercepts"),
    "not valid for ordinal"
  )
})


# ── bf_table three-level FE tests ────────────────────────────────────────────

test_that("bf_table fe: p=2, 2 predictors → 4 cell + 2 joint + 1 joint all = 7", {
  # Build mock with n_fe = 3 (Intercept + x_1 + x_2), p = 2
  p <- 2L; K <- 1L; n_fe <- 3L; n_re <- 0L
  beta_nm <- sprintf("beta[%d,%d]", rep(1:n_fe, times = p), rep(1:p, each = n_fe))
  phi_nm  <- c("phi[1,1]", "phi[2,1]", "phi[1,2]", "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm, "sigma[1]", "sigma[2]")

  set.seed(99L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, c("Intercept", "x_1", "x_2")))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "fe")
  expect_equal(nrow(res), 7)

  # Per-cell rows
  cell_rows <- res[res$type == "Fixed Effect", ]
  expect_equal(nrow(cell_rows), 4)  # 2 predictors * 2 outcomes

  # Per-predictor joint rows
  joint_rows <- res[res$type == "Fixed Effect (joint)", ]
  expect_equal(nrow(joint_rows), 2)
  expect_true(all(joint_rows$method == "mvn"))

  # Global joint-all row
  all_rows <- res[res$type == "Fixed Effect (joint all)", ]
  expect_equal(nrow(all_rows), 1)
  expect_equal(all_rows$predictor, "all_fe")
})

test_that("bf_table fe: p=1, 2 predictors → 2 cell + 0 joint + 1 joint all = 3", {
  # p=1 means Phase B is skipped (per-cell is already the joint)
  p <- 1L; K <- 1L; n_fe <- 3L; n_re <- 0L
  par_nms <- c("beta[1,1]", "beta[2,1]", "beta[3,1]",
               "phi[1,1]", "sigma[1]")

  set.seed(98L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, "y_1"))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, c("Intercept", "x_1", "x_2")))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, "lag1_y_1"))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "fe")
  # 2 per-cell + 0 per-predictor joint (p=1) + 1 joint-all (2 predictors) = 3
  expect_equal(nrow(res), 3)
  expect_equal(sum(res$type == "Fixed Effect (joint)"), 0)
  expect_equal(sum(res$type == "Fixed Effect (joint all)"), 1)
})

test_that("bf_table intercepts: p=2, 1 predictor row → 2 + 1 + 0 = 3", {
  # Phase C skipped because only 1 predictor row
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = "intercepts")
  expect_equal(nrow(res), 3)
  expect_equal(sum(res$type == "Intercept (joint)"), 1)
  expect_equal(sum(grepl("joint all", res$type)), 0)
})


# ── lag_fe type tests ────────────────────────────────────────────────────────

test_that("bf_table lag_fe: p=2, K=1, 1 interaction → 2 rows", {
  # Build a mock with lag×x_1 interaction
  p <- 2L; K <- 1L
  # X: Intercept, x_1, lag1_y_1:x_1, lag1_y_2:x_1
  n_fe <- 4L; n_re <- 0L
  fe_names <- c("Intercept", "x_1", "lag1_y_1:x_1", "lag1_y_2:x_1")
  beta_nm <- sprintf("beta[%d,%d]", rep(1:n_fe, times = p), rep(1:p, each = n_fe))
  phi_nm  <- c("phi[1,1]", "phi[2,1]", "phi[1,2]", "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm, "sigma[1]", "sigma[2]")

  set.seed(97L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, fe_names))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(
      p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z,
      fe_interaction_terms    = list(c("lag", "x_1")),
      fe_interaction_colnames = c("lag1_y_1:x_1", "lag1_y_2:x_1")
    ),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "lag_fe")
  # 1 per-lag joint + 1 omnibus joint = 2
  expect_equal(nrow(res), 2)
  expect_true("Lag Interaction (per lag)" %in% res$type)
  expect_true("Lag Interaction (joint)"   %in% res$type)
  expect_true(all(is.finite(res$BF01) & res$BF01 > 0))
})

test_that("bf_table lag_fe: p=2, K=2 → 3 rows (2 per-lag + 1 omnibus)", {
  p <- 2L; K <- 2L
  fe_names <- c("Intercept", "x_1",
                "lag1_y_1:x_1", "lag1_y_2:x_1",
                "lag2_y_1:x_1", "lag2_y_2:x_1")
  n_fe <- length(fe_names); n_re <- 0L
  beta_nm <- sprintf("beta[%d,%d]", rep(1:n_fe, times = p), rep(1:p, each = n_fe))
  phi_nm  <- sprintf("phi[%d,%d]",
                      rep(1:(p * K), times = p),
                      rep(1:p, each = p * K))
  par_nms <- c(beta_nm, phi_nm, "sigma[1]", "sigma[2]")

  set.seed(96L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, fe_names))
  B <- matrix(0, 10, p * K, dimnames = list(NULL,
      c("lag1_y_1", "lag1_y_2", "lag2_y_1", "lag2_y_2")))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(
      p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z,
      fe_interaction_terms    = list(c("lag", "x_1")),
      fe_interaction_colnames = c("lag1_y_1:x_1", "lag1_y_2:x_1",
                                  "lag2_y_1:x_1", "lag2_y_2:x_1")
    ),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "lag_fe")
  expect_equal(nrow(res), 3)
  per_lag <- res[res$type == "Lag Interaction (per lag)", ]
  expect_equal(nrow(per_lag), 2)
  omnibus <- res[res$type == "Lag Interaction (joint)", ]
  expect_equal(nrow(omnibus), 1)
})

test_that("bf_table lag_fe errors when no lag interactions present", {
  mock <- make_mock_bvarnet("gaussian")
  expect_error(
    bf_table(mock, type = "lag_fe"),
    "No lag interaction columns"
  )
})




# ── to_stan_data metadata persistence ────────────────────────────────────────

test_that("to_stan_data stores fe_interaction_terms and colnames", {
  df <- make_test_df(N = 3, T_obs = 20, p = 2, q = 1, family = "bernoulli")
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), "x_1",
                     fe_interactions = list(c("lag", "x_1")), K = 1)

  expect_true(!is.null(sd$fe_interaction_terms))
  expect_true(is.list(sd$fe_interaction_terms))
  expect_length(sd$fe_interaction_terms, 1)
  expect_equal(sd$fe_interaction_terms[[1]], c("lag", "x_1"))

  expect_true(!is.null(sd$fe_interaction_colnames))
  expect_true(all(grepl(":x_1$", sd$fe_interaction_colnames)))
})

test_that("to_stan_data sets empty metadata when no fe_interactions", {
  df <- make_test_df(N = 3, T_obs = 20, p = 2, q = 1, family = "bernoulli")
  sd <- to_stan_data(df, "bernoulli", "id", "t",
                     paste0("y_", 1:2), "x_1", K = 1)
  expect_equal(sd$fe_interaction_terms, list())
  expect_null(sd$fe_interaction_colnames)
})


# ── bf_table() temporal joint BF ─────────────────────────────────────────────

test_that("bf_table temporal returns 3 rows (joint + AR + CL) for p=2, K=1", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = "temporal")

  expect_s3_class(res, "data.frame")
  # 1 combined + 1 AR + 1 CL = 3
  expect_equal(nrow(res), 3)
  expect_true("Temporal (joint)"    %in% res$type)
  expect_true("Temporal AR (joint)" %in% res$type)
  expect_true("Temporal CL (joint)" %in% res$type)

  combined <- res[res$type == "Temporal (joint)", ]
  expect_equal(combined$predictor, "all_phi")
  expect_equal(combined$outcome, "\u2014")
  expect_true(combined$BF01 > 0)
  expect_true(is.finite(combined$BF01))
  expect_equal(combined$BF01 * combined$BF10, 1)
  expect_equal(combined$log_BF01, log(combined$BF01))
  expect_equal(combined$method, "mvn")

  ar_row <- res[res$type == "Temporal AR (joint)", ]
  expect_equal(ar_row$predictor, "all_ar")
  expect_true(ar_row$BF01 > 0 & is.finite(ar_row$BF01))

  cl_row <- res[res$type == "Temporal CL (joint)", ]
  expect_equal(cl_row$predictor, "all_cl")
  expect_true(cl_row$BF01 > 0 & is.finite(cl_row$BF01))
})

test_that("bf_table temporal collects all phi params (p=2, K=1)", {
  mock <- make_mock_bvarnet("gaussian")
  sd   <- mock$standata
  # p=2, K=1 → 4 phi params total (2 AR + 2 CL)
  all_phi <- get_phi_indices(sd, lag = 1, effect = "all")
  expect_length(all_phi, 4)

  res <- bf_table(mock, type = "temporal")
  combined <- res[res$type == "Temporal (joint)", ]
  # The method should be "mvn" since 4 > 1
  expect_equal(combined$method, "mvn")
})

test_that("bf_table temporal can be combined with ar and cl", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock, type = c("ar", "cl", "temporal"))

  # p=2, K=1: ar(2+1) + cl(2+1) + temporal(3: combined+AR+CL) = 9
  expect_equal(nrow(res), 9)
  expect_true("Temporal (joint)" %in% res$type)
  expect_true("Temporal AR (joint)" %in% res$type)
  expect_true("Temporal CL (joint)" %in% res$type)
  expect_true(any(grepl("Autoregressive", res$type)))
  expect_true(any(grepl("Cross-lagged", res$type)))
})

test_that("bf_table temporal works for all families", {
  for (fam in c("gaussian", "bernoulli", "ordinal")) {
    mock <- make_mock_bvarnet(fam)
    res  <- bf_table(mock, type = "temporal")
    # 3 rows: combined + AR + CL
    expect_equal(nrow(res), 3, info = paste("family:", fam))
    expect_true(all(is.finite(res$BF01)), info = paste("family:", fam))
  }
})

test_that("bf_table temporal without lag interactions emits 3 rows only", {
  mock <- make_mock_bvarnet("gaussian")
  res  <- bf_table(mock, type = "temporal")
  # No lag interactions → combined + AR + CL = 3
  expect_equal(nrow(res), 3)
  expect_setequal(res$type,
    c("Temporal (joint)", "Temporal AR (joint)", "Temporal CL (joint)"))
  expect_false("Temporal Interaction (joint)" %in% res$type)
  expect_false("Temporal + Interactions (joint)" %in% res$type)
})

test_that("bf_table temporal with lag interaction emits interaction rows", {
  p <- 2L; K <- 1L
  fe_names <- c("Intercept", "x_1", "lag1_y_1:x_1", "lag1_y_2:x_1")
  n_fe <- length(fe_names); n_re <- 0L
  beta_nm <- sprintf("beta[%d,%d]", rep(1:n_fe, times = p), rep(1:p, each = n_fe))
  phi_nm  <- c("phi[1,1]", "phi[2,1]", "phi[1,2]", "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm, "sigma[1]", "sigma[2]")

  set.seed(123L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, fe_names))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(
      p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z,
      fe_interaction_terms    = list(c("lag", "x_1")),
      fe_interaction_colnames = c("lag1_y_1:x_1", "lag1_y_2:x_1")
    ),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "temporal")

  # 1 combined + 1 AR + 1 CL + 1 per-term interaction + 1 ARx + 1 CLx + 1 omnibus = 7
  expect_equal(nrow(res), 7)
  expect_true("Temporal (joint)" %in% res$type)
  expect_true("Temporal AR (joint)" %in% res$type)
  expect_true("Temporal CL (joint)" %in% res$type)
  expect_true("Temporal Interaction (joint)" %in% res$type)
  expect_true("Temporal AR × Interaction (joint)" %in% res$type)
  expect_true("Temporal CL × Interaction (joint)" %in% res$type)
  expect_true("Temporal + Interactions (joint)" %in% res$type)
  expect_equal(res$predictor[res$type == "Temporal (joint)"], "all_phi")
  expect_equal(res$predictor[res$type == "Temporal AR (joint)"], "all_ar")
  expect_equal(res$predictor[res$type == "Temporal CL (joint)"], "all_cl")
  expect_equal(res$predictor[res$type == "Temporal Interaction (joint)"], "x_1")
  expect_equal(res$predictor[res$type == "Temporal AR × Interaction (joint)"], "x_1")
  expect_equal(res$predictor[res$type == "Temporal CL × Interaction (joint)"], "x_1")
  expect_equal(res$predictor[res$type == "Temporal + Interactions (joint)"],
               "all_temporal")
  expect_true(all(is.finite(res$BF01) & res$BF01 > 0))
})

test_that("bf_table temporal with 2 interaction terms emits per-term rows", {
  p <- 2L; K <- 1L
  fe_names <- c("Intercept", "x_1", "x_2",
                "lag1_y_1:x_1", "lag1_y_2:x_1",
                "lag1_y_1:x_2", "lag1_y_2:x_2")
  n_fe <- length(fe_names); n_re <- 0L
  beta_nm <- sprintf("beta[%d,%d]", rep(1:n_fe, times = p), rep(1:p, each = n_fe))
  phi_nm  <- c("phi[1,1]", "phi[2,1]", "phi[1,2]", "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm, "sigma[1]", "sigma[2]")

  set.seed(456L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, fe_names))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(
      p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z,
      fe_interaction_terms    = list(c("lag", "x_1"), c("lag", "x_2")),
      fe_interaction_colnames = c("lag1_y_1:x_1", "lag1_y_2:x_1",
                                  "lag1_y_1:x_2", "lag1_y_2:x_2")
    ),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "temporal")

  # 1 combined + 1 AR + 1 CL
  # + 2 per-term interaction + 2 AR×interaction + 2 CL×interaction
  # + 1 omnibus = 10
  expect_equal(nrow(res), 10)
  int_rows <- res[res$type == "Temporal Interaction (joint)", ]
  expect_equal(nrow(int_rows), 2)
  expect_setequal(int_rows$predictor, c("x_1", "x_2"))
  expect_true("Temporal AR (joint)" %in% res$type)
  expect_true("Temporal CL (joint)" %in% res$type)
  arx_rows <- res[res$type == "Temporal AR \u00d7 Interaction (joint)", ]
  expect_equal(nrow(arx_rows), 2)
  expect_setequal(arx_rows$predictor, c("x_1", "x_2"))
  clx_rows <- res[res$type == "Temporal CL \u00d7 Interaction (joint)", ]
  expect_equal(nrow(clx_rows), 2)
  expect_setequal(clx_rows$predictor, c("x_1", "x_2"))
})

test_that("bf_table temporal p=1: no CL row emitted", {
  # p=1 → only 1 AR param, 0 CL params → combined + AR only = 2
  p <- 1L; K <- 1L; n_fe <- 2L; n_re <- 0L
  par_nms <- c("beta[1,1]", "beta[2,1]", "phi[1,1]", "sigma[1]")

  set.seed(789L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, "y_1"))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, c("Intercept", "x_1")))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, "lag1_y_1"))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(p = p, K = K, n_fe = n_fe, n_re = n_re,
                     Y = Y, X = X, B = B, Z = Z,
                     fe_interaction_terms = list()),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock, type = "temporal")
  expect_equal(nrow(res), 2)
  expect_true("Temporal (joint)"    %in% res$type)
  expect_true("Temporal AR (joint)" %in% res$type)
  expect_false("Temporal CL (joint)" %in% res$type)
  # combined and AR should have logspline since only 1 param each
  expect_true(all(res$method == "logspline"))
})


# ── bf_table() type = "all" auto-detection ───────────────────────────────────

test_that("bf_table type='all' works for gaussian (no interactions)", {
  mock <- make_mock_bvarnet("gaussian")
  res <- bf_table(mock)  # default type = "all"

  expect_s3_class(res, "data.frame")
  # intercepts + ar + cl + fe + temporal (no lag_fe)
  # intercepts: 2 + 1 joint = 3
  # ar: 2 + 1 joint = 3
  # cl: 2 + 1 joint = 3
  # fe: 2 + 1 joint = 3 (1 predictor, p=2 → joint; no global joint-all since 1 pred)
  # temporal: 3 (combined + AR + CL)
  expect_true(any(grepl("Intercept", res$type)))
  expect_true(any(grepl("Autoregressive", res$type)))
  expect_true(any(grepl("Cross-lagged", res$type)))
  expect_true(any(grepl("Fixed Effect", res$type)))
  expect_true(any(grepl("Temporal", res$type)))
  # No lag interaction rows
  expect_false(any(grepl("Lag Interaction", res$type)))
  expect_true(all(res$BF01 > 0))
})

test_that("bf_table type='all' skips intercepts for ordinal", {
  mock <- make_mock_bvarnet("ordinal")
  res <- bf_table(mock)
  expect_false(any(grepl("^Intercept", res$type)))
  expect_true(any(grepl("Autoregressive", res$type)))
})

test_that("bf_table type='all' includes lag_fe when interactions exist", {
  p <- 2L; K <- 1L
  fe_names <- c("Intercept", "x_1", "lag1_y_1:x_1", "lag1_y_2:x_1")
  n_fe <- length(fe_names); n_re <- 0L
  beta_nm <- sprintf("beta[%d,%d]", rep(1:n_fe, times = p), rep(1:p, each = n_fe))
  phi_nm  <- c("phi[1,1]", "phi[2,1]", "phi[1,2]", "phi[2,2]")
  par_nms <- c(beta_nm, phi_nm, "sigma[1]", "sigma[2]")

  set.seed(111L)
  draws <- array(rnorm(40 * 2 * length(par_nms)),
                 dim = c(40L, 2L, length(par_nms)),
                 dimnames = list(NULL, NULL, par_nms))

  Y <- matrix(0, 10, p, dimnames = list(NULL, c("y_1", "y_2")))
  X <- matrix(0, 10, n_fe, dimnames = list(NULL, fe_names))
  B <- matrix(0, 10, p * K, dimnames = list(NULL, c("lag1_y_1", "lag1_y_2")))
  Z <- matrix(0, 10, 0)

  mock <- structure(list(
    draws     = draws,
    standata  = list(
      p = p, K = K, n_fe = n_fe, n_re = n_re, Y = Y, X = X, B = B, Z = Z,
      fe_interaction_terms    = list(c("lag", "x_1")),
      fe_interaction_colnames = c("lag1_y_1:x_1", "lag1_y_2:x_1")
    ),
    priors    = set_priors(),
    family    = "gaussian"
  ), class = "bvarnet")

  res <- bf_table(mock)
  expect_true(any(grepl("Lag Interaction", res$type)))
  expect_true(any(grepl("Temporal", res$type)))
  expect_true(any(grepl("Intercept", res$type)))
  expect_true(all(res$BF01 > 0))
})
