with_internal_locations <- function(object) {
  internal <- object$draws[, , grep("^beta\\[", dimnames(object$draws)[[3]]), drop = FALSE]
  dimnames(internal)[[3]] <- sub("^beta", "beta_centered", dimnames(internal)[[3]])
  internal[] <- internal + 50
  object$draws <- abind_simple(object$draws, internal)
  object
}

test_that("public extraction, summaries and BFs ignore centered blocks and column order", {
  obj <- make_mock_bvarnet(family = "gaussian", n_iter = 500L, n_chains = 2L)
  eff <- set_priors(intercept = prior("normal", 1, 2), beta = prior("normal", -.3, .8))
  obj$priors_effective <- rep(list(eff), 2)
  changed <- with_internal_locations(obj)
  changed$draws <- changed$draws[, , rev(seq_len(dim(changed$draws)[3])), drop = FALSE]
  for (block in c("beta", "phi")) {
    actual <- extract_draws(changed, block)
    expected <- extract_draws(obj, block)
    expect_equal(actual[, colnames(expected), drop = FALSE], expected)
  }
  a <- extract_param(obj); b <- extract_param(changed)
  key <- function(x) paste(x$type, x$predictor, x$outcome)
  expect_equal(b[match(key(a), key(b)), ], a, ignore_attr = TRUE)
  expect_false(any(grepl("centered", colnames(extract_draws(changed)))))
  expect_equal(.extract_draws_raw(changed, "beta[1,1]"), .extract_draws_raw(obj, "beta[1,1]"))
  for (params in list("beta[1,1]", c("beta[1,1]", "beta[2,1]"))) {
    before <- savage_dickey(obj, params, null_value = .2, method = "auto")
    after <- savage_dickey(changed, params, null_value = .2, method = "auto")
    expect_equal(after, before)
    expected_prior <- dnorm(.2, 1, 2)
    if (length(params) == 2) expected_prior <- expected_prior * dnorm(.2, -.3, .8)
    expect_equal(after$prior_density, expected_prior)
  }
  restored <- unserialize(serialize(changed, NULL))
  expect_equal(extract_param(restored), extract_param(changed))
  tab <- bf_table(changed, type = "intercepts", null_value = .2, round = 10)
  estimates <- extract_param(changed, bayes_factor = TRUE, null_value = .2,
                             type = "Intercept")
  scalar <- tab[tab$type == "Intercept", ]
  expect_equal(scalar$BF10[match(estimates$outcome, scalar$outcome)],
               estimates$BF10, tolerance = 1e-8)
})

test_that("mixed nodes remove internal coordinates from draws and diagnostics", {
  families <- c(a = "gaussian", b = "gaussian", c = "ordinal", d = "ordinal")
  shared <- list(p = 4L, K = 2L, n_fe = 2L, n_re = 0L, n_obs = 12L,
                 Y = cbind(rep(0, 12), rep(0, 12), rep(1:3, 4), rep(1:4, 3)))
  make_node <- function(family, C) {
    nms <- c("lp__", sprintf("phi[%d,1]", 1:8),
             if (family == "ordinal") c("beta[1,1]", sprintf("kappa[1,%d]", seq_len(C - 1)),
                                          sprintf("kappa_centered[1,%d]", seq_len(C - 1)))
             else c("beta[1,1]", "beta[2,1]", "beta_centered[1,1]", "beta_centered[2,1]", "sigma[1]"))
    values <- array(seq_len(40 * 2 * length(nms)) / 100, c(40, 2, length(nms)),
                    dimnames = list(NULL, NULL, nms))
    values[, , grep("centered", nms)] <- 99
    list(draws = function(format) posterior::as_draws_array(values),
         diagnostic_summary = function() data.frame(num_divergent = c(0, 0), num_max_treedepth = c(0, 0), ebfmi = c(1, 1)),
         time = function() list(total = 1), metadata = function() list(), return_codes = function() c(0, 0))
  }
  fits <- Map(make_node, families, c(0, 0, 3, 4))
  obj <- .combine_nodewise_fits(fits, families, shared, set_priors(), rep(list(set_priors()), 4),
                                c("beta", "phi", "kappa", "sigma"), 40, 2)
  nms <- dimnames(obj$draws)[[3]]
  expect_false(anyDuplicated(nms) > 0)
  expect_false(anyDuplicated(obj$convergence$variable) > 0)
  expect_false(any(grepl("centered", c(nms, obj$convergence$variable))))
  expect_true(all(is.na(obj$draws[, , "beta[1,3]"])))
  expect_equal(unname(obj$draws[, , "beta[2,3]"]), unname(unclass(fits[[3]]$draws("array"))[, , "beta[1,1]"]))
  expect_true(all(c("phi[8,4]", "kappa[3,2]", "kappa[4,3]") %in% nms))
  expect_false("kappa[3,3]" %in% nms)
  expected <- posterior::summarise_draws(fits[[1]]$draws("array"), posterior::ess_bulk)
  expect_equal(obj$convergence$ess_bulk[match("beta[1,1]", obj$convergence$variable)],
               expected$`posterior::ess_bulk`[match("beta[1,1]", expected$variable)])
  expect_false("backend" %in% names(obj))
})
