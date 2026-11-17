test_that("day/beep lag windows match exact scheduled predecessors at every order", {
  d <- expand.grid(beep = 1:6, day = 1:2, id = c(2, 10))
  d$y <- seq_len(nrow(d))
  for (K in 1:3) {
    sd <- suppressMessages(to_stan_data(d, "gaussian", "id", "day", "y",
                                       character(), K = K, beep_col = "beep", save_data = TRUE))
    expect_equal(sd$n_obs, 4 * (6 - K))
    expect_equal(sd$J, 2)
    expect_equal(sd$index_obs, sd$data_used[sd$row_map, c("id", "day", "beep")], ignore_attr = TRUE)
    for (r in seq_len(sd$n_obs)) for (lag in seq_len(K)) {
      ix <- sd$index_obs[r, ]
      expected <- d$y[d$id == ix$id & d$day == ix$day & d$beep == ix$beep - lag]
      expect_equal(unname(sd$B[r, lag]), expected)
    }
    expect_false("skip_lag" %in% names(sd$design_spec))
    expect_false("gap_action" %in% names(sd$design_spec))
  }
})

test_that("a boundary is enforced even with consecutive numbers across days", {
  d <- data.frame(id = 1, day = c(1,1,1,2,2,2), beep = 1:6, y = 11:16)
  sd <- suppressMessages(to_stan_data(d, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "beep"))
  expect_equal(sd$index_obs$beep, c(2,3,5,6))
  expect_equal(unname(sd$B[,1]), c(11,12,14,15))
  classic <- suppressMessages(to_stan_data(d, "gaussian", "id", "beep", "y", NULL, K = 1))
  expect_equal(classic$time_obs, 2:6)
})

test_that("missing occasions exclude windows without recursively deleting predictors", {
  d <- data.frame(id = 1, day = 1, beep = 1:8, y = c(0,0,3,NA,5,6,7,8), x = 0)
  for (K in 1:3) {
    sd <- suppressMessages(to_stan_data(d, "gaussian", "id", "day", "y", "x", K = K, beep_col = "beep"))
    expect_equal(sd$index_obs$beep, if (K == 1) c(2,3,6,7,8) else if (K == 2) c(3,7,8) else 8)
    absent <- d[!is.na(d$y), ]
    other <- suppressMessages(to_stan_data(absent, "gaussian", "id", "day", "y", "x", K = K, beep_col = "beep"))
    expect_identical(sd$B, other$B)
  }
  # Missing a covariate retains the existing conservative listwise contract.
  d$y[4] <- 4; d$x[4] <- NA_real_
  sd <- suppressMessages(to_stan_data(d, "gaussian", "id", "day", "y", "x", K = 1, beep_col = "beep"))
  expect_equal(sd$index_obs$beep, c(2,3,6,7,8))
  expect_identical(unname(sd$B[1,1]), 0)
})

test_that("typed segment labels and composite duplicate checks preserve metadata", {
  base <- data.frame(id = 1, day = rep(c(1,3), each = 4), beep = rep(1:4, 2), y = 1:8)
  for (days in list(as.Date("2026-01-01") + base$day, as.character(base$day), factor(base$day))) {
    d <- base; d$day <- days
    sd <- suppressMessages(to_stan_data(d, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "beep"))
    expect_equal(sd$time_obs, days[c(2:4,6:8)])
    expect_identical(sd$design_spec$beep_col, "beep")
    expect_error(to_stan_data(rbind(d, d[2,]), "gaussian", "id", "day", "y", NULL,
                              K = 1, beep_col = "beep"), "Duplicate \\(id, day, beep\\)")
  }
})

test_that("invalid keys fail clearly and the removed API cannot reach sampling", {
  d <- data.frame(id = 1, day = 1, beep = 1:4, y = 1:4)
  for (key in c("id", "day", "beep")) {
    bad <- d; bad[[key]][2] <- NA
    expect_error(to_stan_data(bad, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "beep"), "Index column")
  }
  bad <- d; bad$beep[2] <- Inf
  expect_error(to_stan_data(bad, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "beep"), "finite")
  bad$beep[2] <- 2.2
  expect_error(to_stan_data(bad, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "beep"), "integer-valued")
  expect_error(to_stan_data(d, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "missing"), "not found")
  expect_error(to_stan_data(d, "gaussian", "id", "day", "y", NULL, K = 1, beep_col = "day"), "distinct")
  expect_error(to_stan_data(d, "gaussian", "id", "day", "y", NULL, K = 1.5), "positive integer")
  for (value in list(TRUE, FALSE, NULL, NA, "drop"))
    expect_error(bvar(skip_lag = value), "skip_lag.*removed")
  expect_false("skip_lag" %in% names(formals(bvar)))
  expect_false("skip_lag" %in% names(formals(.to_stan_data_shared)))
})

test_that("all families share the same retained windows and original random-effect IDs", {
  d <- expand.grid(beep = 1:5, day = 1:2, id = 1:2)
  d$y1 <- rep(1:2, 10); d$y2 <- rep(2:1, 10); d$x <- seq_len(nrow(d))/10
  d$y1[d$id == 1 & d$day == 1 & d$beep == 3] <- NA
  for (family in c("gaussian", "bernoulli", "ordinal")) {
    x <- d
    if (family == "bernoulli") x[c("y1", "y2")] <- x[c("y1", "y2")] - 1
    sd <- suppressMessages(to_stan_data(x, family, "id", "day", c("y1","y2"), "x", K = 1,
      beep_col = "beep", re_cols = "Intercept", re_temporal = TRUE))
    expect_equal(sd$n_obs, 14)
    expect_equal(sd$n_re, 3)
    expect_equal(sd$Z[,1], rep(1,14))
    expect_equal(sd$id, match(sd$index_obs$id, sd$id_levels))
    expect_equal(unname(sd$Z[,-1]), unname(sd$B))
  }
})

test_that("R-only index metadata never reaches the homogeneous Stan sampler", {
  probe <- bvar
  environment(probe) <- list2env(list(.bvarnet_stan_model = function(...) {
    list(cache_dir = NULL, model = list(sample = function(data, ...) {
      expect_false(any(c("index_obs", "time_obs", "row_map", "n_rows_data",
                         "design_spec", "data_used") %in% names(data)))
      stop("SAMPLER_CAPTURED", call. = FALSE)
    }))
  }), parent = environment(bvar))
  d <- data.frame(id = 1, day = rep(as.Date("2026-11-17"), 5), beep = 1:5, y = 1:5)
  expect_error(suppressMessages(probe(id_col = "id", time_col = "day", beep_col = "beep",
    y_cols = "y", data = d, family = "gaussian")), "SAMPLER_CAPTURED")
})
