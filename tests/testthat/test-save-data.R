test_that("save_data requires a single nonmissing logical", {
  expect_identical(formals(bvar)$save_data, TRUE)
  for (value in list(NULL, NA, logical(), c(TRUE, FALSE), 1, "FALSE"))
    expect_error(bvar(save_data = value), "save_data must be TRUE or FALSE")
})

test_that("data-free model metadata preserves posterior operations", {
  # These comparisons test data retention, not BF precision from short runs.
  bf <- function(...) suppressWarnings(bf_table(...))
  for (family in list("gaussian", "bernoulli", "ordinal",
                      c("gaussian", "bernoulli"), c("gaussian", "ordinal"))) {
    original <- make_mock_bvarnet(family = family, n_re = 2L, n_iter = 300L)
    original$priors_effective <- rep(list(original$priors), original$standata$p)
    stripped <- original
    stripped$model_info <- .bvarnet_model_info(original$standata)
    stripped$standata <- NULL
    stripped <- unserialize(serialize(stripped, NULL))

    expect_false(any(c("Y", "X", "B", "Z", "id", "id_levels", "index_obs",
                       "time_obs", "row_map", "data_used", "x_center_means") %in%
                     names(stripped$model_info)))
    expect_equal(extract_draws(stripped), extract_draws(original))
    expect_equal(extract_param(stripped), extract_param(original))
    expect_equal(summary(stripped), summary(original))
    expect_equal(capture.output(print(stripped)), capture.output(print(original)))
    expect_equal(extract_network_matrix(stripped), extract_network_matrix(original))
    for (what in c("sd", "mean_u", "draws_u"))
      expect_equal(extract_random_effects(stripped, what),
                   extract_random_effects(original, what))
    expect_equal(bf(stripped), bf(original))
    expect_equal(bf(stripped, variable = "y_1"), bf(original, variable = "y_1"))
  }
})

test_that("real fits keep Stan data by default and omit it when requested", {
  bf <- function(...) suppressWarnings(bf_table(...))
  skip_if_not_installed("cmdstanr")
  for (family in c("gaussian", "bernoulli", "ordinal"))
    skip_if_no_bvarnet_models(family)

  d <- expand.grid(beep = 1:8, day = 1:2, id = c("private-a", "private-b"))
  d$g <- sin(seq_len(nrow(d)))
  d$b <- seq_len(nrow(d)) %% 2
  d$o <- seq_len(nrow(d)) %% 3 + 1
  d$x <- seq_len(nrow(d)) / 20
  d$unused_sensitive <- "private-unused-column"
  for (family in list("gaussian", "bernoulli", "ordinal",
                      c("gaussian", "bernoulli", "ordinal"))) {
    cols <- if (length(family) > 1) c("g", "b", "o") else
      switch(family, gaussian = "g", bernoulli = "b", ordinal = "o")
    args <- list(id_col = "id", time_col = "day", beep_col = "beep",
                 y_cols = cols, x_cols = "x", data = d, family = family,
                 K = 2, center_x = TRUE, re_cols = "Intercept",
                 fe_interactions = list(c("lag", "x")),
                 chains = 2, cores = 2, iter = 200, warmup = 100, seed = 1711,
                 init = 0, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    kept <- suppressWarnings(suppressMessages(do.call(bvar, args)))
    removed <- suppressWarnings(suppressMessages(do.call(bvar, c(args, list(save_data = FALSE)))))
    removed <- unserialize(serialize(removed, NULL))

    expect_type(kept$standata, "list")
    expect_true(all(c("Y", "X", "B", "Z", "id") %in% names(kept$standata)))
    expect_null(removed$standata)
    expect_false("data_used" %in% names(kept))
    expect_false("data_used" %in% names(kept$standata))
    expect_false("data_used" %in% names(removed))
    expect_false(grepl("private-", rawToChar(serialize(removed, NULL, ascii = TRUE)), fixed = TRUE))
    expect_equal(removed$model_info, kept$model_info)
    expect_identical(removed$draws, kept$draws)
    expect_true(all(removed$return_codes == 0))
    expect_equal(extract_param(removed), extract_param(kept))
    expect_equal(summary(removed), summary(kept))
    expect_equal(bf(removed, type = "lag_fe"), bf(kept, type = "lag_fe"))
    expect_equal(extract_random_effects(removed, "draws_u"),
                 extract_random_effects(kept, "draws_u"))
  }
})
