# ──────────────────────────────────────────────────────────────────────────────
# test-draw-ordering.R — consumers of the draws array must resolve parameters by
# parsed index, never by column position.
#
# CmdStan flattens every index column-major (first index fastest), while
# .combine_nodewise_fits() binds node-outer. The two orders disagree for sd_u,
# so any positional reshape is wrong on at least one path. These tests poke
# values into named cells and check where they land.
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 .stan_indices()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".stan_indices parses one row per name", {
  expect_equal(.stan_indices(c("u[1,2,3]", "u[10,1,2]")),
               matrix(c(1L, 2L, 3L, 10L, 1L, 2L), nrow = 2, byrow = TRUE))
  expect_equal(.stan_indices("sigma[2]"), matrix(2L, nrow = 1))
})

test_that(".stan_indices rejects unparseable or ragged names", {
  expect_error(.stan_indices("lp__"), "Cannot parse")
  expect_error(.stan_indices(c("u[1,1]", "u[1,1,1]")), "inconsistent index")
})


# ═══════════════════════════════════════════════════════════════════════════════
# §2 The fixture reproduces true CmdStan order
# ═══════════════════════════════════════════════════════════════════════════════

test_that("mock draws are emitted first-index-fastest, like CmdStan", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 2L, J = 3L)
  nms <- dimnames(obj$draws)[[3]]

  # u: array[p] matrix[J, n_re] -> node fastest, then subject, then re
  expect_equal(head(grep("^u\\[", nms, value = TRUE), 4L),
               c("u[1,1,1]", "u[2,1,1]", "u[1,2,1]", "u[2,2,1]"))

  # sd_u: matrix[p, n_re] -> node (row) fastest
  expect_equal(grep("^sd_u\\[", nms, value = TRUE),
               c("sd_u[1,1]", "sd_u[2,1]", "sd_u[1,2]", "sd_u[2,2]"))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §3 .extract_u_draws() maps by name
# ═══════════════════════════════════════════════════════════════════════════════

test_that("extract_random_effects maps u draws by name, not position", {
  obj <- make_mock_bvarnet(family = "gaussian", n_re = 2L, J = 3L,
                            n_iter = 4L, n_chains = 1L)

  # Cells chosen so the positional reshape sends them somewhere else entirely.
  obj$draws[1, 1, "u[2,3,1]"] <- 99
  obj$draws[1, 1, "u[1,1,2]"] <- 77

  arr <- extract_random_effects(obj, what = "draws_u")   # [draw, node, subj, re]
  expect_equal(arr[1, "y_2", "3", 1], 99)
  expect_equal(arr[1, "y_1", "1", 2], 77)
})

test_that(".extract_u_draws round-trips every named cell", {
  obj <- make_mock_bvarnet(family = "gaussian", n_re = 2L, J = 3L,
                            n_iter = 4L, n_chains = 1L)
  arr <- .extract_u_draws(obj)

  for (nm in grep("^u\\[", dimnames(obj$draws)[[3]], value = TRUE)) {
    ix <- .stan_indices(nm)
    expect_equal(as.vector(arr[, ix[1], ix[2], ix[3]]),
                 as.vector(obj$draws[, , nm]),
                 info = nm)
  }
})

test_that(".extract_u_draws errors when a u cell is missing from the draws", {
  obj <- make_mock_bvarnet(family = "gaussian", n_re = 2L, J = 3L)
  keep <- dimnames(obj$draws)[[3]] != "u[2,3,1]"
  obj$draws <- obj$draws[, , keep, drop = FALSE]
  expect_error(.extract_u_draws(obj), "Expected .* u parameters")
})

test_that(".extract_u_draws handles degenerate dimensions", {
  # Single RE, single subject: every flattening order coincides here, which is
  # exactly why positional bugs survive such cases.
  obj <- make_mock_bvarnet(family = "gaussian", n_re = 1L, J = 1L,
                            n_iter = 4L, n_chains = 1L)
  arr <- .extract_u_draws(obj)
  expect_equal(dim(arr), c(4L, 2L, 1L, 1L))
  expect_equal(as.vector(arr[, 2, 1, 1]),
               as.vector(obj$draws[, , "u[2,1,1]"]))
})


# ═══════════════════════════════════════════════════════════════════════════════
# §4 sd_u labels survive the nodewise column order
# ═══════════════════════════════════════════════════════════════════════════════

# Rewrite the sd_u block into the order .combine_nodewise_fits() produces:
# each node's row is emitted contiguously, so re varies fastest.
as_nodewise_sd_u_order <- function(obj) {
  nms    <- dimnames(obj$draws)[[3]]
  sd_idx <- grep("^sd_u\\[", nms)
  ix     <- .stan_indices(nms[sd_idx])
  perm   <- order(ix[, 1L], ix[, 2L])        # node slowest, re fastest

  obj$draws[, , sd_idx] <- obj$draws[, , sd_idx[perm], drop = FALSE]
  dimnames(obj$draws)[[3]][sd_idx] <- nms[sd_idx][perm]
  obj
}

test_that("Random Effect SD rows are labelled from names, not column position", {
  for (obj in list(make_mock_bvarnet(family = "gaussian", n_re = 2L, J = 3L),
                   as_nodewise_sd_u_order(
                     make_mock_bvarnet(family = c("gaussian", "ordinal"),
                                       n_re = 2L, J = 3L)))) {
    tab <- extract_param(obj, type = "Random Effect SD")

    d  <- extract_draws(obj, "sd_u")
    ix <- .stan_indices(colnames(d))
    nm <- get_param_names(obj$standata)

    expect_equal(tab$predictor, nm$y[ix[, 1L]])
    expect_equal(tab$outcome,   nm$re[ix[, 2L]])
    # The statistics must belong to the row they are labelled with.
    expect_equal(tab$mean, unname(colMeans(d)))
  }
})

test_that("extract_random_effects(what='sd') agrees with extract_param", {
  obj <- as_nodewise_sd_u_order(
    make_mock_bvarnet(family = "gaussian", n_re = 2L, J = 3L))

  from_extractor <- extract_random_effects(obj, what = "sd")
  from_param     <- extract_param(obj, type = "Random Effect SD")
  rownames(from_param) <- NULL

  expect_equal(from_extractor$predictor, from_param$predictor)
  expect_equal(from_extractor$outcome,   from_param$outcome)
  expect_equal(from_extractor$mean,      from_param$mean)
})


# ═══════════════════════════════════════════════════════════════════════════════
# §5 A real fit really does emit first-index-fastest
# ═══════════════════════════════════════════════════════════════════════════════

test_that("CmdStan emits u and sd_u with the first index varying fastest", {
  skip_if_no_bvarnet_models()

  sim <- sim_var(N = 3, T_obs = 12, p = 2, K = 1, family = "gaussian", seed = 3)
  fit <- bvar(id_col = "id", time_col = "t", y_cols = c("y_1", "y_2"),
              K = 1, family = "gaussian", data = sim$data, re_cols = c("Intercept"),
              re_temporal = TRUE, iter = 200, warmup = 150,
              chains = 1, seed = 3)
  skip_if(fit$standata$n_re < 2, "need n_re > 1 to discriminate orderings")

  nms <- dimnames(fit$draws)[[3]]
  expect_equal(head(grep("^u\\[", nms, value = TRUE), 2L),
               c("u[1,1,1]", "u[2,1,1]"))
  expect_equal(head(grep("^sd_u\\[", nms, value = TRUE), 2L),
               c("sd_u[1,1]", "sd_u[2,1]"))

  # And the extractor agrees with the named cells.
  arr <- extract_random_effects(fit, what = "draws_u")
  expect_equal(as.vector(arr[, 2, 3, 1]), as.vector(fit$draws[, , "u[2,3,1]"]))
})
