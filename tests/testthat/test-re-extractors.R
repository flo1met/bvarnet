# ──────────────────────────────────────────────────────────────────────────────
# test-re-extractors.R — tests for .extract_u_draws() and .posterior_mean_u()
# ──────────────────────────────────────────────────────────────────────────────

# ═══════════════════════════════════════════════════════════════════════════════
# §1 .extract_u_draws()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".extract_u_draws returns correct 4D array dimensions", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 2L, J = 5L,
                            n_iter = 20L, n_chains = 2L)
  u <- .extract_u_draws(obj)

  expect_true(is.array(u))
  expect_equal(length(dim(u)), 4L)
  # [draw, node, subject, re]
  S <- 20L * 2L
  expect_equal(dim(u), c(S, 2L, 5L, 2L))
})

test_that(".extract_u_draws has informative dimnames", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 1L, J = 3L)
  u <- .extract_u_draws(obj)

  dn <- dimnames(u)
  expect_named(dn, c("draw", "node", "subject", "re"))
  expect_equal(dn$node, c("y_1", "y_2"))
  expect_equal(dn$subject, as.character(1:3))
  # RE names come from Z colnames
  expect_equal(dn$re, "z_1")
})

test_that(".extract_u_draws errors when n_re = 0", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 0L)
  expect_error(.extract_u_draws(obj), "n_re = 0")
})

test_that(".extract_u_draws recovers deterministic values", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 1L, J = 2L,
                            n_iter = 4L, n_chains = 1L)
  u <- .extract_u_draws(obj)

  # Manually extract u[1,1,1] from draws to verify
  raw <- obj$draws[, , "u[1,1,1]", drop = TRUE]
  flat <- as.vector(raw)
  expect_equal(u[, 1, 1, 1], flat)

  # u[2,2,1]
  raw2 <- obj$draws[, , "u[2,2,1]", drop = TRUE]
  flat2 <- as.vector(raw2)
  expect_equal(u[, 2, 2, 1], flat2)
})

test_that(".extract_u_draws works for all families", {
  for (fam in c("bernoulli", "gaussian", "ordinal")) {
    obj <- make_mock_bvarnet(family = fam, n_re = 1L, J = 3L)
    u <- .extract_u_draws(obj)
    S <- 20L * 2L
    expect_equal(dim(u), c(S, 2L, 3L, 1L),
                 info = paste("family:", fam))
  }
})

test_that(".extract_u_draws works with multiple REs", {
  obj <- make_mock_bvarnet(family = "gaussian", n_re = 3L, J = 4L)
  u <- .extract_u_draws(obj)
  S <- 20L * 2L
  expect_equal(dim(u), c(S, 2L, 4L, 3L))
  expect_equal(dimnames(u)$re, c("z_1", "z_2", "z_3"))
})

# ═══════════════════════════════════════════════════════════════════════════════
# §2 .posterior_mean_u()
# ═══════════════════════════════════════════════════════════════════════════════

test_that(".posterior_mean_u returns correct 3D array", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 2L, J = 5L,
                            n_iter = 20L, n_chains = 2L)
  u_mean <- .posterior_mean_u(obj)

  expect_true(is.array(u_mean))
  expect_equal(length(dim(u_mean)), 3L)
  # [node, subject, re]
  expect_equal(dim(u_mean), c(2L, 5L, 2L))
})

test_that(".posterior_mean_u equals manual column mean", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 1L, J = 2L,
                            n_iter = 10L, n_chains = 1L)
  u_mean <- .posterior_mean_u(obj)

  # Manually compute for u[1,1,1]
  raw <- as.vector(obj$draws[, , "u[1,1,1]", drop = TRUE])
  expect_equal(u_mean[1, 1, 1], mean(raw))

  # u[2,2,1]
  raw2 <- as.vector(obj$draws[, , "u[2,2,1]", drop = TRUE])
  expect_equal(u_mean[2, 2, 1], mean(raw2))
})

test_that(".posterior_mean_u errors when n_re = 0", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 0L)
  expect_error(.posterior_mean_u(obj), "n_re = 0")
})

test_that(".posterior_mean_u has informative dimnames", {
  obj <- make_mock_bvarnet(family = "gaussian", n_re = 2L, J = 3L)
  u_mean <- .posterior_mean_u(obj)

  dn <- dimnames(u_mean)
  expect_equal(dn[[1]], c("y_1", "y_2"))
  expect_equal(dn[[2]], as.character(1:3))
  expect_equal(dn[[3]], c("z_1", "z_2"))
})

# ═══════════════════════════════════════════════════════════════════════════════
# §3 Backward compatibility: make_mock_bvarnet with n_re = 0
# ═══════════════════════════════════════════════════════════════════════════════

test_that("make_mock_bvarnet with n_re=0 still works for existing tests", {
  obj <- make_mock_bvarnet(family = "bernoulli")
  expect_s3_class(obj, "bvarnet")
  expect_equal(obj$standata$n_re, 0L)
  # No u draws
  u_idx <- grep("^u\\[", dimnames(obj$draws)[[3]])
  expect_length(u_idx, 0L)
})

test_that("make_mock_bvarnet with n_re>0 has sd_u draws", {
  obj <- make_mock_bvarnet(family = "bernoulli", n_re = 2L, J = 3L)
  sd_u_idx <- grep("^sd_u\\[", dimnames(obj$draws)[[3]])
  # p=2 nodes, n_re=2 → 4 sd_u parameters

  expect_length(sd_u_idx, 4L)
  # All sd_u draws should be positive
  expect_true(all(obj$draws[, , sd_u_idx] > 0))
})
