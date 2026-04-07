# test-parse_family.R — unit tests for .parse_family() and family helpers
# ──────────────────────────────────────────────────────────────────────────────

test_that(".parse_family() recycles scalar to named vector", {
  y_cols <- c("y_1", "y_2", "y_3")
  result <- .parse_family("gaussian", y_cols)
  expect_length(result, 3L)
  expect_equal(names(result), y_cols)
  expect_true(all(result == "gaussian"))
})

test_that(".parse_family() accepts named vector (reordered to y_cols)", {
  y_cols <- c("y_1", "y_2")
  fam <- c(y_2 = "ordinal", y_1 = "bernoulli")
  result <- .parse_family(fam, y_cols)
  expect_equal(as.character(result), c("bernoulli", "ordinal"))
  expect_equal(names(result), y_cols)
})

test_that(".parse_family() accepts unnamed vector of length p", {
  y_cols <- c("y_1", "y_2")
  result <- .parse_family(c("gaussian", "bernoulli"), y_cols)
  expect_equal(as.character(result), c("gaussian", "bernoulli"))
  expect_equal(names(result), y_cols)
})

test_that(".parse_family() errors on invalid family string", {
  expect_error(.parse_family("poisson", c("y_1", "y_2")))
})

test_that(".parse_family() errors on wrong length", {
  expect_error(.parse_family(c("gaussian", "bernoulli", "ordinal"), c("y_1", "y_2")))
})

test_that(".parse_family() errors on mismatched names", {
  fam <- c(y_1 = "gaussian", y_3 = "bernoulli")
  expect_error(.parse_family(fam, c("y_1", "y_2")))
})

# ── family helpers ──────────────────────────────────────────────────────────

test_that(".family_has() detects families in vector", {
  obj <- list(family = c(y_1 = "gaussian", y_2 = "ordinal"))
  expect_true(.family_has(obj, "gaussian"))
  expect_true(.family_has(obj, "ordinal"))
  expect_false(.family_has(obj, "bernoulli"))
})

test_that(".family_which() returns correct indices", {
  obj <- list(family = c(y_1 = "gaussian", y_2 = "ordinal", y_3 = "gaussian"))
  expect_equal(unname(.family_which(obj, "gaussian")), c(1L, 3L))
  expect_equal(unname(.family_which(obj, "ordinal")), 2L)
})

test_that(".is_mixed() distinguishes mixed from homogeneous", {
  mixed <- list(family = c(y_1 = "gaussian", y_2 = "bernoulli"))
  homo  <- list(family = c(y_1 = "gaussian", y_2 = "gaussian"))
  expect_true(.is_mixed(mixed))
  expect_false(.is_mixed(homo))
})

test_that(".format_family() formats correctly", {
  fam_scalar <- c(y_1 = "gaussian", y_2 = "gaussian")
  expect_equal(.format_family(fam_scalar), "gaussian")

  fam_mixed <- c(y_1 = "gaussian", y_2 = "ordinal")
  result <- .format_family(fam_mixed)
  expect_true(grepl("gaussian", result))
  expect_true(grepl("ordinal", result))
})
