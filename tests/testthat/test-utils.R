# ──────────────────────────────────────────────────────────────────────────────
# test-utils.R — internal utilities.
# ──────────────────────────────────────────────────────────────────────────────

test_that("%||% returns the fallback only for NULL", {
  expect_equal(NULL %||% 2, 2)
  expect_equal(1 %||% 2, 1)
  expect_equal(NA %||% 2, NA)
  expect_equal(list() %||% 2, list())
})

test_that("%||% is defined by the package, not inherited from base R >= 4.4", {
  # DESCRIPTION declares R (>= 3.5), where base has no `%||%`. If the package
  # stops defining it, this passes on a new R and fails for the users we claim
  # to support — so assert on the namespace, not on the call working.
  expect_true(exists("%||%", envir = asNamespace("bvarnet"), inherits = FALSE))
})

test_that("no package code calls an undefined global function", {
  skip_on_cran()
  skip_if_not_installed("codetools")

  reports <- character(0)
  codetools::checkUsagePackage(
    "bvarnet",
    report = function(msg) reports <<- c(reports, msg)
  )
  undefined <- grep("no visible global function definition", reports, value = TRUE)
  expect_equal(undefined, character(0))
})
