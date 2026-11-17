location_init_fixture <- function(family = "gaussian") {
  ordinal <- family == "ordinal"
  list(p = 1L, K = 1L, n_fe = if (ordinal) 1L else 2L, C = 4L,
       X = if (ordinal) matrix(c(8, 12)) else cbind(1, c(8, 12)),
       B = matrix(if (family == "bernoulli") c(0, 1) else c(1, 3)))
}

test_that("complete public initializers map location and preserve other blocks", {
  sd <- location_init_fixture()
  init <- list(sigma = 1.1, phi = matrix(.4), beta = matrix(c(1.2, .3), 2),
               sd_u = matrix(.5), z_u = array(.2, c(1, 1, 1)))
  mapped <- .bvarnet_centered_init(init, sd, "gaussian", 1L)[[1]]
  expect_equal(mapped$beta_centered, matrix(c(5, .3), 2))
  expect_null(mapped[["beta", exact = TRUE]])
  for (nm in c("phi", "sigma", "sd_u", "z_u"))
    expect_identical(mapped[[nm]], init[[nm]])

  sd <- location_init_fixture("bernoulli")
  init$beta[1, 1] <- 1.8
  init$sigma <- NULL
  expect_equal(.bvarnet_centered_init(init, sd, "bernoulli", 1L)[[1]]$beta_centered[1, 1], 5)

  sd <- location_init_fixture("ordinal")
  init <- list(kappa = list(c(5.4, 4, 4.9)), phi = matrix(.4), beta = matrix(.3))
  mapped <- .bvarnet_centered_init(init, sd, "ordinal", 1L)[[1]]
  expect_equal(mapped$kappa_centered, matrix(c(1.6, .2, 1.1), 1))
  expect_identical(mapped$beta, init$beta)
  expect_null(mapped[["kappa", exact = TRUE]])
})

test_that("lists, functions and JSON initializers agree by name for every chain", {
  sd <- location_init_fixture()
  init <- list(phi = matrix(.4), beta = matrix(c(1.2, .3), 2))
  expected <- .bvarnet_centered_init(init, sd, "gaussian", 2L)
  expect_identical(.bvarnet_centered_init(list(init, init), sd, "gaussian", 2L), expected)
  expect_identical(.bvarnet_centered_init(function() init, sd, "gaussian", 2L), expected)
  expect_identical(.bvarnet_centered_init(function(chain_id) init, sd, "gaussian", 2L), expected)
  file <- tempfile(fileext = ".json")
  on.exit(unlink(file))
  jsonlite::write_json(init, file, auto_unbox = TRUE, digits = NA)
  expect_equal(.bvarnet_centered_init(file, sd, "gaussian", 2L), expected)
  expect_equal(.bvarnet_centered_init(rep(file, 2), sd, "gaussian", 2L), expected)
})

test_that("partial, malformed, mixed-reference and fit-based starts error clearly", {
  sd <- location_init_fixture()
  init <- list(beta = matrix(c(1.2, .3), 2), phi = matrix(.4))
  expect_error(.bvarnet_centered_init(list(beta = init$beta), sd, "gaussian", 1), "complete public")
  expect_error(.bvarnet_centered_init(list(phi = init$phi), sd, "gaussian", 1), "complete public")
  bad <- init; bad$beta <- t(bad$beta)
  expect_error(.bvarnet_centered_init(bad, sd, "gaussian", 1), "dimensions 2 x 1")
  bad <- init; bad$beta[1] <- NA_real_
  expect_error(.bvarnet_centered_init(bad, sd, "gaussian", 1), "finite matrix")
  bad <- init; bad$beta_centered <- bad$beta
  expect_error(.bvarnet_centered_init(bad, sd, "gaussian", 1), "internal centered")
  expect_error(.bvarnet_centered_init(list(init), sd, "gaussian", 2), "per chain")
  expect_error(.bvarnet_centered_init(new.env(), sd, "gaussian", 1), "Fit/draw-based")
  expect_error(.bvarnet_centered_init(function(unused) init, sd, "gaussian", 1), "chain_id")
  expect_error(.bvarnet_centered_init("missing-init.json", sd, "gaussian", 1), "existing JSON")
  expect_identical(.bvarnet_centered_init(0, sd, "gaussian", 1), 0)
  expect_identical(.bvarnet_centered_init(2, sd, "gaussian", 1), 2)
  expect_null(.bvarnet_centered_init(NULL, sd, "gaussian", 1))
})

test_that("ordinal zero-width beta and mixed node-specific functions are supported", {
  sd <- location_init_fixture("ordinal")
  sd$n_fe <- 0L; sd$X <- matrix(numeric(), 2, 0)
  init <- list(phi = matrix(.4), kappa = matrix(c(2.4, 1, 1.9), 1))
  mapped <- .bvarnet_centered_init(init, sd, "ordinal", 1)[[1]]
  expect_equal(mapped$kappa_centered, matrix(c(1.6, .2, 1.1), 1))
  expect_null(mapped[["beta", exact = TRUE]])
  f <- function(chain_id, node_id) {
    list(beta = matrix(c(node_id + chain_id, .3), 2), phi = matrix(.4))
  }
  mapped <- .bvarnet_centered_init(f, location_init_fixture(), "gaussian", 2, node_id = 3)
  expect_equal(vapply(mapped, function(x) x$beta_centered[1, 1], numeric(1)), c(7.8, 8.8))
})

test_that("partial matching cannot bypass initializer adaptation", {
  skip_if_not_installed("cmdstanr")
  expect_identical(.check_sampler_dots(init = 0), list(init = 0))
  expect_error(.check_sampler_dots(ini = 0), "Spell 'init' in full")
})
