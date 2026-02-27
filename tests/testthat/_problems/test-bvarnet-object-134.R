# Extracted from test-bvarnet-object.R:134

# test -------------------------------------------------------------------------
n_iter <- 20L
n_chains <- 2L
obj    <- make_mock_bvarnet("bernoulli", n_iter = n_iter, n_chains = n_chains)
result <- bvarnet:::extract_draws(obj, "beta")
