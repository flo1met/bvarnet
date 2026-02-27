# Extracted from test-bvarnet-object.R:124

# test -------------------------------------------------------------------------
obj    <- make_mock_bvarnet("bernoulli")
result <- bvarnet:::extract_draws(obj, "beta")
