# Extracted from test-bvarnet-object.R:177

# test -------------------------------------------------------------------------
obj  <- make_mock_bvarnet("bernoulli")
res  <- bvarnet:::extract_draws(obj, "beta")
