# Extracted from test-bvarnet-object.R:150

# test -------------------------------------------------------------------------
obj    <- make_mock_bvarnet("gaussian")
result <- bvarnet:::extract_draws(obj, "beta")
