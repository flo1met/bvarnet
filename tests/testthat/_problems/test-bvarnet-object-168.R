# Extracted from test-bvarnet-object.R:168

# test -------------------------------------------------------------------------
obj    <- make_mock_bvarnet("ordinal")
result <- bvarnet:::extract_draws(obj, "kappa")
