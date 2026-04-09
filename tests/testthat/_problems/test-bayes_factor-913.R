# Extracted from test-bayes_factor.R:913

# test -------------------------------------------------------------------------
mock <- make_mock_bvarnet("bernoulli")
expect_error(
    bf_table(mock, variable = "y_1", type = "fe"),
    "covariate effects"
  )
