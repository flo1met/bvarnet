bvar <- function(id_col,
                 time_col,
                 y_cols,
                 x_cols,
                 center_x,
                 interactions = NULL,
                 fe_interactions = NULL,
                 re_interactions = NULL,
                 re_cols = NULL,
                 re_temporal = FALSE,
                 K,
                 data,
                 family = "bernoulli",
                 priors,
                 iter = 4000,
                 warmup = 1000,
                 chains = 4,
                 cores = 1,
                 seed = NULL

  ) {
  stanmodel <- instantiate::stan_package_model(name = "model_binary", package = "bvarnet")

  standata <- to_stan_data(data = data,
                            id_col = id_col,
                            time_col = time_col,
                            y_cols = y_cols,
                            x_cols = x_cols,
                            center_x = center_x,
                            fe_interactions = fe_interactions,
                            re_interactions = re_interactions,
                            re_cols = re_cols,
                            re_temporal = re_temporal,
                            K = K)

  stanfit <- stanmodel$sample(data = standata,
                              seed = seed,
                              iter_warmup = warmup,
                              iter_sampling = iter,
                              chains = chains,
                              parallel_chains = cores)

  return(stanfit)
}
