#' Fit a Bayesian multilevel VAR network model
#'
#' Compiles and samples the appropriate Stan model for the chosen family,
#' extracts all results into plain base-R objects, and returns a \code{bvarnet}
#' object.
#'
#' @param id_col Character. Name of the subject/group identifier column.
#' @param time_col Character. Name of the time column.
#' @param y_cols Character vector. Names of the outcome columns.
#' @param x_cols Character vector. Names of the covariate columns.
#' @param center_x Logical. Grand-mean centre covariates before fitting?
#'   Default \code{FALSE}.
#' @param fe_interactions List or NULL. Fixed-effect interaction terms to add
#'   to the design matrix. Each element is a character vector of column names
#'   to interact, or \code{c("lag", "x")} to interact all lag columns with
#'   a covariate.
#' @param re_interactions List or NULL. Random-effect interaction terms.
#' @param re_cols Character vector. Columns from X to include as random slopes.
#' @param re_temporal Logical. Include random slopes on lag predictors?
#'   Default \code{FALSE}.
#' @param K Integer. AR order. Default 1.
#' @param na_action Character. Missing-data strategy; currently only
#'   \code{"listwise"}.
#' @param skip_lag Logical. If \code{TRUE} (default), rows with irregular time
#'   gaps have their lag set to zero rather than being dropped.
#' @param data Data frame in long format.
#' @param family Character. Observation model: \code{"bernoulli"},
#'   \code{"ordinal"}, or \code{"gaussian"}.
#' @param priors A \code{bvarnet_priors} object from \code{set_priors()}.
#'   Defaults to \code{set_priors()} (package defaults).
#' @param iter Integer. Number of post-warmup iterations per chain. Default 4000.
#' @param warmup Integer. Number of warmup iterations per chain. Default 1000.
#' @param chains Integer. Number of MCMC chains. Default 4.
#' @param cores Integer. Number of chains to run in parallel. Default 1.
#' @param seed Integer or NULL. RNG seed.
#' @param keep_fit Logical. If \code{TRUE}, the raw \code{CmdStanMCMC} object
#'   is retained in the returned object as \code{$fit}. Useful for accessing
#'   Stan profiling results (\code{$fit$profiles()}), LOO/WAIC helpers, or
#'   other CmdStanR methods. Default \code{FALSE} to keep memory usage low.
#'
#' @return A \code{bvarnet} object (a named list) with slots:
#'   \code{draws}, \code{summary}, \code{diagnostics}, \code{timing},
#'   \code{metadata}, \code{return_codes}, \code{family}, \code{standata},
#'   \code{priors}. If \code{keep_fit = TRUE}, also \code{fit} (the raw
#'   \code{CmdStanMCMC} object).
#'
#' @export
bvar <- function(id_col,
                 time_col,
                 y_cols,
                 x_cols,
                 center_x = FALSE,
                 fe_interactions = NULL,
                 re_interactions = NULL,
                 re_cols = NULL,
                 re_temporal = FALSE,
                 K = 1,
                 na_action = c("listwise"),
                 skip_lag = TRUE,
                 data,
                 family = c("bernoulli", "ordinal", "gaussian"),
                 priors = set_priors(),
                 iter = 4000,
                 warmup = 1000,
                 chains = 4,
                 cores = 1,
                 seed = NULL,
                 keep_fit = FALSE

  ) {

  family <- match.arg(family)
  model_name <- switch(family,
                       bernoulli = "model_binary",
                       ordinal   = "model_ordinal",
                       gaussian  = "model_gaussian",
                       stop("Unknown family: ", family)
  )
  stanmodel <- instantiate::stan_package_model(name = model_name, package = "bvarnet")

  standata <- to_stan_data(data = data,
                            family = family,
                            id_col = id_col,
                            time_col = time_col,
                            y_cols = y_cols,
                            x_cols = x_cols,
                            center_x = center_x,
                            fe_interactions = fe_interactions,
                            re_interactions = re_interactions,
                            re_cols = re_cols,
                            re_temporal = re_temporal,
                            K = K,
                            na_action = na_action,
                            skip_lag = skip_lag,
                            priors = priors
                           )

  stanfit <- stanmodel$sample(data = standata,
                              seed = seed,
                              iter_warmup = warmup,
                              iter_sampling = iter,
                              chains = chains,
                              parallel_chains = cores)

  # Extract everything from CmdStanMCMC into plain base-R objects, then discard
  # the fit object (CSV refs, compiled binary, lazy draws) to keep memory lean.
  raw_draws   <- stanfit$draws(format = "array")
  draws       <- unclass(raw_draws)         # strip draws_array class; dimnames preserved
  attr(draws, "class") <- NULL              # ensure it is a plain array
  summary     <- as.data.frame(stanfit$summary())
  diagnostics <- as.data.frame(stanfit$diagnostic_summary())
  timing      <- stanfit$time()
  metadata    <- stanfit$metadata()
  return_codes <- stanfit$return_codes()

  out <- structure(
    list(
      draws        = draws,
      summary      = summary,
      diagnostics  = diagnostics,
      timing       = timing,
      metadata     = metadata,
      return_codes = return_codes,
      family       = family,
      standata     = standata,
      priors       = priors
    ),
    class = "bvarnet"
  )
  out
}
