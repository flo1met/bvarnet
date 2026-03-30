#' Extract temporal (VAR lag) effects
#'
#' Returns a data frame of autoregressive and/or cross-lagged parameter
#' summaries with convergence diagnostics, filtered by lag and effect type.
#'
#' @param object A \code{bvarnet} object returned by \code{\link{bvar}}.
#' @param lag Integer or \code{NULL}. If specified, only effects from this
#'   lag are returned. Default \code{NULL} (all lags).
#' @param effect Character. One of \code{"all"} (default), \code{"ar"}
#'   (autoregressive only), or \code{"cl"} (cross-lagged only).
#' @param bayes_factor Logical; if \code{TRUE}, append BF columns.
#'   Default \code{FALSE}.
#' @param null_value Numeric; null hypothesis for BF. Default 0.
#'
#' @return A data frame with columns \code{type}, \code{predictor},
#'   \code{outcome}, \code{mean}, \code{median}, \code{q5}, \code{q95},
#'   \code{rhat}, \code{ess_bulk}, \code{ess_tail}, and optionally
#'   \code{BF01}, \code{BF10}.
#'
#' @export
extract_temporal <- function(object,
                             lag = NULL,
                             effect = c("all", "ar", "cl"),
                             bayes_factor = FALSE,
                             null_value = 0) {
  stopifnot(inherits(object, "bvarnet"))
  effect <- match.arg(effect)

  tab <- extract_param(object,
                        bayes_factor = bayes_factor,
                        null_value   = null_value)

  # Filter to temporal types
  type_filter <- switch(effect,
    ar  = "Autoregressive",
    cl  = "Cross-lagged",
    all = c("Autoregressive", "Cross-lagged")
  )
  out <- tab[tab$type %in% type_filter, , drop = FALSE]

  # Optionally filter by lag
  if (!is.null(lag)) {
    stopifnot(is.numeric(lag), length(lag) == 1L, lag >= 1L)
    lag_pattern <- paste0("^lag", as.integer(lag), "_")
    out <- out[grepl(lag_pattern, out$predictor), , drop = FALSE]
  }

  rownames(out) <- NULL
  out
}


#' Extract random-effect summaries
#'
#' Returns random-effect standard deviations (group-level variance),
#' subject-level posterior means, or the full posterior draws of the
#' subject-level random effects \code{u}.
#'
#' @param object A \code{bvarnet} object returned by \code{\link{bvar}}.
#' @param what Character. What to extract:
#'   \describe{
#'     \item{\code{"sd"}}{Data frame of random-effect SD summaries
#'       (from \code{extract_param}).}
#'     \item{\code{"mean_u"}}{3D array \code{[node, subject, re]} of
#'       posterior means of subject-level effects.}
#'     \item{\code{"draws_u"}}{4D array \code{[draw, node, subject, re]}
#'       of full posterior draws.}
#'   }
#'
#' @return Depends on \code{what}; see above.
#'
#' @export
extract_random_effects <- function(object,
                                   what = c("sd", "mean_u", "draws_u")) {
  stopifnot(inherits(object, "bvarnet"))
  what <- match.arg(what)

  if (object$standata$n_re == 0L)
    stop("No random effects in this model (n_re = 0).", call. = FALSE)

  switch(what,
    sd = {
      tab <- extract_param(object)
      out <- tab[tab$type == "Random Effect SD", , drop = FALSE]
      rownames(out) <- NULL
      out
    },
    mean_u = .posterior_mean_u(object),
    draws_u = .extract_u_draws(object)
  )
}


#' Extract a network matrix of temporal coefficients
#'
#' Returns a named \code{p x p} matrix of posterior summary statistics for
#' the VAR lag coefficients at a chosen lag, suitable for network
#' visualisation (e.g., with \pkg{igraph} or \pkg{qgraph}).
#'
#' @param object A \code{bvarnet} object returned by \code{\link{bvar}}.
#' @param lag Integer. Which lag block. Default 1.
#' @param stat Character. Summary statistic to fill the matrix with:
#'   \code{"mean"} (default), \code{"median"}, \code{"q5"}, or \code{"q95"}.
#'
#' @return A named \code{p x p} numeric matrix. Element \code{[i, j]}
#'   gives the effect of variable \code{i} (lagged) on variable \code{j}
#'   (outcome). Row and column names are the outcome variable names.
#'
#' @export
extract_network_matrix <- function(object,
                                   lag = 1L,
                                   stat = c("mean", "median", "q5", "q95")) {
  stopifnot(inherits(object, "bvarnet"))
  stat <- match.arg(stat)

  sd <- object$standata
  p  <- sd$p
  K  <- sd$K
  stopifnot(lag >= 1L, lag <= K)

  nm <- get_param_names(sd)

  # Extract phi draws and build summary
  draws_phi <- extract_draws(object, "phi")
  phi_tab   <- build_summary_table(draws_phi, nm$b, nm$y, "phi")

  # Filter to the requested lag block
  lag_pattern <- paste0("^lag", as.integer(lag), "_")
  lag_rows <- grepl(lag_pattern, phi_tab$predictor)
  sub <- phi_tab[lag_rows, , drop = FALSE]

  # Build p x p matrix: sub rows cycle through predictors (lagged vars)
  # for each outcome column.
  mat <- matrix(sub[[stat]], nrow = p, ncol = p)
  y_names <- nm$y
  rownames(mat) <- y_names
  colnames(mat) <- y_names
  mat
}
