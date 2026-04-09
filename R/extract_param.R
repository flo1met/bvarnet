#' Extract labelled parameter summaries from a fitted bvarnet model
#'
#' Returns a single flat data frame with posterior summaries (mean, median,
#' 5th/95th percentiles) and convergence diagnostics (Rhat, ESS) for all
#' model parameters.
#'
#' @param object A \code{bvarnet} object returned by \code{bvar()}.
#' @param bayes_factor Logical; if \code{TRUE}, append \code{BF01} and
#'   \code{BF10} columns computed via the Savage-Dickey density ratio for
#'   beta and phi parameters.  Default \code{FALSE}.
#' @param null_value Numeric scalar; the null hypothesis value for Bayes
#'   factor computation (default 0).  Only used when \code{bayes_factor = TRUE}.
#'
#' @return A data frame with columns: \code{type}, \code{predictor},
#'   \code{outcome}, \code{mean}, \code{median}, \code{q5}, \code{q95},
#'   \code{rhat}, \code{ess_bulk}, \code{ess_tail}, and optionally
#'   \code{BF01}, \code{BF10}.
#'
#' @export
extract_param <- function(object, bayes_factor = FALSE, null_value = 0) {
  stopifnot(inherits(object, "bvarnet"))

  sd   <- object$standata
  nm   <- get_param_names(sd)
  smry <- object$convergence   # data.frame: variable, rhat, ess_bulk, ess_tail

  # Join Rhat + ESS from object$convergence by Stan parameter name.
  # stan_colnames: character vector aligned with rows of tab.
  join_convergence <- function(tab, stan_colnames) {
    idx          <- match(stan_colnames, smry$variable)
    tab$rhat     <- smry$rhat[idx]
    tab$ess_bulk <- smry$ess_bulk[idx]
    tab$ess_tail <- smry$ess_tail[idx]
    tab
  }

  # ---------- Intercepts & fixed effects (beta) ----------
  draws_beta <- extract_draws(object, "beta")

  # Remove ordinal-intercept NA sentinel columns (D4) BEFORE summary.
  # Sentinels only exist in mixed-family models where X retains the Intercept
  # row; pure ordinal strips the Intercept so beta[1,j] is a real covariate.
  beta_cols <- colnames(draws_beta)
  has_intercept <- "Intercept" %in% colnames(sd$X)
  if (has_intercept) {
    ord_indices <- which(object$family == "ordinal")
    sentinel_cols <- paste0("beta[1,", ord_indices, "]")
    keep_beta_cols <- !beta_cols %in% sentinel_cols
  } else {
    keep_beta_cols <- rep(TRUE, length(beta_cols))
  }
  draws_beta_clean <- draws_beta[, keep_beta_cols, drop = FALSE]

  # Build filtered name lists for row/col mapping
  # Original order: (fe1,y1), (fe2,y1), ..., (fe_nfe,y1), (fe1,y2), ...
  all_beta_indices <- expand.grid(fe = seq_len(length(nm$fe)),
                                   node = seq_len(length(nm$y)))
  keep_idx <- keep_beta_cols
  kept_fe <- all_beta_indices$fe[keep_idx]
  kept_node <- all_beta_indices$node[keep_idx]

  # Build summary table manually for the filtered columns
  d_mean   <- colMeans(draws_beta_clean)
  d_median <- apply(draws_beta_clean, 2L, stats::median)
  d_q5     <- apply(draws_beta_clean, 2L, stats::quantile, probs = 0.05)
  d_q95    <- apply(draws_beta_clean, 2L, stats::quantile, probs = 0.95)
  beta_tab <- data.frame(
    type      = ifelse(nm$fe[kept_fe] == "Intercept", "Intercept", "Fixed Effect"),
    predictor = nm$fe[kept_fe],
    outcome   = nm$y[kept_node],
    mean      = as.numeric(d_mean),
    median    = as.numeric(d_median),
    q5        = as.numeric(d_q5),
    q95       = as.numeric(d_q95),
    stringsAsFactors = FALSE
  )
  beta_stan_names <- colnames(draws_beta_clean)
  beta_tab <- join_convergence(beta_tab, beta_stan_names)

  # ---------- Autoregressive & Cross-lagged effects (phi) ----------
  draws_phi <- extract_draws(object, "phi")
  phi_tab   <- build_summary_table(draws_phi, nm$b, nm$y, "placeholder")

  # Classify each phi row as Autoregressive or Cross-lagged
  p_  <- sd$p
  K_  <- sd$K
  nr_ <- p_ * K_
  nc_ <- p_
  col_indices <- rep(seq_len(nc_), each = nr_)
  row_indices <- rep(seq_len(nr_), times = nc_)
  row_within  <- ((row_indices - 1L) %% p_) + 1L
  phi_tab$type <- ifelse(row_within == col_indices,
                         "Autoregressive", "Cross-lagged")
  phi_tab <- join_convergence(phi_tab, colnames(draws_phi))

  # ---------- Random-effect SDs (sd_u) ----------
  re_sd_tab <- if (sd$n_re > 0) {
    draws_sd <- extract_draws(object, "sd_u")
    tab      <- build_summary_table(draws_sd, nm$y, nm$re, "Random Effect SD")
    join_convergence(tab, colnames(draws_sd))
  } else NULL

  # ---------- Residual SD (sigma, gaussian only) ----------
  gauss_idx <- .family_which(object, "gaussian")
  sigma_tab <- if (length(gauss_idx) > 0) {
    draws_sigma <- extract_draws(object, "sigma")
    gauss_names <- nm$y[gauss_idx]
    tab <- build_summary_table(draws_sigma, gauss_names, "sigma", "Residual SD")
    join_convergence(tab, colnames(draws_sigma))
  } else NULL

  # ---------- Thresholds (kappa, ordinal only) ----------
  ord_idx <- .family_which(object, "ordinal")
  kappa_tab <- if (length(ord_idx) > 0) {
    draws_kappa <- extract_draws(object, "kappa")
    cn    <- colnames(draws_kappa)
    parts <- strsplit(gsub("kappa\\[|\\]", "", cn), ",")
    j_idx <- as.integer(vapply(parts, `[[`, character(1L), 1L))
    c_idx <- as.integer(vapply(parts, `[[`, character(1L), 2L))
    tab <- data.frame(
      type      = "Threshold",
      predictor = paste0("kappa(", nm$y[j_idx], ", c", c_idx, ")"),
      outcome   = "\u2014",
      mean      = colMeans(draws_kappa),
      median    = apply(draws_kappa, 2L, stats::median),
      q5        = apply(draws_kappa, 2L, stats::quantile, probs = 0.05),
      q95       = apply(draws_kappa, 2L, stats::quantile, probs = 0.95),
      stringsAsFactors = FALSE
    )
    join_convergence(tab, cn)
  } else NULL

  # ---------- Combine into a single flat data.frame ----------
  out <- do.call(rbind, Filter(Negate(is.null),
                               list(beta_tab, phi_tab, re_sd_tab, sigma_tab, kappa_tab)))
  rownames(out) <- NULL

  # ---------- Append BFs for beta and phi params ----------
  if (isTRUE(bayes_factor)) {
    out$BF01 <- NA_real_
    out$BF10 <- NA_real_

    # Identify rows with Stan param names we can compute BFs for
    # beta and phi rows have matching Stan colnames from their draws
    bf_types <- c("Intercept", "Fixed Effect", "Autoregressive", "Cross-lagged")
    bf_stan_names <- c(beta_stan_names, colnames(draws_phi))

    bf_idx <- which(out$type %in% bf_types)
    for (i in seq_along(bf_idx)) {
      row_i    <- bf_idx[i]
      stan_nm  <- bf_stan_names[i]
      res <- tryCatch(
        savage_dickey(object, params = stan_nm, null_value = null_value,
                      method = "logspline"),
        error = function(e) NULL
      )
      if (!is.null(res)) {
        out$BF01[row_i] <- res$BF01
        out$BF10[row_i] <- res$BF10
      }
    }
  }

  out
}
