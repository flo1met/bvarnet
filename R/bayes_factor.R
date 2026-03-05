# -----------------------------------------------------------------------------
# Savage-Dickey density ratio Bayes factor estimation
#
# Core functions:
#   - bf_table()          : primary export — tidy data frame of BFs
#   - savage_dickey()     : internal workhorse (named list per test)
#   - eval_prior_density(): analytical prior density evaluation
#   - get_phi_indices()   : resolve AR / CL / all phi parameter names
# -----------------------------------------------------------------------------


# ---- Prior density evaluation ------------------------------------------------

#' Evaluate the prior density at a point (internal)
#'
#' Computes the analytical density of a \code{bvarnet_prior} at \code{x}.
#'
#' @param prior A \code{bvarnet_prior} object.
#' @param x     Numeric scalar — the point at which to evaluate the density.
#'
#' @return Numeric scalar — the density.
#' @keywords internal
eval_prior_density <- function(prior, x) {
  # TODO: if stationarity constraints are added to Phi, the effective prior
	#       becomes truncated and this analytical density must be corrected.
  switch(prior$family,
    normal    = stats::dnorm(x,  mean = prior$loc, sd = prior$scale),
    student_t = stats::dt((x - prior$loc) / prior$scale, df = prior$df) / prior$scale,
    cauchy    = stats::dcauchy(x, location = prior$loc, scale = prior$scale),
    stop("Unknown prior family: ", prior$family, call. = FALSE)
  )
}


#' Evaluate the joint prior density for independent priors (internal)
#'
#' @param prior_list A named list of \code{bvarnet_prior} objects keyed by
#'   parameter type (e.g., \code{"phi"}, \code{"beta"}).
#' @param param_types Character vector of prior types for each parameter.
#' @param null_vec Numeric vector of null values (same length as
#'   \code{param_types}).
#'
#' @return Numeric scalar — the product of marginal densities.
#' @keywords internal
eval_joint_prior_density <- function(prior_list, param_types, null_vec) {
  prod(mapply(function(ptype, x) {
    eval_prior_density(prior_list[[ptype]], x)
  }, param_types, null_vec, USE.NAMES = FALSE))
}


# ---- Phi index resolution ----------------------------------------------------

#' Resolve Stan parameter names for phi sub-matrices (internal)
#'
#' Returns a character vector of Stan parameter names (e.g.,
#' \code{"phi[1,1]"}) for autoregressive, cross-lagged, or all effects.
#'
#' @param sd      The \code{standata} list from a \code{bvarnet} object.
#' @param lag     Integer; which lag block (1 to K). Default 1.
#' @param effect  One of \code{"ar"}, \code{"cl"}, \code{"all"}.
#'
#' @return Character vector of Stan parameter names.
#' @keywords internal
get_phi_indices <- function(sd, lag = 1L, effect = c("ar", "cl", "all")) {
  effect <- match.arg(effect)
  p <- sd$p
  K <- sd$K
  if (lag < 1L || lag > K)
    stop(sprintf("'lag' must be between 1 and K = %d.", K), call. = FALSE)

  # phi is matrix[p*K, p] — phi[row, col]
  # For lag block k: rows (k-1)*p + 1 .. k*p
  # AR: row = (k-1)*p + j, col = j  (self-loops)
  # CL: everything else in that lag block
  offset <- (lag - 1L) * p
  nms    <- character(0)

  for (col in seq_len(p)) {
    for (row_within in seq_len(p)) {
      row <- offset + row_within
      is_ar <- (row_within == col)
      include <- switch(effect,
        ar  = is_ar,
        cl  = !is_ar,
        all = TRUE
      )
      if (include)
        nms <- c(nms, sprintf("phi[%d,%d]", row, col))
    }
  }
  nms
}


#' Resolve Stan parameter names for beta sub-groups (internal)
#'
#' @param sd   The \code{standata} list.
#' @param type One of \code{"intercepts"} or \code{"fe"} (non-intercept fixed effects).
#' @return Character vector of Stan parameter names.
#' @keywords internal
get_beta_indices <- function(sd, type = c("intercepts", "fe")) {
  type  <- match.arg(type)
  p     <- sd$p
  n_fe  <- sd$n_fe

  nms <- character(0)
  if (type == "intercepts") {
    # Beta row 1 = intercept, across all p columns
    for (col in seq_len(p))
      nms <- c(nms, sprintf("beta[1,%d]", col))
  } else {
    # All beta rows except row 1
    if (n_fe < 2L)
      stop("No non-intercept fixed effects available.", call. = FALSE)
    for (col in seq_len(p)) {
      for (row in 2:n_fe) {
        nms <- c(nms, sprintf("beta[%d,%d]", row, col))
      }
    }
  }
  nms
}


# ---- Density estimation internals -------------------------------------------

#' Compute the SDDR using logspline (internal)
#'
#' @param draws Numeric vector of posterior draws.
#' @param prior A \code{bvarnet_prior} object.
#' @param null  Numeric scalar — null value.
#'
#' @return Named list with elements \code{BF01}, \code{post_density},
#'   \code{prior_density}.
#' @keywords internal
.compute_sddr_logspline <- function(draws, prior, null = 0) {
  # logspline with tryCatch — fall back to KDE on failure
  fit <- tryCatch(
    logspline::logspline(draws),
    error = function(e) NULL
  )

  if (!is.null(fit)) {
    post_den <- logspline::dlogspline(null, fit)
  } else {
    warning("logspline::logspline() failed; falling back to KDE (stats::density).",
            call. = FALSE)
    kde <- stats::density(draws, n = 4096)
    post_den <- stats::approx(kde$x, kde$y, xout = null, rule = 2)$y
  }

  prior_den <- eval_prior_density(prior, null)

  list(
    BF01          = post_den / prior_den,
    post_density  = post_den,
    prior_density = prior_den
  )
}


#' Compute the SDDR using a multivariate normal approximation (internal)
#'
#' @param draws_mat  Numeric matrix — S rows x d columns of posterior draws.
#' @param prior_list Named list of \code{bvarnet_prior} objects.
#' @param param_types Character vector of prior types for each column.
#' @param null_vec   Numeric vector of null values.
#' @param ridge      Ridge regularisation for the covariance matrix.
#'
#' @return Named list with elements \code{BF01}, \code{post_density},
#'   \code{prior_density}.
#' @keywords internal
.compute_sddr_mvn <- function(draws_mat, prior_list, param_types, null_vec,
                              ridge = 1e-8) {
  d <- ncol(draws_mat)
  S <- nrow(draws_mat)

  mu_hat    <- colMeans(draws_mat)
  Sigma_hat <- stats::cov(draws_mat) + ridge * diag(d)

  post_den  <- mvtnorm::dmvnorm(null_vec, mean = mu_hat, sigma = Sigma_hat)
  prior_den <- eval_joint_prior_density(prior_list, param_types, null_vec)

  list(
    BF01          = post_den / prior_den,
    post_density  = post_den,
    prior_density = prior_den
  )
}


# ---- Main workhorse ----------------------------------------------------------

#' Determine the prior type for a Stan parameter name
#' @keywords internal
.param_type <- function(param_name) {
  if (grepl("^phi\\[",   param_name)) return("phi")
  if (grepl("^beta\\[",  param_name)) return("beta")
  if (grepl("^sigma\\[", param_name)) return("sigma")
  if (grepl("^sd_u\\[",  param_name)) return("sd_u")
  if (grepl("^kappa\\[", param_name)) return("kappa")
  stop(sprintf("Cannot determine prior type for parameter '%s'.", param_name),
       call. = FALSE)
}

#' Determine whether a parameter uses a half-prior
#' @keywords internal
.is_half_prior <- function(param_name) {
  grepl("^(sigma|sd_u)\\[", param_name)
}


#' Compute the Savage-Dickey density ratio Bayes factor (internal)
#'
#' This is the internal workhorse called by \code{bf_table()}.
#' It returns a plain named list for a single test (one parameter or one joint
#' set of parameters).
#'
#' @param object     A \code{bvarnet} object.
#' @param params     Character vector of Stan parameter names.
#' @param null_value Numeric scalar or vector of null values (recycled if scalar).
#' @param method     One of \code{"auto"}, \code{"logspline"}, \code{"mvn"}.
#'
#' @return A named list with elements: \code{BF01}, \code{BF10}, \code{log_BF01},
#'   \code{post_density}, \code{prior_density}, \code{method}, \code{params},
#'   \code{null_value}.
#' @keywords internal
savage_dickey <- function(object, params, null_value = 0,
                          method = c("auto", "logspline", "mvn")) {
  stopifnot(inherits(object, "bvarnet"))
  method    <- match.arg(method)

  # --- Validate params exist in draws ---
  all_par <- dimnames(object$draws)[[3]]
  missing <- setdiff(params, all_par)
  if (length(missing) > 0)
    stop(sprintf("Parameter(s) not found in draws: %s",
                 paste(missing, collapse = ", ")), call. = FALSE)

  # --- Recycle null_value ---
  d <- length(params)
  if (length(null_value) == 1L) null_value <- rep(null_value, d)
  stopifnot(length(null_value) == d)

  # --- Block half-prior params (sigma, sd_u) ---
  # The SDDR is not valid for parameters with half-priors: boundary at 0

  # violates the interior-null requirement at the natural null, and testing
  # at positive interior values has limited practical utility.
  half_flags <- vapply(params, .is_half_prior, logical(1), USE.NAMES = FALSE)
  if (any(half_flags)) {
    bad <- params[half_flags]
    stop(sprintf(
      paste0("Cannot compute SDDR for half-prior parameter(s): %s. ",
             "Parameters with half-priors (sigma, sd_u) are not supported ",
             "because the Savage-Dickey density ratio is not valid for ",
             "distributions bounded at zero.\n\n"),
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  # --- Check for low draw count ---
  draws_3d <- object$draws
  S_total  <- dim(draws_3d)[1] * dim(draws_3d)[2]
  if (S_total < 1000)
    warning("Very few posterior draws (S = ", S_total,
            "); BF estimates may be unreliable.", call. = FALSE)

  # --- Auto-select method ---
  if (method == "auto")
    method <- if (d == 1L) "logspline" else "mvn"

  # --- Extract draws ---
  param_types <- vapply(params, .param_type, character(1), USE.NAMES = FALSE)
  priors      <- object$priors   # bvarnet_priors object

  if (method == "logspline" && d == 1L) {

    draws_vec <- .extract_draws_raw(object, params)
    res <- .compute_sddr_logspline(
      draws_vec, prior = priors[[param_types]], null = null_value
    )

  } else if (method == "mvn" || (method == "logspline" && d > 1L)) {

    if (method == "logspline" && d > 1L) {
      warning("logspline is univariate; switching to MVN for joint test.",
              call. = FALSE)
      method <- "mvn"
    }
    draws_mat <- .extract_draws_raw(object, params)
    res <- .compute_sddr_mvn(
      draws_mat, priors, param_types, null_vec = null_value
    )

  } else {
    stop("Invalid method/dimension combination.", call. = FALSE)
  }

  list(
    BF01          = res$BF01,
    BF10          = 1 / res$BF01,
    log_BF01      = log(res$BF01),
    post_density  = res$post_density,
    prior_density = res$prior_density,
    method        = method,
    params        = params,
    null_value    = null_value
  )
}


#' Extract raw posterior draws for given Stan parameter names (internal)
#'
#' Returns a numeric vector (single param) or matrix (multiple params), flattening
#' the iter x chains dimensions.
#'
#' @param object A \code{bvarnet} object.
#' @param params Character vector of Stan parameter names.
#' @return Numeric vector or matrix.
#' @keywords internal
.extract_draws_raw <- function(object, params) {
  draws_3d <- object$draws   # iter x chains x params
  idx <- match(params, dimnames(draws_3d)[[3]])
  chunk <- draws_3d[, , idx, drop = FALSE]
  # Flatten iter x chains
  dim(chunk) <- c(prod(dim(chunk)[1:2]), length(params))
  colnames(chunk) <- params
  if (length(params) == 1L) return(as.numeric(chunk))
  chunk
}


# ---- Primary export: bf_table() ----------------------------------------------

#' Compute Bayes factor table for a bvarnet model
#'
#' Computes Savage-Dickey density ratio Bayes factors for each parameter in the
#' requested subset and returns a tidy data frame.  For joint tests
#' (\code{type = "ar"}, \code{"cl"}), both per-parameter
#' (logspline) and joint (MVN) rows are returned.
#'
#' @param object     A \code{bvarnet} object returned by \code{bvar()}.
#' @param type       Character vector. Which parameter groups to test.
#'   Options: \code{"ar"} (autoregressive), \code{"cl"} (cross-lagged),
#'   \code{"intercepts"}, \code{"fe"}
#'   (non-intercept fixed effects).
#' @param lag        Integer; which lag block to use (default 1). Applies to
#'   \code{"ar"}, \code{"cl"}, and \code{"phi"}.
#' @param null_value Numeric scalar; the null hypothesis value (default 0).
#'
#' @return A data frame with columns: \code{type}, \code{predictor},
#'   \code{outcome}, \code{BF01}, \code{BF10}, \code{log_BF01},
#'   \code{post_density}, \code{prior_density}, \code{method}.
#'
#' @export
bf_table <- function(object,
                     type = c("ar", "cl", "intercepts", "fe"),
                     lag = 1L,
                     null_value = 0) {
  stopifnot(inherits(object, "bvarnet"))

  sd <- object$standata
  nm <- get_param_names(sd)
  rows <- list()

  for (tp in type) {
    tp <- match.arg(tp, c("ar", "cl", "intercepts", "fe"))

    # --- Resolve parameter names ---
    param_names <- switch(tp,
      ar          = get_phi_indices(sd, lag = lag, effect = "ar"),
      cl          = get_phi_indices(sd, lag = lag, effect = "cl"),
      intercepts  = get_beta_indices(sd, type = "intercepts"),
      fe          = get_beta_indices(sd, type = "fe")
    )

    # --- Per-parameter BFs (logspline) ---
    for (pnm in param_names) {
      res <- savage_dickey(object, params = pnm, null_value = null_value,
                           method = "logspline")
      lab <- .param_label(pnm, nm)
      rows[[length(rows) + 1L]] <- data.frame(
        type          = .type_label(tp),
        predictor     = lab$predictor,
        outcome       = lab$outcome,
        BF01          = res$BF01,
        BF10          = res$BF10,
        log_BF01      = res$log_BF01,
        post_density  = res$post_density,
        prior_density = res$prior_density,
        method        = res$method,
        stringsAsFactors = FALSE
      )
    }

    # --- Joint BF (MVN) for multi-parameter groups ---
    if (length(param_names) > 1L) {
      joint_res <- savage_dickey(object, params = param_names,
                                 null_value = null_value, method = "mvn")
      rows[[length(rows) + 1L]] <- data.frame(
        type          = paste0(.type_label(tp), " (joint)"),
        predictor     = paste0("all_", tp),
        outcome       = "\u2014",
        BF01          = joint_res$BF01,
        BF10          = joint_res$BF10,
        log_BF01      = joint_res$log_BF01,
        post_density  = joint_res$post_density,
        prior_density = joint_res$prior_density,
        method        = joint_res$method,
        stringsAsFactors = FALSE
      )
    }
  }

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}


# ---- Label helpers -----------------------------------------------------------

#' Produce human-readable labels for a Stan parameter name (internal)
#' @keywords internal
.param_label <- function(stan_name, nm) {
  # Parse "phi[r,c]" or "beta[r,c]"
  parts <- regmatches(stan_name, regexec("^(\\w+)\\[(\\d+),(\\d+)\\]$", stan_name))[[1]]
  if (length(parts) == 4L) {
    param <- parts[2]
    row_i <- as.integer(parts[3])
    col_i <- as.integer(parts[4])

    if (param == "phi") {
      predictor <- if (row_i <= length(nm$b)) nm$b[row_i] else stan_name
      outcome   <- if (col_i <= length(nm$y)) nm$y[col_i] else stan_name
    } else if (param == "beta") {
      predictor <- if (row_i <= length(nm$fe)) nm$fe[row_i] else stan_name
      outcome   <- if (col_i <= length(nm$y))  nm$y[col_i]  else stan_name
    } else {
      predictor <- stan_name
      outcome   <- stan_name
    }
  } else {
    predictor <- stan_name
    outcome   <- "\u2014"
  }

  list(predictor = predictor, outcome = outcome)
}

#' Human-readable type label (internal)
#' @keywords internal
.type_label <- function(type) {
  switch(type,
    ar         = "Autoregressive",
    cl         = "Cross-lagged",
    intercepts = "Intercept",
    fe         = "Fixed Effect",
    type
  )
}
