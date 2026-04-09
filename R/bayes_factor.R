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

  # Guard: ordinal models have no intercept in X

  if (type == "intercepts" && !("Intercept" %in% colnames(sd$X)))
    stop(
      '`type = "intercepts"` is not valid for ordinal models: ',
      "beta row 1 is a covariate, not an intercept. ",
      "Cutpoints are stored in `kappa`.",
      call. = FALSE
    )

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


#' Group Stan beta parameter names by predictor row (internal)
#'
#' Returns a named list — one element per predictor row — each element being a
#' character vector of Stan names \code{beta[row, 1:p]} for that row.
#'
#' @param sd   The \code{standata} list.
#' @param type One of \code{"fe"} or \code{"intercepts"}.
#' @return Named list of character vectors keyed by predictor name.
#' @keywords internal
get_beta_indices_by_predictor <- function(sd, type = c("fe", "intercepts")) {
  type  <- match.arg(type)
  p     <- sd$p
  n_fe  <- sd$n_fe
  fe_names <- colnames(sd$X)
  if (is.null(fe_names)) fe_names <- paste0("fe", seq_len(n_fe))

  # Guard: ordinal models have no intercept in X
  if (type == "intercepts" && !("Intercept" %in% fe_names))
    stop(
      '`type = "intercepts"` is not valid for ordinal models: ',
      "beta row 1 is a covariate, not an intercept. ",
      "Cutpoints are stored in `kappa`.",
      call. = FALSE
    )

  rows <- if (type == "fe") seq(2L, n_fe) else 1L

  stats::setNames(
    lapply(rows, function(r) {
      vapply(seq_len(p), function(col) sprintf("beta[%d,%d]", r, col),
             character(1L))
    }),
    fe_names[rows]
  )
}


#' Group Stan beta names by lag×predictor interaction term (internal)
#'
#' @param sd The \code{standata} list.
#' @return Named list; each element has \code{full} (all params),
#'   \code{by_lag} (list of per-lag-block param vectors), \code{ar}
#'   (AR-like interaction params where lagged outcome == target outcome),
#'   and \code{cl} (CL-like interaction params where lagged outcome !=
#'   target outcome).
#' @keywords internal
get_lag_interaction_indices_by_term <- function(sd) {
  terms    <- sd$fe_interaction_terms
  fe_ic    <- sd$fe_interaction_colnames
  fe_names <- colnames(sd$X)
  p <- sd$p; K <- sd$K

  if (is.null(terms)) {
    stop("fe_interaction_terms metadata is missing from standata.",
         call. = FALSE)
  }
  if (length(terms) == 0L) return(list())

  out <- list()
  for (fac in terms) {
    has_lag <- any(fac == "lag")
    if (!has_lag) next   # plain x:x interactions handled by per-predictor plan

    suffix <- paste(fac[fac != "lag"], collapse = ":")

    # Use stored colnames; filter by suffix token
    term_ic   <- fe_ic[startsWith(fe_ic, "lag") & endsWith(fe_ic, paste0(":", suffix))]
    beta_rows <- match(term_ic, fe_names)

    if (length(term_ic) == 0L || any(is.na(beta_rows)))
      stop(sprintf("Lag-interaction columns for term '%s' not found in X.", suffix),
           call. = FALSE)

    lag_groups <- lapply(seq_len(K), function(k) {
      row_slice <- beta_rows[((k - 1L) * p + 1L):(k * p)]
      unlist(lapply(row_slice, function(r)
        vapply(seq_len(p), function(col) sprintf("beta[%d,%d]", r, col), character(1L))
      ))
    })
    names(lag_groups) <- paste0("lag", seq_len(K), ":", suffix)

    # AR/CL split: within each lag block, row_slice[j] maps to lagged
    # outcome y_j.  beta[row_for_y_j, col] is AR-like when j == col.
    ar_params <- character(0)
    cl_params <- character(0)
    for (k in seq_len(K)) {
      row_slice <- beta_rows[((k - 1L) * p + 1L):(k * p)]
      for (j in seq_len(p)) {
        r <- row_slice[j]
        for (col in seq_len(p)) {
          pname <- sprintf("beta[%d,%d]", r, col)
          if (j == col) ar_params <- c(ar_params, pname)
          else          cl_params <- c(cl_params, pname)
        }
      }
    }

    out[[suffix]] <- list(
      full   = unlist(lag_groups, use.names = FALSE),
      by_lag = lag_groups,
      ar     = ar_params,
      cl     = cl_params
    )
  }
  out
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
#'
#' @return Named list with elements \code{BF01}, \code{post_density},
#'   \code{prior_density}.
#' @keywords internal
.compute_sddr_mvn <- function(draws_mat, prior_list, param_types, null_vec) {
  d <- ncol(draws_mat)
  S <- nrow(draws_mat)

  mu_hat    <- colMeans(draws_mat)
  Sigma_hat <- stats::cov(draws_mat)

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


# ---- Parameter-block helpers for `parameter` mode ----------------------------

#' Classify a Stan parameter name into a human-readable type (internal)
#'
#' @param stan_name Single Stan parameter name (e.g. \code{"phi[1,2]"}).
#' @param sd        The \code{standata} list from a \code{bvarnet} object.
#' @return Character scalar: one of \code{"Autoregressive"},
#'   \code{"Cross-lagged"}, \code{"Intercept"}, \code{"Fixed Effect"}.
#' @keywords internal
.classify_param_type <- function(stan_name, sd) {
  parts <- regmatches(stan_name,
                      regexec("^(\\w+)\\[(\\d+),(\\d+)\\]$", stan_name))[[1]]
  if (length(parts) != 4L)
    stop(sprintf("Cannot classify parameter '%s'.", stan_name), call. = FALSE)

  param <- parts[2]
  row_i <- as.integer(parts[3])
  col_i <- as.integer(parts[4])

  if (param == "phi") {
    row_within <- ((row_i - 1L) %% sd$p) + 1L
    if (row_within == col_i) "Autoregressive" else "Cross-lagged"
  } else if (param == "beta") {
    fe_names <- colnames(sd$X)
    if (!is.null(fe_names) && row_i <= length(fe_names) &&
        fe_names[row_i] == "Intercept") {
      "Intercept"
    } else {
      "Fixed Effect"
    }
  } else {
    stop(sprintf("Unsupported parameter block '%s' in .classify_param_type().",
                 param), call. = FALSE)
  }
}


#' Build per-cell + joint BF rows for a parameter block (internal)
#'
#' @param object     A \code{bvarnet} object.
#' @param block      Character scalar: \code{"phi"} or \code{"beta"}.
#' @param null_value Numeric scalar; null hypothesis value.
#' @return A list of 1-row data frames (per-cell + optional joint).
#' @keywords internal
.bf_block_rows <- function(object, block, null_value) {
  sd <- object$standata
  nm <- get_param_names(sd)

  draws <- extract_draws(object, block)
  param_names <- colnames(draws)

  # For beta: filter out ordinal sentinel columns (beta[1,j] for ordinal nodes)
  if (block == "beta") {
    ord_indices <- which(object$family == "ordinal")
    if (length(ord_indices) > 0L) {
      sentinel_cols <- paste0("beta[1,", ord_indices, "]")
      param_names <- param_names[!param_names %in% sentinel_cols]
    }
  }

  param_names <- unique(param_names)
  rows <- list()

  # Per-cell BFs (logspline)
  for (pnm in param_names) {
    res <- savage_dickey(object, params = pnm, null_value = null_value,
                         method = "logspline")
    lab  <- .param_label(pnm, nm)
    tp   <- .classify_param_type(pnm, sd)
    rows[[length(rows) + 1L]] <- data.frame(
      type          = tp,
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

  # Joint BF (MVN) when more than one parameter
  if (length(param_names) > 1L) {
    joint_res <- savage_dickey(object, params = param_names,
                               null_value = null_value, method = "mvn")
    block_label <- switch(block,
      phi  = "Phi (joint)",
      beta = "Beta (joint)"
    )
    rows[[length(rows) + 1L]] <- data.frame(
      type          = block_label,
      predictor     = paste0("all_", block),
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

  rows
}


# ---- Primary export: bf_table() ----------------------------------------------

#' Compute Bayes factor table for a bvarnet model
#'
#' Computes Savage-Dickey density ratio Bayes factors for each parameter in the
#' requested subset and returns a tidy data frame.
#'
#' **Type mode** (\code{parameter = NULL}, default): uses semantic groupings.
#' For \code{type = "fe"} and \code{"intercepts"}, the table contains three
#' levels: per-cell (logspline), per-predictor joint (MVN), and a global
#' joint-all (MVN).  For \code{type = "ar"} and \code{"cl"}, the existing
#' two-level structure (per-cell + per-type joint) is unchanged.
#'
#' \code{type = "lag_fe"} emits only grouped joint rows for lag × predictor
#' interaction terms: per-lag-block and full-term omnibus.  Per-cell rows for
#' these parameters are already included when \code{type = "fe"} is requested.
#'
#' **Parameter mode** (\code{parameter} non-NULL): returns per-cell BFs plus
#' a joint BF for each requested Stan parameter block.  Mutually exclusive
#' with non-default \code{type} or \code{lag}.
#'
#' @param object     A \code{bvarnet} object returned by \code{bvar()}.
#' @param type       Character vector or \code{"all"} (default).  Which
#'   parameter groups to test.
#'   Options: \code{"ar"} (autoregressive), \code{"cl"} (cross-lagged),
#'   \code{"intercepts"}, \code{"fe"} (non-intercept fixed effects),
#'   \code{"lag_fe"} (lag × predictor interaction joint tests),
#'   \code{"temporal"} (joint test of all phi parameters across all lags,
#'   i.e. the entire temporal structure AR + CL, excluding covariates;
#'   additionally emits separate joint rows for AR-only and CL-only
#'   components; when lag \ifelse{html}{\out{&times;}}{\eqn{\times}}
#'   covariate interactions are present, additional rows are emitted for
#'   per-interaction-term and AR-only / CL-only interaction sub-tests,
#'   plus a full temporal + interactions omnibus).
#'   \code{"all"} auto-selects all applicable types (skips
#'   \code{"intercepts"} for ordinal models and \code{"lag_fe"} when
#'   no lag interactions exist).
#'   Ignored when \code{parameter} is non-NULL.
#' @param lag        Integer; which lag block to use (default 1). Applies to
#'   \code{"ar"}, \code{"cl"}, and \code{"phi"} in type mode only.
#'   Ignored (and must be default) when \code{parameter} is non-NULL.
#' @param null_value Numeric scalar; the null hypothesis value (default 0).
#' @param parameter  Character vector or \code{NULL} (default).  Stan
#'   parameter block names to compute per-cell + joint BFs for.  Currently
#'   supported: \code{"phi"} and \code{"beta"}.  \code{"sigma"} and
#'   \code{"sd_u"} are rejected (half-prior makes SDDR invalid).
#'   \code{"kappa"} is not yet supported (ordered constraint requires
#'   specialised density evaluation).  Mutually exclusive with non-default
#'   \code{type} or \code{lag}.
#'
#' @return A data frame with columns: \code{type}, \code{predictor},
#'   \code{outcome}, \code{BF01}, \code{BF10}, \code{log_BF01},
#'   \code{post_density}, \code{prior_density}, \code{method}.
#'
#' @export
bf_table <- function(object,
                     type = "all",
                     lag = 1L,
                     null_value = 0,
                     parameter = NULL) {
  stopifnot(inherits(object, "bvarnet"))

  # ---- Parameter mode (early return) ----
  if (!is.null(parameter)) {
    # Mutual exclusivity
    if (!identical(type, "all"))
      stop("'parameter' and 'type' are mutually exclusive. ",
           "Use one or the other, not both.", call. = FALSE)
    if (!identical(lag, 1L))
      stop("'parameter' and 'lag' are mutually exclusive. ",
           "In parameter mode all lags are included automatically.",
           call. = FALSE)

    parameter <- unique(parameter)

    # Validate requested blocks
    half_prior_blocks <- c("sigma", "sd_u")
    bad_half <- intersect(parameter, half_prior_blocks)
    if (length(bad_half) > 0L)
      stop(sprintf(
        "Cannot compute SDDR for half-prior parameter(s): %s. ",
        paste(bad_half, collapse = ", ")
      ), "Parameters with half-priors (sigma, sd_u) are not supported ",
      "because the Savage-Dickey density ratio is not valid for ",
      "distributions bounded at zero.", call. = FALSE)

    if ("kappa" %in% parameter)
      stop("Cannot compute SDDR for 'kappa'. ",
           "Kappa has an ordered constraint in Stan, so the unconstrained ",
           "prior density used by the Savage-Dickey method is incorrect. ",
           "This feature may be added in a future version.", call. = FALSE)

    valid_blocks <- c("phi", "beta")
    bad <- setdiff(parameter, valid_blocks)
    if (length(bad) > 0L)
      stop(sprintf("Unsupported parameter block(s): %s. Supported: %s.",
                   paste(bad, collapse = ", "),
                   paste(valid_blocks, collapse = ", ")), call. = FALSE)

    rows <- list()
    for (blk in parameter)
      rows <- c(rows, .bf_block_rows(object, blk, null_value))

    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  # ---- Type mode (existing logic) ----
  sd <- object$standata
  nm <- get_param_names(sd)

  # Resolve "all" to applicable types
  if (identical(type, "all")) {
    type <- c("ar", "cl", "temporal")
    if (sd$n_fe >= 2L)
      type <- c(type, "fe")
    if (any(object$family != "ordinal") && "Intercept" %in% colnames(sd$X))
      type <- c("intercepts", type)
    lag_terms <- get_lag_interaction_indices_by_term(sd)
    if (length(lag_terms) > 0L)
      type <- c(type, "lag_fe")
  }

  rows <- list()

  for (tp in type) {
    tp <- match.arg(tp, c("ar", "cl", "intercepts", "fe", "lag_fe", "temporal"))

    # ------------------------------------------------------------------
    # AR / CL — unchanged two-level structure
    # ------------------------------------------------------------------
    if (tp %in% c("ar", "cl")) {
      param_names <- switch(tp,
        ar = get_phi_indices(sd, lag = lag, effect = "ar"),
        cl = get_phi_indices(sd, lag = lag, effect = "cl")
      )

      # Per-parameter BFs (logspline)
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

      # Joint BF (MVN)
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

    # ------------------------------------------------------------------
    # FE / INTERCEPTS — three-phase structure
    # ------------------------------------------------------------------
    } else if (tp %in% c("fe", "intercepts")) {
      param_names <- get_beta_indices(sd, type = tp)
      by_pred     <- get_beta_indices_by_predictor(sd, type = tp)
      type_base   <- .type_label(tp)

      # Phase A — Per-cell BFs (logspline)
      for (pnm in param_names) {
        res <- savage_dickey(object, params = pnm, null_value = null_value,
                             method = "logspline")
        lab <- .param_label(pnm, nm)
        rows[[length(rows) + 1L]] <- data.frame(
          type          = type_base,
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

      # Phase B — Per-predictor joint BFs (MVN), only when p > 1
      for (pred_name in names(by_pred)) {
        pvec <- by_pred[[pred_name]]
        if (length(pvec) > 1L) {
          jres <- savage_dickey(object, params = pvec,
                                null_value = null_value, method = "mvn")
          rows[[length(rows) + 1L]] <- data.frame(
            type          = paste0(type_base, " (joint)"),
            predictor     = pred_name,
            outcome       = "\u2014",
            BF01          = jres$BF01,
            BF10          = jres$BF10,
            log_BF01      = jres$log_BF01,
            post_density  = jres$post_density,
            prior_density = jres$prior_density,
            method        = jres$method,
            stringsAsFactors = FALSE
          )
        }
      }

      # Phase C — Global joint-all BF (MVN), only when > 1 predictor row
      if (length(by_pred) > 1L) {
        all_params <- unlist(by_pred, use.names = FALSE)
        gres <- savage_dickey(object, params = all_params,
                              null_value = null_value, method = "mvn")
        rows[[length(rows) + 1L]] <- data.frame(
          type          = paste0(type_base, " (joint all)"),
          predictor     = paste0("all_", tp),
          outcome       = "\u2014",
          BF01          = gres$BF01,
          BF10          = gres$BF10,
          log_BF01      = gres$log_BF01,
          post_density  = gres$post_density,
          prior_density = gres$prior_density,
          method        = gres$method,
          stringsAsFactors = FALSE
        )
      }

    # ------------------------------------------------------------------
    # LAG_FE — grouped joint rows only
    # ------------------------------------------------------------------
    } else if (tp == "lag_fe") {
      term_groups <- get_lag_interaction_indices_by_term(sd)

      if (length(term_groups) == 0L)
        stop(
          "No lag interaction columns found in standata$X. ",
          "Was the model fitted with fe_interactions containing 'lag' terms?",
          call. = FALSE
        )

      for (suffix in names(term_groups)) {
        tg <- term_groups[[suffix]]

        # Per-lag-block joint BFs
        for (lag_name in names(tg$by_lag)) {
          lg_params <- tg$by_lag[[lag_name]]
          if (length(lg_params) > 0L) {
            lres <- savage_dickey(object, params = lg_params,
                                  null_value = null_value,
                                  method = if (length(lg_params) == 1L) "logspline" else "mvn")
            rows[[length(rows) + 1L]] <- data.frame(
              type          = "Lag Interaction (per lag)",
              predictor     = lag_name,
              outcome       = "\u2014",
              BF01          = lres$BF01,
              BF10          = lres$BF10,
              log_BF01      = lres$log_BF01,
              post_density  = lres$post_density,
              prior_density = lres$prior_density,
              method        = lres$method,
              stringsAsFactors = FALSE
            )
          }
        }

        # Full-term omnibus joint BF
        if (length(tg$full) > 1L) {
          ores <- savage_dickey(object, params = tg$full,
                                null_value = null_value, method = "mvn")
          rows[[length(rows) + 1L]] <- data.frame(
            type          = "Lag Interaction (joint)",
            predictor     = suffix,
            outcome       = "\u2014",
            BF01          = ores$BF01,
            BF10          = ores$BF10,
            log_BF01      = ores$log_BF01,
            post_density  = ores$post_density,
            prior_density = ores$prior_density,
            method        = ores$method,
            stringsAsFactors = FALSE
          )
        }
      }

    # ------------------------------------------------------------------
    # TEMPORAL — joint BF over entire phi (AR + CL, all lags)
    #   plus separate AR-only and CL-only joint BFs
    #   plus lag × covariate interactions if present
    # ------------------------------------------------------------------
    } else if (tp == "temporal") {
      K <- sd$K
      all_phi <- character(0)
      all_ar  <- character(0)
      all_cl  <- character(0)
      for (k in seq_len(K)) {
        all_phi <- c(all_phi, get_phi_indices(sd, lag = k, effect = "all"))
        all_ar  <- c(all_ar,  get_phi_indices(sd, lag = k, effect = "ar"))
        all_cl  <- c(all_cl,  get_phi_indices(sd, lag = k, effect = "cl"))
      }

      # Base temporal BF (phi only — AR + CL combined)
      if (length(all_phi) == 1L) {
        tres <- savage_dickey(object, params = all_phi,
                              null_value = null_value, method = "logspline")
      } else {
        tres <- savage_dickey(object, params = all_phi,
                              null_value = null_value, method = "mvn")
      }
      rows[[length(rows) + 1L]] <- data.frame(
        type          = "Temporal (joint)",
        predictor     = "all_phi",
        outcome       = "\u2014",
        BF01          = tres$BF01,
        BF10          = tres$BF10,
        log_BF01      = tres$log_BF01,
        post_density  = tres$post_density,
        prior_density = tres$prior_density,
        method        = tres$method,
        stringsAsFactors = FALSE
      )

      # AR-only joint BF
      if (length(all_ar) > 0L) {
        if (length(all_ar) == 1L) {
          ar_res <- savage_dickey(object, params = all_ar,
                                  null_value = null_value, method = "logspline")
        } else {
          ar_res <- savage_dickey(object, params = all_ar,
                                  null_value = null_value, method = "mvn")
        }
        rows[[length(rows) + 1L]] <- data.frame(
          type          = "Temporal AR (joint)",
          predictor     = "all_ar",
          outcome       = "\u2014",
          BF01          = ar_res$BF01,
          BF10          = ar_res$BF10,
          log_BF01      = ar_res$log_BF01,
          post_density  = ar_res$post_density,
          prior_density = ar_res$prior_density,
          method        = ar_res$method,
          stringsAsFactors = FALSE
        )
      }

      # CL-only joint BF (empty when p = 1)
      if (length(all_cl) > 0L) {
        if (length(all_cl) == 1L) {
          cl_res <- savage_dickey(object, params = all_cl,
                                  null_value = null_value, method = "logspline")
        } else {
          cl_res <- savage_dickey(object, params = all_cl,
                                  null_value = null_value, method = "mvn")
        }
        rows[[length(rows) + 1L]] <- data.frame(
          type          = "Temporal CL (joint)",
          predictor     = "all_cl",
          outcome       = "\u2014",
          BF01          = cl_res$BF01,
          BF10          = cl_res$BF10,
          log_BF01      = cl_res$log_BF01,
          post_density  = cl_res$post_density,
          prior_density = cl_res$prior_density,
          method        = cl_res$method,
          stringsAsFactors = FALSE
        )
      }

      # Lag × covariate interaction terms (if any)
      term_groups <- get_lag_interaction_indices_by_term(sd)
      if (length(term_groups) > 0L) {
        lag_int_params <- unlist(
          lapply(term_groups, `[[`, "full"), use.names = FALSE
        )

        # Per-term joint BFs
        for (suffix in names(term_groups)) {
          tg_params <- term_groups[[suffix]]$full
          if (length(tg_params) > 0L) {
            # Full interaction (AR+CL)
            ires <- savage_dickey(
              object, params = tg_params,
              null_value = null_value,
              method = if (length(tg_params) == 1L) "logspline" else "mvn"
            )
            rows[[length(rows) + 1L]] <- data.frame(
              type          = "Temporal Interaction (joint)",
              predictor     = suffix,
              outcome       = "\u2014",
              BF01          = ires$BF01,
              BF10          = ires$BF10,
              log_BF01      = ires$log_BF01,
              post_density  = ires$post_density,
              prior_density = ires$prior_density,
              method        = ires$method,
              stringsAsFactors = FALSE
            )

            # AR × interaction: beta params where lagged outcome == target
            ar_int <- term_groups[[suffix]]$ar
            if (length(ar_int) > 0L) {
              arx_res <- savage_dickey(
                object, params = ar_int,
                null_value = null_value,
                method = if (length(ar_int) == 1L) "logspline" else "mvn"
              )
              rows[[length(rows) + 1L]] <- data.frame(
                type          = "Temporal AR \u00d7 Interaction (joint)",
                predictor     = suffix,
                outcome       = "\u2014",
                BF01          = arx_res$BF01,
                BF10          = arx_res$BF10,
                log_BF01      = arx_res$log_BF01,
                post_density  = arx_res$post_density,
                prior_density = arx_res$prior_density,
                method        = arx_res$method,
                stringsAsFactors = FALSE
              )
            }
            # CL × interaction: beta params where lagged outcome != target
            cl_int <- term_groups[[suffix]]$cl
            if (length(cl_int) > 0L) {
              clx_res <- savage_dickey(
                object, params = cl_int,
                null_value = null_value,
                method = if (length(cl_int) == 1L) "logspline" else "mvn"
              )
              rows[[length(rows) + 1L]] <- data.frame(
                type          = "Temporal CL \u00d7 Interaction (joint)",
                predictor     = suffix,
                outcome       = "\u2014",
                BF01          = clx_res$BF01,
                BF10          = clx_res$BF10,
                log_BF01      = clx_res$log_BF01,
                post_density  = clx_res$post_density,
                prior_density = clx_res$prior_density,
                method        = clx_res$method,
                stringsAsFactors = FALSE
              )
            }
          }
        }

        # Full temporal + interactions omnibus
        all_temporal <- c(all_phi, lag_int_params)
        ares <- savage_dickey(object, params = all_temporal,
                              null_value = null_value, method = "mvn")
        rows[[length(rows) + 1L]] <- data.frame(
          type          = "Temporal + Interactions (joint)",
          predictor     = "all_temporal",
          outcome       = "\u2014",
          BF01          = ares$BF01,
          BF10          = ares$BF10,
          log_BF01      = ares$log_BF01,
          post_density  = ares$post_density,
          prior_density = ares$prior_density,
          method        = ares$method,
          stringsAsFactors = FALSE
        )
      }
    }  }

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
