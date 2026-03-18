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
#' @return Named list; each element has \code{full} (all params) and
#'   \code{by_lag} (list of per-lag-block param vectors).
#' @keywords internal
get_lag_interaction_indices_by_term <- function(sd) {
  terms    <- sd$fe_interaction_terms
  fe_ic    <- sd$fe_interaction_colnames
  fe_names <- colnames(sd$X)
  p <- sd$p; K <- sd$K

  if (is.null(terms) || length(terms) == 0) {
    # Fallback: parse column names for pre-implementation fits (dev only)
    return(.parse_lag_interaction_cols(sd))
  }

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

    out[[suffix]] <- list(
      full   = unlist(lag_groups, use.names = FALSE),
      by_lag = lag_groups
    )
  }
  out
}


#' Parse lag×predictor interaction columns from X colnames (dev fallback)
#'
#' Temporary internal helper for pre-implementation fit objects that lack
#' \code{fe_interaction_terms} metadata. Scheduled for removal before alpha.
#'
#' @param sd The \code{standata} list.
#' @return Same structure as \code{get_lag_interaction_indices_by_term()}.
#' @keywords internal
.parse_lag_interaction_cols <- function(sd) {
  fe_names <- colnames(sd$X)
  lag_pat  <- "^(lag(\\d+)_[^:]+):(.+)$"
  matches  <- regmatches(fe_names, regexec(lag_pat, fe_names))

  hit <- vapply(matches, length, integer(1)) == 4L
  if (!any(hit)) return(list())

  m <- do.call(rbind, matches[hit])
  suffixes <- unique(m[, 4L])
  p <- sd$p; K <- sd$K

  out <- list()
  for (suf in suffixes) {
    suf_idx <- which(hit & vapply(matches, function(x)
      length(x) == 4L && x[4L] == suf, logical(1)))
    ks <- as.integer(vapply(matches[suf_idx], `[[`, character(1), 3L))
    beta_rows <- suf_idx

    lag_groups <- lapply(seq_len(K), function(k) {
      r_slice <- beta_rows[ks == k]
      unlist(lapply(r_slice, function(r)
        vapply(seq_len(p), function(col) sprintf("beta[%d,%d]", r, col), character(1L))
      ))
    })
    names(lag_groups) <- paste0("lag", seq_len(K), ":", suf)
    out[[suf]] <- list(full = unlist(lag_groups, use.names = FALSE), by_lag = lag_groups)
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


# ---- Primary export: bf_table() ----------------------------------------------

#' Compute Bayes factor table for a bvarnet model
#'
#' Computes Savage-Dickey density ratio Bayes factors for each parameter in the
#' requested subset and returns a tidy data frame.
#'
#' For \code{type = "fe"} and \code{"intercepts"}, the table contains three
#' levels: per-cell (logspline), per-predictor joint (MVN), and a global
#' joint-all (MVN).  For \code{type = "ar"} and \code{"cl"}, the existing
#' two-level structure (per-cell + per-type joint) is unchanged.
#'
#' \code{type = "lag_fe"} emits only grouped joint rows for lag × predictor
#' interaction terms: per-lag-block and full-term omnibus.  Per-cell rows for
#' these parameters are already included when \code{type = "fe"} is requested.
#'
#' @param object     A \code{bvarnet} object returned by \code{bvar()}.
#' @param type       Character vector. Which parameter groups to test.
#'   Options: \code{"ar"} (autoregressive), \code{"cl"} (cross-lagged),
#'   \code{"intercepts"}, \code{"fe"} (non-intercept fixed effects),
#'   \code{"lag_fe"} (lag × predictor interaction joint tests),
#'   \code{"temporal"} (joint test of all phi parameters across all lags,
#'   i.e. the entire temporal structure AR + CL, excluding covariates; when
#'   lag \ifelse{html}{\out{&times;}}{\eqn{\times}} covariate interactions
#'   are present, additional rows are emitted for per-interaction-term and
#'   full temporal + interactions omnibus joint tests).
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
    #   plus lag × covariate interactions if present
    # ------------------------------------------------------------------
    } else if (tp == "temporal") {
      K <- sd$K
      all_phi <- character(0)
      for (k in seq_len(K)) {
        all_phi <- c(all_phi, get_phi_indices(sd, lag = k, effect = "all"))
      }

      # Base temporal BF (phi only)
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
