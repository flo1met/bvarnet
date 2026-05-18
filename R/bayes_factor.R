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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
get_beta_indices <- function(sd, type = c("intercepts", "fe")) {
  type  <- match.arg(type)
  p     <- sd$p
  n_fe  <- sd$n_fe
  has_intercept <- "Intercept" %in% colnames(sd$X)

  # Guard: ordinal models have no intercept in X

  if (type == "intercepts" && !has_intercept)
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
    # FE rows start at 2 when Intercept present, 1 otherwise (pure ordinal)
    fe_start <- if (has_intercept) 2L else 1L
    if (n_fe < fe_start)
      stop("No non-intercept fixed effects available.", call. = FALSE)
    for (col in seq_len(p)) {
      for (row in fe_start:n_fe) {
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
#' @noRd
get_beta_indices_by_predictor <- function(sd, type = c("fe", "intercepts")) {
  type  <- match.arg(type)
  p     <- sd$p
  n_fe  <- sd$n_fe
  fe_names <- colnames(sd$X)
  if (is.null(fe_names)) fe_names <- paste0("fe", seq_len(n_fe))
  has_intercept <- "Intercept" %in% fe_names

  # Guard: ordinal models have no intercept in X
  if (type == "intercepts" && !has_intercept)
    stop(
      '`type = "intercepts"` is not valid for ordinal models: ',
      "beta row 1 is a covariate, not an intercept. ",
      "Cutpoints are stored in `kappa`.",
      call. = FALSE
    )

  fe_start <- if (has_intercept) 2L else 1L
  if (type == "fe" && n_fe < fe_start)
    stop("No non-intercept fixed effects available.", call. = FALSE)
  rows <- if (type == "fe") seq(fe_start, n_fe) else 1L

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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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
#' @noRd
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

# ---- Variable-filtering helpers ----------------------------------------------

#' Validate variable names and classify into network/covariate (internal)
#'
#' @param sd       The \code{standata} list from a \code{bvarnet} object.
#' @param variable Character vector of variable names to look up.
#' @return A list with components \code{y_idx} (integer indices into Y columns,
#'   or \code{NULL}), \code{x_names} (character vector of matched covariate
#'   names, or \code{NULL}), and \code{has_y}/\code{has_x} logicals.
#' @keywords internal
#' @noRd
.classify_variable <- function(sd, variable) {
  y_names <- colnames(sd$Y)
  x_names <- colnames(sd$X)
  # Base covariates: exclude Intercept and interaction columns
  x_base <- x_names[!grepl(":", x_names) & x_names != "Intercept"]

  var_y <- intersect(variable, y_names)
  var_x <- intersect(variable, x_base)
  bad   <- setdiff(variable, c(y_names, x_base))
  if (length(bad) > 0L)
    stop(sprintf(
      "Unknown variable(s): %s. Available network variables: %s; covariates: %s.",
      paste(bad, collapse = ", "),
      paste(y_names, collapse = ", "),
      if (length(x_base)) paste(x_base, collapse = ", ") else "(none)"
    ), call. = FALSE)

  list(
    y_idx  = if (length(var_y)) match(var_y, y_names) else NULL,
    x_names = if (length(var_x)) var_x else NULL,
    has_y  = length(var_y) > 0L,
    has_x  = length(var_x) > 0L
  )
}

#' Filter phi indices to effects FROM specific variables (internal)
#'
#' Keeps only \code{phi[row, col]} entries where the lagged predictor
#' (\code{row_within}) is in \code{var_idx}.
#'
#' @param sd      The \code{standata} list.
#' @param lag     Integer; which lag block.
#' @param effect  One of \code{"ar"}, \code{"cl"}, \code{"all"}.
#' @param var_idx Integer vector of variable column positions.
#' @return Character vector of filtered Stan parameter names.
#' @keywords internal
#' @noRd
.filter_phi_by_variable <- function(sd, lag, effect, var_idx) {
  all_params <- get_phi_indices(sd, lag = lag, effect = effect)
  p <- sd$p
  keep <- vapply(all_params, function(pnm) {
    parts <- regmatches(pnm, regexec("^phi\\[(\\d+),(\\d+)\\]$", pnm))[[1]]
    row_i <- as.integer(parts[2])
    row_within <- ((row_i - 1L) %% p) + 1L
    row_within %in% var_idx
  }, logical(1L))
  all_params[keep]
}

#' Filter lag-interaction indices to effects FROM specific variables (internal)
#'
#' @param sd      The \code{standata} list.
#' @param var_idx Integer vector of variable column positions.
#' @return Same structure as \code{get_lag_interaction_indices_by_term()} but
#'   with only parameters where the lagged variable is in \code{var_idx}.
#' @keywords internal
#' @noRd
.filter_lag_interaction_by_variable <- function(sd, var_idx) {
  term_groups <- get_lag_interaction_indices_by_term(sd)
  p <- sd$p; K <- sd$K

  for (suffix in names(term_groups)) {
    tg <- term_groups[[suffix]]

    # Within each lag block of p rows, position j maps to lagged variable j.
    # Each position j generates p params (one per outcome col).
    # In by_lag[[k]] (length p*p), entries ((j-1)*p + 1):(j*p) → variable j.
    new_by_lag <- lapply(seq_len(K), function(k) {
      lg <- tg$by_lag[[k]]
      keep_idx <- integer(0)
      for (j in var_idx) {
        start <- (j - 1L) * p + 1L
        end   <- j * p
        if (end <= length(lg))
          keep_idx <- c(keep_idx, start:end)
      }
      lg[keep_idx]
    })
    names(new_by_lag) <- names(tg$by_lag)

    new_full <- unlist(new_by_lag, use.names = FALSE)
    # Reuse original AR/CL classification
    new_ar <- intersect(new_full, tg$ar)
    new_cl <- intersect(new_full, tg$cl)

    term_groups[[suffix]] <- list(
      full   = new_full,
      by_lag = new_by_lag,
      ar     = new_ar,
      cl     = new_cl
    )
  }

  # Remove terms that are now empty after filtering
  term_groups[vapply(term_groups, function(tg) length(tg$full) > 0L, logical(1L))]
}

#' Filter lag-interaction term groups to a specific set of covariates (internal)
#'
#' @param sd      The \code{standata} list.
#' @param x_names Character vector of covariate names to keep.
#' @return Filtered version of \code{get_lag_interaction_indices_by_term()},
#'   retaining only terms whose non-lag suffix matches \code{x_names}.
#' @keywords internal
#' @noRd
.filter_lag_interaction_by_covariate <- function(sd, x_names) {
  term_groups <- get_lag_interaction_indices_by_term(sd)
  # Term keys are the covariate suffix (e.g., "x_1" from c("lag","x_1"))
  term_groups[names(term_groups) %in% x_names]
}

#' Filter FE by_pred list to specific covariates (internal)
#'
#' @param by_pred Named list from \code{get_beta_indices_by_predictor()}.
#' @param x_names Character vector of covariate names to keep.
#' @return Filtered list retaining only predictors in \code{x_names}.
#' @keywords internal
#' @noRd
.filter_fe_by_covariate <- function(by_pred, x_names) {
  by_pred[names(by_pred) %in% x_names]
}


# ---- Primary export: bf_table() ----------------------------------------------#' Compute Bayes factor table for a bvarnet model
#'
#' Compute Savege-Dickey Bayes factors
#'
#' Computes Savage-Dickey density ratio Bayes factors for each (requested set of) parameter in the model.  
#' By default, all applicable parameters are tested and returned in a tidy data frame.  
#' The \code{type} argument controls which parameter groups are included; the \code{variable} argument can be used to filter to effects involving specific variables.
#' The \code{log_BF10} argument allows including the natural log of the Bayes factor in the output, and \code{round} controls numeric rounding of the results.
#'
#' @param object     A \code{bvarnet} object returned by \code{bvar()}.
#' @param type Character vector specifying which parameter groups to test.
#'   Use \code{"all"} (default) to include all applicable groups automatically.
#'   Available options:
#'   \describe{
#'     \item{\code{"ar"}}{Autoregressive effects (self-loops). Per-cell BFs
#'       for the lag specified by \code{lag}, plus a joint BF.}
#'     \item{\code{"cl"}}{Cross-lagged effects. Same structure as \code{"ar"}.}
#'     \item{\code{"intercepts"}}{Intercept parameters. Skipped automatically
#'       for ordinal outcomes.}
#'     \item{\code{"fe"}}{Non-intercept fixed effects (covariates).}
#'     \item{\code{"lag_fe"}}{Joint BFs for lag \eqn{\times} covariate
#'       interaction terms. Only available when the model was fitted with
#'       \code{fe_interactions} containing lag terms.}
#'     \item{\code{"temporal"}}{Joint BF for the entire temporal structure
#'       (all AR + CL parameters across all lags). When lag \eqn{\times}
#'       covariate interactions are present, additional omnibus rows are
#'       included.}
#'   }
#' @param lag        Integer; which lag block to use (default 1). Applies to
#'   \code{"ar"} and \code{"cl"} types.
#' @param null_value Numeric scalar; the null hypothesis value (default 0).
#' @param variable   Character vector or \code{NULL} (default).  One or more
#'   variable names.  When set, only effects involving these variables are
#'   included.
#' @param log_BF10 Logical; if \code{TRUE}, an additional \code{log_BF10}
#'   column (natural log of \code{BF10}) is appended to the output.
#'   Default is \code{FALSE}.
#' @param round Integer or \code{NULL}; number of decimal places to round
#'   numeric output columns.  Default is \code{5}.  Set to \code{NULL} to
#'   disable rounding.
#'
#' @return A data frame with columns: \code{type}, \code{predictor},
#'   \code{outcome}, \code{BF10} (and optionally \code{log_BF10}).
#'
#' @export
bf_table <- function(object,
                     type = "all",
                     lag = 1L,
                     null_value = 0,
                     variable = NULL,
                     log_BF10 = FALSE,
                     round = 5L) {
  stopifnot(inherits(object, "bvarnet"))

  # ---- Variable validation ----
  var_idx <- NULL
  var_x   <- NULL
  if (!is.null(variable)) {
    stopifnot(is.character(variable), length(variable) >= 1L)
    variable <- unique(variable)
    sd_tmp <- object$standata
    vclass <- .classify_variable(sd_tmp, variable)
    var_idx <- vclass$y_idx
    var_x   <- vclass$x_names

    # Cannot combine covariate-only variable with intercepts
    explicit_types <- if (!identical(type, "all")) type else character(0)
    if ("intercepts" %in% explicit_types)
      stop("'variable' cannot be combined with type = \"intercepts\".",
           call. = FALSE)
    # Cannot combine network-only variable with fe (no covariate to filter)
    if ("fe" %in% explicit_types && !vclass$has_x)
      stop("'variable' with type = \"fe\" requires at least one covariate name ",
           "(from x_cols). The variables you supplied are network variables; ",
           "use type = \"ar\", \"cl\", or \"temporal\" instead.",
           call. = FALSE)
    # Cannot combine covariate-only variable with ar/cl (phi is not filtered)
    if (any(c("ar", "cl") %in% explicit_types) && !vclass$has_y)
      stop("'variable' with type = \"ar\" or \"cl\" requires at least one ",
           "network variable (from colnames(standata$Y)). The variables you ",
           "supplied are covariates; use type = \"fe\", \"lag_fe\", or ",
           "\"temporal\" instead.",
           call. = FALSE)
  }

  # ---- Type mode (existing logic) ----
  sd <- object$standata
  nm <- get_param_names(sd)
  has_intercept <- "Intercept" %in% colnames(sd$X)
  # Minimum n_fe to have non-intercept FE: 2 with intercept, 1 without
  fe_min <- if (has_intercept) 2L else 1L

  # Resolve "all" to applicable types
  if (identical(type, "all")) {
    if (!is.null(var_idx) || !is.null(var_x)) {
      # Variable mode: start with types relevant to what was requested
      type <- character(0)
      if (!is.null(var_idx))
        type <- c("ar", "cl", "temporal")
      if (!is.null(var_x) && sd$n_fe >= fe_min)
        type <- c(type, "fe")
      if (!("temporal" %in% type) && !is.null(var_x))
        type <- c(type, "temporal")  # for interaction rows
    } else {
      type <- c("ar", "cl", "temporal")
      if (sd$n_fe >= fe_min)
        type <- c(type, "fe")
      if (has_intercept && any(object$family != "ordinal"))
        type <- c("intercepts", type)
    }
    lag_terms <- get_lag_interaction_indices_by_term(sd)
    if (length(lag_terms) > 0L) {
      # In variable mode, only add lag_fe if filtered terms would be non-empty
      if (!is.null(var_idx) || !is.null(var_x)) {
        filtered <- lag_terms
        if (!is.null(var_idx))
          filtered <- .filter_lag_interaction_by_variable(sd, var_idx)
        if (!is.null(var_x))
          filtered <- filtered[names(filtered) %in% var_x]
        if (length(filtered) > 0L)
          type <- c(type, "lag_fe")
      } else {
        type <- c(type, "lag_fe")
      }
    }
    type <- unique(type)
  }

  rows <- list()

  for (tp in type) {
    tp <- match.arg(tp, c("ar", "cl", "intercepts", "fe", "lag_fe", "temporal"))

    # ------------------------------------------------------------------
    # AR / CL — unchanged two-level structure
    # ------------------------------------------------------------------
    if (tp %in% c("ar", "cl")) {
      param_names <- if (!is.null(var_idx)) {
        .filter_phi_by_variable(sd, lag = lag, effect = tp, var_idx = var_idx)
      } else {
        switch(tp,
          ar = get_phi_indices(sd, lag = lag, effect = "ar"),
          cl = get_phi_indices(sd, lag = lag, effect = "cl")
        )
      }

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

      # Filter to requested covariates when variable has x_names
      if (tp == "fe" && !is.null(var_x)) {
        by_pred <- .filter_fe_by_covariate(by_pred, var_x)
        param_names <- unlist(by_pred, use.names = FALSE)
      }

      # Filter out ordinal beta[1,j] sentinels for intercepts in mixed-family
      if (tp == "intercepts") {
        ord_indices <- which(object$family == "ordinal")
        if (length(ord_indices) > 0L) {
          sentinel_cols <- paste0("beta[1,", ord_indices, "]")
          param_names <- param_names[!param_names %in% sentinel_cols]
          by_pred <- lapply(by_pred, function(pvec) {
            pvec[!pvec %in% sentinel_cols]
          })
          by_pred <- by_pred[lengths(by_pred) > 0L]
        }
      }

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
      # Filter by network variable (source variable within interactions)
      if (!is.null(var_idx))
        term_groups <- .filter_lag_interaction_by_variable(sd, var_idx)
      # Filter by covariate (which interaction terms to include)
      if (!is.null(var_x))
        term_groups <- term_groups[names(term_groups) %in% var_x]

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

        # Full-term omnibus joint BF (skip at K=1 — identical to per-lag row)
        if (length(tg$full) > 1L && sd$K > 1L) {
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
      # When only covariates (no network vars), skip phi-only temporal rows
      covariate_only <- !is.null(var_x) && is.null(var_idx)

      all_phi <- character(0)
      all_ar  <- character(0)
      all_cl  <- character(0)
      if (!covariate_only) {
        for (k in seq_len(K)) {
          if (!is.null(var_idx)) {
            all_phi <- c(all_phi, .filter_phi_by_variable(sd, lag = k, effect = "all", var_idx = var_idx))
            all_ar  <- c(all_ar,  .filter_phi_by_variable(sd, lag = k, effect = "ar",  var_idx = var_idx))
            all_cl  <- c(all_cl,  .filter_phi_by_variable(sd, lag = k, effect = "cl",  var_idx = var_idx))
          } else {
            all_phi <- c(all_phi, get_phi_indices(sd, lag = k, effect = "all"))
            all_ar  <- c(all_ar,  get_phi_indices(sd, lag = k, effect = "ar"))
            all_cl  <- c(all_cl,  get_phi_indices(sd, lag = k, effect = "cl"))
          }
        }
      }

      # Base temporal BF (phi only — AR + CL combined)
      if (length(all_phi) > 0L) {
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
      }

      # AR-only joint BF (skip at K=1 to avoid duplicating ar/cl type rows)
      if (length(all_ar) > 0L && K > 1L) {
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

      # CL-only joint BF (skip at K=1 or p=1)
      if (length(all_cl) > 0L && K > 1L) {
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
      if (!is.null(var_idx))
        term_groups <- .filter_lag_interaction_by_variable(sd, var_idx)
      if (!is.null(var_x))
        term_groups <- term_groups[names(term_groups) %in% var_x]
      if (length(term_groups) > 0L) {
        lag_int_params <- unlist(
          lapply(term_groups, `[[`, "full"), use.names = FALSE
        )

        # Per-term sub-tests (AR×interaction, CL×interaction)
        for (suffix in names(term_groups)) {
          tg_params <- term_groups[[suffix]]$full
          if (length(tg_params) > 0L) {
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

        # Full temporal + interactions omnibus (only when both phi and
        # interaction params exist; otherwise it duplicates another row)
        all_temporal <- c(all_phi, lag_int_params)
        if (length(all_phi) > 0L && length(lag_int_params) > 0L) {
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
      }
    }  }

  if (length(rows) == 0L) {
    out <- data.frame(
      type = character(0), predictor = character(0), outcome = character(0),
      BF10 = numeric(0), stringsAsFactors = FALSE
    )
    if (isTRUE(log_BF10)) out$log_BF10 <- numeric(0)
    return(out)
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out <- out[, c("type", "predictor", "outcome", "BF10")]
  if (isTRUE(log_BF10)) out$log_BF10 <- log(out$BF10)
  if (!is.null(round)) {
    num_cols <- sapply(out, is.numeric)
    out[num_cols] <- lapply(out[num_cols], round, digits = round)
  }
  out
}


# ---- Label helpers -----------------------------------------------------------

#' Produce human-readable labels for a Stan parameter name (internal)
#' @keywords internal
#' @noRd
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
#' @noRd
.type_label <- function(type) {
  switch(type,
    ar         = "Autoregressive",
    cl         = "Cross-lagged",
    intercepts = "Intercept",
    fe         = "Fixed Effect",
    type
  )
}
