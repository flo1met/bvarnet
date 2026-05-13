# -----------------------------------------------------------------------------
# Prior system for bvarnet
#
# Phase 1 supports three prior families:
#   1 = Normal, 2 = Student-t, 3 = Cauchy
# Families 4 (Laplace), 7 (Exponential), 8 (Flat) are reserved and will be
# enabled once the Stan function consolidation (#include) is available.
# -----------------------------------------------------------------------------

.supported_families <- list(
  normal    = list(int = 1L, needs_scale = TRUE,  needs_df = FALSE),
  student_t = list(int = 2L, needs_scale = TRUE,  needs_df = TRUE),
  cauchy    = list(int = 3L, needs_scale = TRUE,  needs_df = FALSE)
)

#' Construct a single prior distribution
#'
#' Builds a \code{bvarnet_prior} object specifying the prior family and its
#' parameters.  Supported families in Phase 1 are \code{"normal"},
#' \code{"student_t"}, and \code{"cauchy"}.
#'
#' @param family Character. One of \code{"normal"}, \code{"student_t"},
#'   \code{"cauchy"}.
#' @param loc Location parameter (default 0).
#' @param scale Scale parameter (default 1). Must be > 0.
#' @param df Degrees of freedom for \code{"student_t"} (default 7). Must be
#'   > 0 when \code{family = "student_t"}.
#'
#' @return A \code{bvarnet_prior} S3 object.
#' @export
prior <- function(family, loc = 0, scale = 1, df = 7) {
  family <- tryCatch(
    match.arg(family, choices = names(.supported_families)),
    error = function(e) {
      stop(sprintf(
        "Unrecognised prior family '%s'. Supported: %s",
        family, paste(names(.supported_families), collapse = ", ")
      ), call. = FALSE)
    }
  )

  fam_info <- .supported_families[[family]]

  if (fam_info$needs_scale) {
    if (!is.numeric(scale) || length(scale) != 1L || is.na(scale) || scale <= 0)
      stop("'scale' must be a single positive number.", call. = FALSE)
    if (scale > 50)
      warning("Prior scale is very large. If you want an uninformative prior, consider increasing scale further or switching to a different family.", call. = FALSE)
    if (scale < 0.01)
      warning("Prior scale is very small -- this is a strongly informative prior.", call. = FALSE)
  }

  if (fam_info$needs_df) {
    if (!is.numeric(df) || length(df) != 1L || is.na(df) || df <= 0)
      stop("'df' must be a single positive number for Student-t.", call. = FALSE)
  } else {
    df <- 0
  }

  structure(
    list(
      family     = family,
      family_int = fam_info$int,
      loc        = loc,
      scale      = scale,
      df         = df,
      is_default = FALSE
    ),
    class = "bvarnet_prior"
  )
}

#' Format a bvarnet_prior for printing
#'
#' @param x A \code{bvarnet_prior} object.
#' @param half Logical; if \code{TRUE} prepends "Half-" to indicate a
#'   half-prior (used for positive-constrained parameters like sd_u and sigma).
#' @param ... Ignored.
#' @return A character string.
#' @export
format.bvarnet_prior <- function(x, half = FALSE, ...) {
  prefix <- if (isTRUE(half)) "Half-" else ""
  switch(x$family,
    normal    = sprintf("%sNormal(%g, %g)",               prefix, x$loc, x$scale),
    student_t = sprintf("%sStudent-t(%g, %g, df = %g)",  prefix, x$loc, x$scale, x$df),
    cauchy    = sprintf("%sCauchy(%g, %g)",               prefix, x$loc, x$scale)
  )
}

#' Print a bvarnet_prior
#'
#' @param x A \code{bvarnet_prior} object.
#' @param ... Passed to \code{format.bvarnet_prior()}.
#' @return \code{x} invisibly.
#' @export
print.bvarnet_prior <- function(x, ...) {
  cat(format(x, ...), "\n")
  invisible(x)
}

# Internal helper: construct a default prior (is_default = TRUE)
.default_prior <- function(family, loc, scale, df = 7) {
  p <- prior(family, loc = loc, scale = scale, df = df)
  p$is_default <- TRUE
  p
}

#' Build a prior specification object for `bvar()`
#'
#' Returns a \code{bvarnet_priors} object containing a \code{bvarnet_prior}
#' for every model parameter type.  Any argument left as \code{NULL} uses the
#' package default. 
#' Available prior distributions are:
#' - normal(loc, scale)
#' - student_t(loc, scale, df)
#' - cauchy(loc, scale)
#' For standart deviations and random effects, the prior is automatically converted to a half-prior (truncated at `loc`) in the Stan code, so the printed format reflects this.
#'
#' @param beta   Prior for fixed-effect regression coefficients.
#' @param phi    Prior for lag coefficients.
#' @param sd_u   Prior for random-effect standard deviations (half-prior).
#' @param kappa  Prior for ordinal cut-points (ordinal models only).
#' @param sigma  Prior for residual standard deviation (gaussian models only;
#'   half-prior).
#'
#' @return A \code{bvarnet_priors} S3 object.
#' @export
set_priors <- function(beta  = NULL,
                       phi   = NULL,
                       sd_u  = NULL,
                       kappa = NULL,
                       sigma = NULL) {

  defaults <- list(
    beta  = .default_prior("normal", 0, 1),
    phi   = .default_prior("normal", 0, 0.5),
    sd_u  = .default_prior("normal", 0, 1),
    kappa = .default_prior("normal", 0, 2),
    sigma = .default_prior("normal", 0, 2.5)
  )

  resolve <- function(user, default) {
    if (is.null(user)) return(default)
    if (!inherits(user, "bvarnet_prior"))
      stop("Each prior must be a bvarnet_prior object created with prior().",
           call. = FALSE)
    user
  }

  structure(
    list(
      beta  = resolve(beta,  defaults$beta),
      phi   = resolve(phi,   defaults$phi),
      sd_u  = resolve(sd_u,  defaults$sd_u),
      kappa = resolve(kappa, defaults$kappa),
      sigma = resolve(sigma, defaults$sigma)
    ),
    class = "bvarnet_priors"
  )
}

#' Print a bvarnet_priors specification
#'
#' Shows only the priors explicitly set by the user.  When no priors have
#' been overridden (all defaults), a compact note is printed instead.
#'
#' @param x A \code{bvarnet_priors} object.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.bvarnet_priors <- function(x, ...) {
  half_pars <- c("sd_u", "sigma")
  nms <- names(x)
  # Only consider actual bvarnet_prior entries (skip NULLs from filtered objects)
  nms <- nms[vapply(x[nms], function(p) inherits(p, "bvarnet_prior"), logical(1))]

  user_set <- nms[!vapply(x[nms], function(p) isTRUE(p$is_default), logical(1))]

  cat("bvarnet prior specification:\n")

  if (length(user_set) == 0L) {
    cat("  (all defaults \u2014 see ?get_default_priors)\n")
  } else {
    for (nm in user_set) {
      cat(sprintf("  %-6s ~ %s\n", nm, format(x[[nm]], half = nm %in% half_pars)))
    }
  }

  invisible(x)
}

#' Get the default prior specification for a given model family
#'
#' Returns a \code{bvarnet_priors} object showing the default priors that
#' apply to a particular model configuration. Parameters irrelevant to
#' the chosen family or model structure are omitted, so the returned object
#' reflects what the sampler will actually use.
#'
#' @param family Character (optional). One of \code{"bernoulli"},
#'   \code{"ordinal"}, \code{"gaussian"}. When \code{NULL} (the default),
#'   all parameter priors are shown.
#' @param has_re Logical. Does the model include random effects?
#'   Default \code{TRUE}. When \code{FALSE}, the \code{sd_u} prior is
#'   omitted.
#' @return A \code{bvarnet_priors} object.
#' @export
get_default_priors <- function(family = NULL, has_re = TRUE) {
  p <- set_priors()

  if (!is.null(family)) {
    family <- match.arg(family, c("bernoulli", "ordinal", "gaussian"))
    if (family != "gaussian") p$sigma <- NULL
    if (family != "ordinal")  p$kappa <- NULL
  }

  if (!isTRUE(has_re)) p$sd_u <- NULL

  p
}


# ---- Internal helpers for prior resolution in bvar() ----

#' Ensure all prior slots required by the model family exist
#'
#' If slots are missing (e.g. from a filtered \code{get_default_priors()}
#' object), they are filled with package defaults and a warning is issued.
#' Called by \code{bvar()} before passing priors to \code{to_stan_data()}.
#'
#' @param priors A \code{bvarnet_priors} object (possibly incomplete).
#' @param family_vec Named character vector of families per node.
#' @return A complete \code{bvarnet_priors} object.
#' @keywords internal
.ensure_prior_slots <- function(priors, family_vec) {
  # sd_u is always passed to Stan (even when n_re == 0)
  required <- c("beta", "phi", "sd_u")
  if (any(family_vec == "gaussian")) required <- c(required, "sigma")
  if (any(family_vec == "ordinal"))  required <- c(required, "kappa")

  present <- names(Filter(Negate(is.null), priors))
  missing <- setdiff(required, present)

  if (length(missing) > 0L) {
    defaults <- list(
      beta  = .default_prior("normal", 0, 1),
      phi   = .default_prior("normal", 0, 0.5),
      sd_u  = .default_prior("normal", 0, 1),
      kappa = .default_prior("normal", 0, 2),
      sigma = .default_prior("normal", 0, 2.5)
    )
    for (nm in missing) priors[[nm]] <- defaults[[nm]]
    warning("Prior(s) missing from input but required by this model: ",
            paste(missing, collapse = ", "), ". Using package defaults.",
            call. = FALSE)
  }
  priors
}


#' Emit user-facing warnings about prior usage and return needed-prior names
#'
#' Warns when user-set priors are not needed by the model, and messages when
#' the model uses default priors that the user did not explicitly set (only
#' if the user set at least one prior).
#'
#' @param priors A \code{bvarnet_priors} object.
#' @param family_vec Named character vector of families per node.
#' @param n_re Integer. Number of random-effect columns from the built design.
#' @return Character vector of prior names the model actually uses.
#' @keywords internal
.prior_warnings <- function(priors, family_vec, n_re) {
  needed <- c("beta", "phi")
  if (n_re > 0L) needed <- c(needed, "sd_u")
  if (any(family_vec == "gaussian")) needed <- c(needed, "sigma")
  if (any(family_vec == "ordinal"))  needed <- c(needed, "kappa")

  # Identify which slots the user explicitly set
  nms <- names(Filter(function(p) inherits(p, "bvarnet_prior"), priors))
  user_set <- nms[!vapply(priors[nms], function(p) isTRUE(p$is_default), logical(1))]

  # User-set priors the model doesn't need
  unused <- setdiff(user_set, needed)
  if (length(unused) > 0L)
    warning("Prior(s) set but not used by this model: ",
            paste(unused, collapse = ", "), ". These will be ignored.",
            call. = FALSE)

  # Defaults auto-filled for model (only warn if user set at least one)
  auto_filled <- setdiff(needed, user_set)
  if (length(user_set) > 0L && length(auto_filled) > 0L)
    message("Using default priors for: ", paste(auto_filled, collapse = ", "))

  needed
}
