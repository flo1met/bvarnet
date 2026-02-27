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

#' Build a complete prior specification for bvarnet
#'
#' Returns a \code{bvarnet_priors} object containing a \code{bvarnet_prior}
#' for every model parameter type.  Any argument left as \code{NULL} uses the
#' package default.
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
#' @param x A \code{bvarnet_priors} object.
#' @param ... Ignored.
#' @return \code{x} invisibly.
#' @export
print.bvarnet_priors <- function(x, ...) {
  cat("bvarnet prior specification:\n")
  half_pars <- c("sd_u", "sigma")
  for (nm in names(x)) {
    half <- nm %in% half_pars
    cat(sprintf("  %-6s ~ %s\n", nm, format(x[[nm]], half = half)))
  }
  invisible(x)
}

#' Get the default prior specification for a given model family
#'
#' A convenience wrapper around \code{set_priors()} for inspecting defaults.
#'
#' @param family One of \code{"bernoulli"}, \code{"ordinal"}, \code{"gaussian"}.
#' @return A \code{bvarnet_priors} object.
#' @export
get_default_priors <- function(family) {
  family <- match.arg(family, c("bernoulli", "ordinal", "gaussian"))
  set_priors()
}
