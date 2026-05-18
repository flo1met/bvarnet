## ---- family vector helpers (D2) ----

#' Does any node have a given family?
#' @keywords internal
#' @noRd
.family_has <- function(object, fam) any(object$family == fam)

#' Which nodes have a given family?
#' @keywords internal
#' @noRd
.family_which <- function(object, fam) which(object$family == fam)

#' Is the model mixed-family?
#' @keywords internal
#' @noRd
.is_mixed <- function(object) length(unique(object$family)) > 1L

#' Format a family vector for display
#' @keywords internal
#' @noRd
.format_family <- function(family_vec) {
  if (length(unique(family_vec)) == 1L) return(unname(family_vec[1]))
  paste0("mixed (", paste(names(family_vec), family_vec,
                          sep = "=", collapse = ", "), ")")
}


## ---- resolve variable names from standata ----
#' @keywords internal
get_param_names <- function(sd) {
  p    <- sd$p
  K    <- sd$K
  n_fe <- sd$n_fe
  n_re <- sd$n_re

  y_names  <- colnames(sd$Y)  %||% paste0("y", seq_len(p))
  fe_names <- colnames(sd$X)  %||% paste0("fe", seq_len(n_fe))
  b_names  <- colnames(sd$B)  %||% paste0("lag", rep(seq_len(K), each = p), "_", rep(y_names, K))
  re_names <- if (n_re > 0) (colnames(sd$Z) %||% paste0("re", seq_len(n_re))) else character(0)

  list(y = y_names, fe = fe_names, b = b_names, re = re_names)
}


## ---- build a labelled data.frame from a draws matrix ----
## draws: iterations x parameters matrix (from extract_draws)
## row_names / col_names map to the [row, col] Stan indices
## column order in draws must follow row-major: (1,1), (2,1), ..., (nr,1), (1,2), ...
#' @keywords internal
build_summary_table <- function(draws, row_names, col_names, type) {
  nr   <- length(row_names)
  nc   <- length(col_names)
  ncol_draws <- ncol(draws)
  stopifnot(ncol_draws == nr * nc)

  d_mean   <- colMeans(draws)
  d_median <- apply(draws, 2, stats::median)
  d_q5     <- apply(draws, 2, stats::quantile, probs = 0.05)
  d_q95    <- apply(draws, 2, stats::quantile, probs = 0.95)

  data.frame(
    type      = rep(type, nr * nc),
    predictor = rep(row_names, times = nc),
    outcome   = rep(col_names, each  = nr),
    mean      = as.numeric(d_mean),
    median    = as.numeric(d_median),
    q5        = as.numeric(d_q5),
    q95       = as.numeric(d_q95),
    stringsAsFactors = FALSE
  )
}


## ---- extract posterior draws as a matrix (Stan column names preserved) ----
#' Extract raw posterior draws for a single parameter block
#'
#' Returns an \code{(iterations * chains)} by \code{params} matrix with
#' Stan-indexed column names (e.g. \code{"beta[1,1]"}, \code{"phi[2,3]"}).
#'
#' @param object A \code{bvarnet} object returned by \code{\link{bvar}}.
#' @param parameter Character. One of \code{"beta"}, \code{"phi"},
#'   \code{"sd_u"}, \code{"sigma"}, or \code{"kappa"}.
#'
#' @return A numeric matrix with one row per posterior draw and one column
#'   per Stan parameter element.
#'
#' @export
extract_draws <- function(object, parameter = c("beta", "phi", "sd_u", "sigma", "kappa")) {
  stopifnot(inherits(object, "bvarnet"))
  parameter <- match.arg(parameter, c("beta", "phi", "sd_u", "sigma", "kappa"))

  if (parameter == "sigma" && !.family_has(object, "gaussian"))
    stop("Parameter 'sigma' only exists for gaussian models.")
  if (parameter == "kappa" && !.family_has(object, "ordinal"))
    stop("Parameter 'kappa' only exists for ordinal models.")
  if (parameter == "sd_u" && object$standata$n_re == 0)
    stop("Parameter 'sd_u' not available \u2014 model has no random effects (n_re = 0).")

  draws <- object$draws                        # 3D array: iter x chains x params
  idx   <- grep(paste0("^", parameter, "\\["), dimnames(draws)[[3]])
  chunk <- draws[, , idx, drop = FALSE]
  # flatten chains into rows
  dim(chunk) <- c(prod(dim(chunk)[1:2]), dim(chunk)[3])
  colnames(chunk) <- dimnames(draws)[[3]][idx]
  chunk
}


# ---- Random-effect draw extraction (u) ----

#' Extract all posterior draws of subject-level random effects \code{u}
#'
#' Returns a 4D array with dimensions
#' \code{[draw, node, subject, re]} matching the Stan declaration
#' \code{array[p] matrix[J, n_re] u}.
#'
#' CmdStan flattens \code{array[p] matrix[J, n_re] u} as
#' \code{u[node, subject, re]} in column-major order within each array
#' element, i.e. \code{u[1,1,1], u[1,1,2], ..., u[1,J,n_re], u[2,1,1], ...}.
#'
#' @param object A \code{bvarnet} object.
#' @return A 4D array with dimensions
#'   \code{[S, p, J, n_re]} where \code{S = n_iter * n_chains}.
#' @keywords internal
#' @noRd
.extract_u_draws <- function(object) {
  stopifnot(inherits(object, "bvarnet"))
  sd <- object$standata
  if (sd$n_re == 0L)
    stop("No random effects in this model (n_re = 0).", call. = FALSE)

  p    <- sd$p
  J    <- sd$J
  n_re <- sd$n_re

  draws <- object$draws   # iter x chains x params
  idx   <- grep("^u\\[", dimnames(draws)[[3]])

  if (length(idx) == 0L)
    stop("No 'u[...]' parameters found in posterior draws.", call. = FALSE)

  # Expected number of u parameters: p * J * n_re
  expected <- p * J * n_re
  if (length(idx) != expected)
    stop(sprintf(
      "Expected %d u parameters (p=%d, J=%d, n_re=%d) but found %d in draws.",
      expected, p, J, n_re, length(idx)
    ), call. = FALSE)

  # Flatten iter x chains
  chunk <- draws[, , idx, drop = FALSE]
  S <- prod(dim(chunk)[1:2])
  dim(chunk) <- c(S, length(idx))

  # CmdStan order for array[p] matrix[J, n_re]:
  #   fastest-varying = re (col of matrix), then subject (row of matrix),
  #   then node (array index).
  # So the flat order is: u[1,1,1], u[1,1,2], ..., u[1,1,n_re],
  #                        u[1,2,1], ..., u[1,J,n_re],
  #                        u[2,1,1], ...
  # Reshape: chunk is S x (p * J * n_re).
  # We want out[draw, node, subject, re].
  out <- array(NA_real_, dim = c(S, p, J, n_re))
  col <- 0L
  for (node in seq_len(p)) {
    for (subj in seq_len(J)) {
      for (re in seq_len(n_re)) {
        col <- col + 1L
        out[, node, subj, re] <- chunk[, col]
      }
    }
  }

  dimnames(out) <- list(
    draw    = NULL,
    node    = colnames(sd$Y),
    subject = seq_len(J),
    re      = if (n_re > 0L && !is.null(colnames(sd$Z))) colnames(sd$Z) else paste0("re", seq_len(n_re))
  )

  out
}


#' Posterior mean of subject-level random effects \code{u}
#'
#' Returns a 3D array \code{[node, subject, re]} containing the
#' posterior mean of each \code{u[node, subject, re]} element.
#'
#' @param object A \code{bvarnet} object.
#' @return A 3D array with dimensions \code{[p, J, n_re]}.
#' @keywords internal
#' @noRd
.posterior_mean_u <- function(object) {
  u_draws <- .extract_u_draws(object)   # [S, p, J, n_re]
  # Average over the draw dimension (dim 1)
  apply(u_draws, c(2, 3, 4), mean)
}


## ---- text formatter for data frames (knitr-safe) ----
# Format a data frame as aligned text lines, bypassing print.data.frame
# which may be intercepted by knitr/rmarkdown in Rmd rendering.
.fmt_df <- function(df, right = FALSE) {
  m    <- format.data.frame(df, na.encode = FALSE)
  flag <- if (right) "" else "-"
  widths <- vapply(seq_along(m), function(j)
    max(nchar(names(m)[j]), max(nchar(m[[j]]))), integer(1))
  hdr <- vapply(seq_along(m), function(j)
    formatC(names(m)[j], width = widths[j], flag = flag), character(1))
  out <- character(nrow(m) + 1L)
  out[1L] <- paste0(" ", paste(hdr, collapse = " "))
  for (i in seq_len(nrow(m))) {
    vals <- vapply(seq_along(m), function(j)
      formatC(m[[j]][i], width = widths[j], flag = flag), character(1))
    out[i + 1L] <- paste0(" ", paste(vals, collapse = " "))
  }
  out
}

## ---- print method for bvarnet objects ----
#' Print a bvarnet model object
#'
#' Displays a brief summary of the fitted model: family, dimensions,
#' Rhat, divergences, chain return codes, priors, and total sampling time.
#'
#' @param x A \code{bvarnet} object.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#' @export
print.bvarnet <- function(x, ...) {
  sd <- x$standata

  cat("BVAR Network fit\n")
  cat(strrep("=", 40), "\n")

  # Family
  cat("Family:      ", .format_family(x$family), "\n", sep = "")

  # Dimensions
  cat("Outcomes (p):", sd$p,    "\n")
  cat("Lags (K):    ", sd$K,    "\n")
  if (!is.null(sd$n_fe)) cat("Fixed eff.:  ", sd$n_fe, "\n")
  if (!is.null(sd$n_re) && sd$n_re > 0)
    cat("Random eff.: ", sd$n_re, "\n")
  cat("Observations:", sd$n_obs, "\n")

  # Convergence
  smry <- x$convergence
  if (!is.null(smry) && "rhat" %in% names(smry)) {
    rhat_max <- max(smry$rhat, na.rm = TRUE)
    cat(sprintf("Rhat max:    %.3f\n", rhat_max))
    if (rhat_max > 1.01)
      cat("  WARNING: Rhat > 1.01 detected \u2014 chains may not have converged.\n")
  }

  diag <- x$diagnostics
  n_div <- if (!is.null(diag) && "num_divergent" %in% names(diag))
    sum(diag$num_divergent) else 0L
  if (n_div > 0)
    cat(sprintf("Divergences: %d  WARNING: check model/priors.\n", n_div))
  else
    cat("Divergences: 0\n")

  # Return codes
  rc <- x$return_codes
  if (!is.null(rc) && any(rc != 0))
    cat(sprintf("Chain status: %d chain(s) returned non-zero exit codes: %s\n",
                sum(rc != 0), paste(which(rc != 0), collapse = ", ")))

  # Priors
  if (!is.null(x$priors) && inherits(x$priors, "bvarnet_priors")) {
    half_pars <- c("sd_u", "sigma")

    # Determine which priors to show
    priors_to_show <- x$priors_needed
    if (is.null(priors_to_show)) {
      # Backward compat: old objects without priors_needed
      has_intercept <- .family_has(x, "gaussian") || .family_has(x, "bernoulli")
      priors_to_show <- if (has_intercept) c("intercept", "beta", "phi") else c("beta", "phi")
      has_re <- !is.null(sd$n_re) && sd$n_re > 0
      if (has_re) priors_to_show <- c(priors_to_show, "sd_u")
      if (.family_has(x, "gaussian")) priors_to_show <- c(priors_to_show, "sigma")
      if (.family_has(x, "ordinal"))  priors_to_show <- c(priors_to_show, "kappa")
    }
    # Filter to slots that actually exist (robustness for serialized objects)
    priors_to_show <- Filter(function(nm) !is.null(x$priors[[nm]]), priors_to_show)

    any_user <- any(!vapply(
      x$priors[priors_to_show],
      function(p) isTRUE(p$is_default), logical(1)
    ))

    if (any_user) {
      # Show each prior on its own line with (default) tags
      cat("Priors:\n")
      for (nm in priors_to_show) {
        half <- nm %in% half_pars
        tag <- if (isTRUE(x$priors[[nm]]$is_default)) "  (default)" else ""
        cat(sprintf("  %-6s ~ %s%s\n", nm, format(x$priors[[nm]], half = half), tag))
      }
    } else {
      # All defaults — compact single line
      prior_str <- paste(
        vapply(priors_to_show, function(nm) {
          half <- nm %in% half_pars
          paste0(nm, " ~ ", format(x$priors[[nm]], half = half))
        }, character(1L)),
        collapse = ", "
      )
      cat("Priors:       ", prior_str, " (all defaults)\n", sep = "")
    }
  }

  # Timing
  t <- x$timing
  if (!is.null(t$total))
    cat(sprintf("Total time:  %.1f sec\n", t$total))

  cat(strrep("=", 40), "\n")
  invisible(x)
}
