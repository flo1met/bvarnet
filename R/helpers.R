## ---- sampler argument passthrough ----

#' Sampler arguments `bvar()` supplies itself, mapped to the `bvar()` argument
#' that controls each one.
#'
#' The first block must stay in sync with the `$sample()` calls in `bvar()` and
#' `.bvar_nodewise()`: every name spliced in there belongs here, otherwise a
#' user could supply it via `...` and reach `do.call()` with a duplicated
#' argument name.
#'
#' The second block holds CmdStanR's deprecated or removed aliases for those
#' same arguments. In CmdStanR versions that still provide them, `$sample()`
#' resolves an alias by overwriting its modern equivalent (e.g.
#' `iter_warmup <- num_warmup`), so letting one through `...` would silently
#' override the value `bvar()` set behind a deprecation warning. Keeping the
#' removed aliases reserved also gives users of newer CmdStanR versions the
#' same informative `bvar()` error. `stepsize`, `validate_csv` and
#' `save_extra_diagnostics` were aliases for arguments `bvar()` does not set,
#' so they are not reserved here.
#' @keywords internal
#' @noRd
.bvarnet_reserved_sampler_args <- c(
  data            = "data",
  seed            = "seed",
  iter_warmup     = "warmup",
  iter_sampling   = "iter",
  chains          = "chains",
  parallel_chains = "cores",
  adapt_delta     = "adapt_delta",
  max_treedepth   = "max_treedepth",

  num_warmup      = "warmup",
  num_samples     = "iter",
  num_chains      = "chains",
  num_cores       = "cores",
  cores           = "cores",
  max_depth       = "max_treedepth"
)

#' The `$sample()` arguments `bvar()` always supplies by exact name.
#'
#' Used to figure out which formals remain available for R's partial
#' argument matching once `bvar()`'s own exactly-named arguments have
#' claimed their formals (see `.bvarnet_sample_partial_pool()`).
#' @keywords internal
#' @noRd
.bvarnet_fixed_sample_args <- c(
  "data", "seed", "iter_warmup", "iter_sampling",
  "chains", "parallel_chains", "adapt_delta", "max_treedepth"
)

#' `$sample()` formals still available for partial matching against `...`.
#'
#' `$sample()` has no `...` of its own, so `do.call()` partial-matches any
#' dot name that isn't an exact formal name against the formals `bvar()`
#' hasn't already claimed by exact name. A dot name like `num_warm` then
#' silently resolves to the deprecated `num_warmup` formal and overrides the
#' `warmup` value `bvar()` set. Returns `character(0)` (disabling the
#' partial-match check) if CmdStanR isn't available or its `$sample()`
#' signature can't be introspected; `.check_sampler_dots()` falls back to
#' exact-name matching in that case.
#' @keywords internal
#' @noRd
.bvarnet_sample_partial_pool <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(character(0))
  sample_method <- tryCatch(
    utils::getFromNamespace("CmdStanModel", "cmdstanr")$public_methods$sample,
    error = function(e) NULL
  )
  if (is.null(sample_method)) return(character(0))
  setdiff(names(formals(sample_method)), .bvarnet_fixed_sample_args)
}

#' Validate the `...` arguments `bvar()` forwards to CmdStanR's `$sample()`.
#'
#' Returns the dots as a named list, ready to splice into a `do.call()`.
#' @keywords internal
#' @noRd
.check_sampler_dots <- function(...) {
  dots <- list(...)
  if (length(dots) == 0L) return(list())

  nms <- names(dots)
  if (is.null(nms) || !all(nzchar(nms)))
    stop("All arguments passed to bvar() via '...' must be named. ",
         "They are forwarded to the CmdStanR $sample() method.", call. = FALSE)

  reserved <- .bvarnet_reserved_sampler_args
  pool <- .bvarnet_sample_partial_pool()

  # Resolve each dot name to the reserved name it would collide with at
  # do.call() time: an exact match, or (via R's partial-matching rules) a
  # name that uniquely abbreviates a still-available reserved alias.
  resolved <- vapply(nms, function(nm) {
    if (nm %in% names(reserved)) return(nm)
    if (!length(pool)) return(NA_character_)
    hit <- pmatch(nm, pool)
    if (is.na(hit)) return(NA_character_)
    target <- pool[hit]
    if (target %in% names(reserved)) target else NA_character_
  }, character(1), USE.NAMES = FALSE)

  clash <- !is.na(resolved)
  if (any(clash))
    stop("These $sample() arguments are set by bvar() and cannot be passed ",
         "via '...':\n",
         paste0("  '", nms[clash], "'",
                ifelse(nms[clash] == resolved[clash], "",
                       paste0(" (matches '", resolved[clash], "')")),
                " - use bvar(", reserved[resolved[clash]], " = ...)",
                collapse = "\n"),
         call. = FALSE)

  dup <- unique(nms[duplicated(nms)])
  if (length(dup))
    stop("Duplicated argument(s) passed to bvar() via '...': ",
         paste0("'", dup, "'", collapse = ", "), call. = FALSE)

  dots
}

#' Per-node `dots` for `.bvar_nodewise()`'s `$sample()` calls.
#'
#' A user-supplied `output_basename` disables CmdStanR's timestamp/random CSV
#' suffixes (see `?cmdstanr::"model-method-sample"`), so every per-node fit
#' would otherwise write to the same output files and overwrite each other's
#' draws before `.combine_nodewise_fits()` reads them. Suffix it with the
#' node index to keep each node's files distinct.
#' @keywords internal
#' @noRd
.bvarnet_node_dots <- function(dots, j) {
  if (is.null(dots$output_basename)) return(dots)
  dots$output_basename <- paste0(dots$output_basename, "_node", j)
  dots
}


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


## ---- credible-interval level -> quantile probabilities ----
## Validates ci_level and returns the equal-tailed lower/upper probabilities.
#' @keywords internal
#' @noRd
.ci_probs <- function(ci_level, arg_name = "ci_level") {
  if (!is.numeric(ci_level) || length(ci_level) != 1L || !is.finite(ci_level) ||
      ci_level <= 0 || ci_level >= 1)
    stop(sprintf(
      "`%s` must be a single number strictly between 0 and 1 (e.g. 0.95 for a 95%% credible interval).",
      arg_name
    ), call. = FALSE)
  alpha <- (1 - ci_level) / 2
  c(alpha, 1 - alpha)
}


## ---- core mean/median/CI computation from a draws matrix ----
## Shared by build_summary_table() and extract_param()'s manual beta/kappa
## blocks. The two CI bounds share one probs vector so they're computed in a
## single quantile() pass; median is kept as a separate stats::median() call
## since quantile(probs=0.5) is not bit-identical to median() for even n.
#' @keywords internal
#' @noRd
.summarize_draws <- function(draws, probs) {
  q <- apply(draws, 2L, stats::quantile, probs = c(probs[1L], probs[2L]))
  list(
    mean     = as.numeric(colMeans(draws)),
    median   = as.numeric(apply(draws, 2L, stats::median)),
    ci_lower = as.numeric(q[1L, ]),
    ci_upper = as.numeric(q[2L, ])
  )
}


## ---- join convergence diagnostics (rhat/ess_bulk/ess_tail) onto a table ----
#' @keywords internal
#' @noRd
.join_convergence <- function(tab, stan_colnames, convergence) {
  idx <- match(stan_colnames, convergence$variable)
  tab$rhat     <- convergence$rhat[idx]
  tab$ess_bulk <- convergence$ess_bulk[idx]
  tab$ess_tail <- convergence$ess_tail[idx]
  tab
}


## ---- build a labelled data.frame from a draws matrix ----
## draws: iterations x parameters matrix (from extract_draws)
## row_names / col_names map to the [row, col] Stan indices
## column order in draws must follow row-major: (1,1), (2,1), ..., (nr,1), (1,2), ...
#' @keywords internal
build_summary_table <- function(draws, row_names, col_names, type,
                                ci_level = 0.95, probs = NULL) {
  nr   <- length(row_names)
  nc   <- length(col_names)
  ncol_draws <- ncol(draws)
  stopifnot(ncol_draws == nr * nc)

  if (is.null(probs)) probs <- .ci_probs(ci_level)
  s <- .summarize_draws(draws, probs)

  data.frame(
    type      = rep(type, nr * nc),
    predictor = rep(row_names, times = nc),
    outcome   = rep(col_names, each  = nr),
    mean      = s$mean,
    median    = s$median,
    ci_lower  = s$ci_lower,
    ci_upper  = s$ci_upper,
    stringsAsFactors = FALSE
  )
}


## ---- random-effect SD summary table (sd_u) ----
## Labels rows by parsing the [node, re] indices out of the Stan names rather
## than by column position. sd_u is declared matrix[p, n_re], so the joint path
## emits it node-fastest while .combine_nodewise_fits() binds node-slowest;
## only a name-based mapping is correct for both.
#' @keywords internal
#' @noRd
.sd_u_summary_table <- function(object, probs) {
  nm       <- get_param_names(object$standata)
  draws_sd <- extract_draws(object, "sd_u")
  cn       <- colnames(draws_sd)

  ix <- .stan_indices(cn)   # columns: node, re
  if (ncol(ix) != 2L)
    stop("Expected two indices in 'sd_u[...]' parameter names.", call. = FALSE)

  s   <- .summarize_draws(draws_sd, probs)
  tab <- data.frame(
    type      = "Random Effect SD",
    predictor = nm$y[ix[, 1L]],
    outcome   = nm$re[ix[, 2L]],
    mean      = s$mean,
    median    = s$median,
    ci_lower  = s$ci_lower,
    ci_upper  = s$ci_upper,
    stringsAsFactors = FALSE
  )
  tab <- .join_convergence(tab, cn, object$convergence)
  rownames(tab) <- NULL
  tab
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
#' Columns are located by parsing the \code{[node, subject, re]} indices out of
#' the Stan parameter names, never by position: CmdStan flattens every index
#' column-major with the first index fastest (\code{u[1,1,1], u[2,1,1],
#' u[1,2,1], ...}), while \code{.combine_nodewise_fits()} binds node-outer.
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

  ix <- .stan_indices(dimnames(draws)[[3]][idx])   # columns: node, subject, re
  if (ncol(ix) != 3L)
    stop("Expected three indices in 'u[...]' parameter names.", call. = FALSE)
  if (any(ix[, 1L] > p) || any(ix[, 2L] > J) || any(ix[, 3L] > n_re) || any(ix < 1L))
    stop("A 'u[...]' index is out of range for (p, J, n_re).", call. = FALSE)

  out <- array(NA_real_, dim = c(S, p, J, n_re))
  for (k in seq_along(idx))
    out[, ix[k, 1L], ix[k, 2L], ix[k, 3L]] <- chunk[, k]

  # A duplicate or missing index would leave a cell unfilled.
  if (anyNA(out))
    stop("Posterior draws do not cover every u[node, subject, re] element.",
         call. = FALSE)

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


## ---- disclose data-dependent prior scaling ----
## Default Gaussian priors are widened by the outcome SD before they reach Stan
## (see .scale_default_priors()). Printing only the requested scale would hide
## that, so report the effective one whenever it differs.
#' @keywords internal
#' @noRd
.effective_scale_note <- function(x, nm) {
  pe <- x$priors_effective
  if (is.null(pe) || is.null(x$priors[[nm]])) return("")

  base <- x$priors[[nm]]$scale
  eff  <- vapply(pe, function(p) p[[nm]]$scale %||% NA_real_, numeric(1L))
  eff  <- unique(eff[!is.na(eff)])
  if (length(eff) == 0L || all(abs(eff - base) < 1e-8)) return("")

  if (length(eff) == 1L)
    sprintf("  [scaled by sd(y) to %g]", eff)
  else
    sprintf("  [scaled by sd(y) to %g-%g across nodes]", min(eff), max(eff))
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
    notes <- vapply(priors_to_show, .effective_scale_note, character(1L), x = x)

    if (any_user || any(nzchar(notes))) {
      # Show each prior on its own line with (default) tags
      cat("Priors:\n")
      for (nm in priors_to_show) {
        half <- nm %in% half_pars
        tag <- if (isTRUE(x$priors[[nm]]$is_default)) "  (default)" else ""
        cat(sprintf("  %-9s ~ %s%s%s\n", nm, format(x$priors[[nm]], half = half),
                    tag, notes[[nm]]))
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
