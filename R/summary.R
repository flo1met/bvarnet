#' Summary method for bvarnet objects
#'
#' Returns a labelled posterior summary table grouped by parameter type,
#' with convergence diagnostics and optional Bayes factors. Wraps
#' \code{\link{extract_param}}.
#'
#' @param object A \code{bvarnet} object returned by \code{\link{bvar}}.
#' @param bayes_factor Logical; if \code{TRUE}, append Savage-Dickey BF
#'   columns. Default \code{FALSE}.
#' @param null_value Numeric scalar; null hypothesis value for BF computation.
#'   Default 0.
#' @param ... Ignored.
#'
#' @return An object of class \code{"summary.bvarnet"} (a list) with elements:
#'   \describe{
#'     \item{table}{Data frame from \code{extract_param()}.}
#'     \item{family}{Model family.}
#'     \item{p}{Number of outcome variables.}
#'     \item{K}{AR order.}
#'     \item{n}{Number of observations.}
#'     \item{rhat_max}{Maximum Rhat across all parameters.}
#'     \item{n_divergences}{Total divergent transitions.}
#'   }
#'
#' @export
summary.bvarnet <- function(object, bayes_factor = FALSE, null_value = 0, ...) {
  stopifnot(inherits(object, "bvarnet"))

  tab <- extract_param(object,
                        bayes_factor = bayes_factor,
                        null_value   = null_value)

  sd   <- object$standata
  conv <- object$convergence
  diag <- object$diagnostics

  rhat_max <- if (!is.null(conv) && "rhat" %in% names(conv))
    max(conv$rhat, na.rm = TRUE) else NA_real_
  n_div <- if (!is.null(diag) && "num_divergent" %in% names(diag))
    sum(diag$num_divergent) else 0L

  structure(
    list(
      table         = tab,
      family        = object$family,
      p             = sd$p,
      K             = sd$K,
      n             = sd$n,
      rhat_max      = rhat_max,
      n_divergences = n_div
    ),
    class = "summary.bvarnet"
  )
}


#' Print a bvarnet summary
#'
#' Pretty-prints the output of \code{\link{summary.bvarnet}}, grouping
#' parameters by type and displaying convergence information.
#' Each group is truncated to \code{max_rows} rows; use \code{extract_param()}
#' or dedicated extractors to see full output.
#'
#' @param x A \code{summary.bvarnet} object.
#' @param digits Number of decimal digits for numeric columns. Default 3.
#' @param max_rows Maximum number of rows to print per parameter group.
#'   Default 10.
#' @param ... Ignored.
#'
#' @return \code{x} invisibly.
#' @export
print.summary.bvarnet <- function(x, digits = 3, max_rows = 10, ...) {
  cat("BVAR Network Summary\n")
  cat(strrep("=", 50), "\n")
  cat(sprintf("Family: %s | p=%d | K=%d | n=%d\n",
              .format_family(x$family), x$p, x$K, x$n))
  cat(sprintf("Rhat max: %.3f | Divergences: %d\n",
              x$rhat_max, x$n_divergences))

  if (x$rhat_max > 1.01)
    cat("  WARNING: Rhat > 1.01 \u2014 chains may not have converged.\n")
  if (x$n_divergences > 0)
    cat("  WARNING: divergent transitions detected \u2014 check model/priors.\n")

  tab <- x$table

  # Determine which columns to print
  display_cols <- c("predictor", "outcome", "mean", "median",
                    "q5", "q95", "rhat", "ess_bulk", "ess_tail")
  if ("BF10" %in% names(tab))
    display_cols <- c(display_cols, "BF10")
  display_cols <- intersect(display_cols, names(tab))

  # Group by type and print each section
  types <- unique(tab$type)
  type_order <- c("Intercept", "Fixed Effect", "Autoregressive",
                  "Cross-lagged", "Random Effect SD", "Residual SD",
                  "Threshold")
  types <- intersect(type_order, types)

  for (tp in types) {
    cat(sprintf("\n--- %s ---\n", tp))
    sub <- tab[tab$type == tp, display_cols, drop = FALSE]
    rownames(sub) <- NULL
    n_total <- nrow(sub)
    truncated <- n_total > max_rows

    if (truncated) sub <- sub[seq_len(max_rows), , drop = FALSE]

    # Round numeric columns
    num_cols <- vapply(sub, is.numeric, logical(1))
    sub[num_cols] <- lapply(sub[num_cols], round, digits = digits)

    print(sub, row.names = FALSE, right = FALSE)

    if (truncated)
      cat(sprintf("... %d more rows. Use extract_param() for full output.\n",
                  n_total - max_rows))
  }

  cat("\n", strrep("=", 50), "\n", sep = "")
  cat("Use extract_param() for the full parameter table.\n")
  cat("Use extract_temporal() for autoregressive / cross-lagged effects.\n")
  cat("Use extract_network_matrix() for p x p network matrices.\n")
  invisible(x)
}
