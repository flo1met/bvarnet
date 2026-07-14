# -----------------------------------------------------------------------------
# Small internal utilities shared across the package.
# -----------------------------------------------------------------------------

#' NULL-coalescing operator
#'
#' Base R gained `%||%` only in 4.4.0; DESCRIPTION declares R (>= 3.5), so the
#' package defines its own.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x


#' Parse the bracket indices out of Stan parameter names
#'
#' Turns names like \code{"u[2,3,1]"} into a matrix of integer indices, one row
#' per name and one column per index position. All names must carry the same
#' number of indices.
#'
#' Every consumer of a draws array must resolve parameters by parsed index
#' rather than by column position: CmdStan flattens column-major (first index
#' fastest), while the nodewise combiner in \code{.combine_nodewise_fits()}
#' emits some blocks node-slowest.
#'
#' @param names Character vector of Stan parameter names.
#' @return An integer matrix with \code{length(names)} rows.
#' @keywords internal
#' @noRd
.stan_indices <- function(names) {
  m <- regmatches(names, regexec("^[A-Za-z_][A-Za-z0-9_]*\\[([0-9,]+)\\]$", names))
  bad <- lengths(m) != 2L
  if (any(bad))
    stop(sprintf("Cannot parse Stan parameter index from: %s",
                 paste(names[bad], collapse = ", ")), call. = FALSE)

  parts <- strsplit(vapply(m, `[[`, character(1L), 2L), ",", fixed = TRUE)
  n_idx <- lengths(parts)
  if (length(unique(n_idx)) != 1L)
    stop("Stan parameter names have inconsistent index counts.", call. = FALSE)

  matrix(as.integer(unlist(parts, use.names = FALSE)),
         nrow = length(names), ncol = n_idx[1L], byrow = TRUE)
}
