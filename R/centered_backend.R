# Internal predictor-centering coordinates. R's public X/B/Z and prior reference
# are unchanged; only Stan's sampled intercepts/thresholds use centered columns.

.bvarnet_public_draws <- function(draws) {
  public <- dimnames(draws)[[3L]]
  public <- public[!grepl("^(beta_centered|kappa_centered)\\[", public)]
  posterior::subset_draws(draws, variable = public)
}

# Complete public location/slope initializers are required: otherwise an
# unspecified slope would make the location shift depend on a random start
# generated later inside CmdStan. Other sampled blocks may remain unspecified.
.bvarnet_centered_init_one <- function(init, standata, family) {
  fail <- function(message) stop("Invalid 'init': ", message, call. = FALSE)
  if (!is.list(init) || is.data.frame(init) || !length(init) ||
      is.null(names(init)) || anyNA(names(init)) ||
      any(!nzchar(names(init))) || anyDuplicated(names(init)))
    fail("each chain must supply a uniquely named parameter list.")

  if (any(c("beta_centered", "kappa_centered") %in% names(init)))
    fail("use public beta/phi/kappa, not internal centered parameter names.")
  allowed <- c("beta", "phi", "z_u", "sd_u",
               if (family == "gaussian") "sigma",
               if (family == "ordinal") "kappa")
  unknown <- setdiff(names(init), allowed)
  if (length(unknown))
    fail(paste0("unsupported parameter(s): ", paste(unknown, collapse = ", "), "."))

  check_matrix <- function(name, dimensions) {
    x <- init[[name]]
    if (length(x) == 0L && prod(dimensions) == 0L)
      return(matrix(numeric(0), dimensions[1L], dimensions[2L]))
    # An ordinal array of vectors may be supplied as a list of node vectors.
    if (name == "kappa" && is.list(x)) x <- do.call(rbind, x)
    if (!is.numeric(x) || !identical(as.integer(dim(x)), as.integer(dimensions)) ||
        any(!is.finite(x)))
      fail(paste0("supply complete public beta/phi",
                  if (family == "ordinal") "/kappa", "; '", name,
                  "' must be a finite matrix with dimensions ",
                  paste(dimensions, collapse = " x "), "."))
    x
  }

  beta <- check_matrix("beta", c(standata$n_fe, standata$p))
  phi <- check_matrix("phi", c(standata$p * standata$K, standata$p))
  x_means <- colMeans(standata$X)
  if (family != "ordinal") x_means[1L] <- 0
  shift <- as.vector(crossprod(x_means, beta) +
                       crossprod(colMeans(standata$B), phi))
  if (any(!is.finite(shift)))
    fail("the location shift is nonfinite; rescale the predictors or initial values.")

  if (family == "ordinal") {
    kappa <- check_matrix("kappa", c(standata$p, standata$C - 1L))
    # JSON cannot encode the second extent of a [0,p] matrix. Let Stan
    # initialize this empty block rather than emitting an ambiguous [].
    init$beta <- if (standata$n_fe > 0L) beta else NULL
    init$kappa <- NULL
    init$kappa_centered <- sweep(kappa, 1L, shift, "-")
  } else {
    beta[1L, ] <- beta[1L, ] + shift
    init$beta <- NULL
    init$beta_centered <- beta
  }
  location <- if (family == "ordinal") init$kappa_centered else init$beta_centered
  if (any(!is.finite(location)))
    fail("the centered location is nonfinite; rescale the predictors or initial values.")
  init
}

.bvarnet_centered_init <- function(init, standata, family, chains, node_id = NULL) {
  if (is.null(init)) return(NULL)
  if (is.numeric(init) && is.null(dim(init)) && length(init) == 1L &&
      is.finite(init) && init >= 0) return(init)

  if (is.function(init)) {
    args <- names(formals(init))
    if (any(!args %in% c("chain_id", "node_id")))
      stop("An 'init' function may take only chain_id and (for mixed fits) node_id.",
           call. = FALSE)
    if ("node_id" %in% args && is.null(node_id))
      stop("The 'init' function's node_id argument is only available for mixed fits.",
           call. = FALSE)
    init <- lapply(seq_len(chains), function(chain_id) {
      supplied <- list(chain_id = chain_id, node_id = node_id)
      do.call(init, supplied[args])
    })
  } else if (is.character(init)) {
    if (!length(init) || anyNA(init) ||
        !length(init) %in% c(1L, chains) || any(!file.exists(init)))
      stop("'init' must name one existing JSON file or one file per chain.", call. = FALSE)
    init <- lapply(rep(init, length.out = chains), function(path)
      jsonlite::fromJSON(path, simplifyVector = TRUE))
  } else if (is.list(init) && !is.object(init)) {
    # Accept one named parameter list as well as CmdStanR's per-chain lists.
    if (length(init) && !is.null(names(init)) &&
        any(names(init) %in% c("beta", "phi", "kappa", "sigma", "sd_u", "z_u",
                               "beta_centered", "kappa_centered")))
      init <- rep(list(init), chains)
  } else {
    stop("Unsupported 'init'. Use a numeric radius, public parameter lists, ",
         "a function, or JSON files. Fit/draw-based warm starts are not supported ",
         "by the centered backend; convert them to complete public beta/phi/kappa lists.",
         call. = FALSE)
  }
  if (length(init) != chains)
    stop("'init' must provide one parameter list per chain.", call. = FALSE)
  lapply(init, .bvarnet_centered_init_one, standata = standata, family = family)
}
