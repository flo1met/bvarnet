#' Apply the data-dependent default-prior scaling
#'
#' Gaussian outcomes are modelled on their raw scale, so the unit-scale default
#' priors on \code{beta}, \code{intercept} and \code{sigma} are widened by the
#' outcome SD. User-supplied priors are taken at face value and never rescaled.
#'
#' This is the single definition of the scaling rule. Both \code{to_stan_data()}
#' and \code{.to_stan_data_node()} call it, and the result is the *effective*
#' prior that Stan actually samples under — the one \code{savage_dickey()} must
#' use as the Savage-Dickey denominator.
#'
#' @param priors A \code{bvarnet_priors} object.
#' @param family Character scalar.
#' @param s_y Numeric scalar. The outcome SD to scale by.
#'
#' @return A \code{bvarnet_priors} object, scaled where applicable.
#' @keywords internal
#' @noRd
.scale_default_priors <- function(priors, family, s_y) {
  if (family != "gaussian" || !is.finite(s_y) || s_y <= 0) return(priors)

  for (par in c("beta", "intercept", "sigma"))
    if (isTRUE(priors[[par]]$is_default))
      priors[[par]]$scale <- priors[[par]]$scale * s_y

  priors
}


#' Build a Stan data list from a long-format data frame
#'
#' Internal function called by \code{bvar()}. Constructs the list passed to
#' the Stan model: design matrices Y, X, B, Z; dimensions; and prior
#' hyperparameters derived from a \code{bvarnet_priors} object.
#'
#' @param data Data frame in long format.
#' @param family Character. \code{"bernoulli"}, \code{"ordinal"}, or
#'   \code{"gaussian"}.
#' @param id_col Character. Subject/group identifier column name.
#' @param time_col Character. Scheduled occasion column when `beep_col` is NULL;
#'   otherwise the day/session label column. Lags never cross day/session labels.
#' @param beep_col Optional scheduled occasion column within each day/session.
#'   With NULL, `time_col` supplies the integer occasion index. With both supplied,
#'   only consecutive integer beeps within the same subject and segment are used.
#' @param y_cols Character vector. Outcome column names.
#' @param x_cols Character vector. Covariate column names.
#' @param center_x Logical. Grand-mean centre X. Default \code{FALSE}.
#' @param fe_interactions List or NULL. Fixed-effect interaction terms.
#' @param re_interactions List or NULL. Random-effect interaction terms.
#' @param re_cols Character vector. Random-slope columns.
#' @param re_temporal Logical. Random slopes on lag predictors. Default FALSE.
#' @param K Integer. AR order.
#' @param na_action Character. Only \code{"listwise"} currently supported.
#' @param priors A \code{bvarnet_priors} object. Defaults to \code{set_priors()}.
#'
#' @return A named list ready to pass to \code{CmdStanModel$sample()}.
#' @keywords internal
#' @noRd
to_stan_data <- function(data,
                         family,
                         id_col,
                         time_col,
                         y_cols,
                         x_cols,
                         center_x = F, # grand-mean centering of covmat X
                         fe_interactions = NULL,
                         re_interactions = NULL,
                         re_cols = character(0),
                         re_temporal = FALSE,
                         K,
                         na_action = c("listwise"), # add automatic LW deletion
                         priors = set_priors(),
                         beep_col = NULL) {

  ## Use shared helper for all family-agnostic data wrangling
  family <- match.arg(family, choices = c("bernoulli", "ordinal", "gaussian"))

  shared <- .to_stan_data_shared(
    data = data, id_col = id_col, time_col = time_col,
    y_cols = y_cols, x_cols = x_cols, center_x = center_x,
    fe_interactions = fe_interactions, re_interactions = re_interactions,
    re_cols = re_cols, re_temporal = re_temporal, K = K,
    na_action = na_action, beep_col = beep_col
  )

  ## Family-specific Y casting
  Y <- shared$Y
  X <- shared$X
  if (family == "ordinal") {
    if (any(Y != floor(Y), na.rm = TRUE)) {
      stop("Ordinal Y must contain integer values. Found non-integer entries.")
    }
    Y <- matrix(as.integer(Y), nrow = nrow(Y), ncol = ncol(Y))
    colnames(Y) <- colnames(shared$Y)
    # Strip intercept from X for ordinal
    if ("Intercept" %in% colnames(X)) {
      X <- X[, colnames(X) != "Intercept", drop = FALSE]
    }
  } else if (family == "bernoulli") {
    Y <- matrix(as.integer(Y), nrow = nrow(Y), ncol = ncol(Y))
    colnames(Y) <- colnames(shared$Y)
  }

  n_fe <- ncol(X)

  if (family == "ordinal") {
    C <- max(Y, na.rm = TRUE)
    if (min(Y, na.rm = TRUE) < 1L) {
      stop("Ordinal Y must be coded as integers starting at 1. Found values < 1.")
    }
  }

  ## Data-dependent scaling for Gaussian. The joint path shares one s_y across
  ## all outcomes; .to_stan_data_node() scales per node.
  priors <- .scale_default_priors(priors, family,
                                  s_y = mean(apply(Y, 2, sd, na.rm = TRUE)))

  ## Rebuild Z from the (possibly stripped) X
  Z <- build_Z(X, shared$B, re_cols = re_cols, re_temporal = re_temporal)
  Z <- add_re_interactions_from_X(Z, X, shared$B, re_interactions)

  out <- list(
    p = shared$p,
    J = shared$J,
    K = shared$K,
    n_obs = shared$n_obs,
    n_fe = n_fe,
    n_re = ncol(Z),
    id = shared$id,
    Y = Y,
    X = X,
    B = shared$B,
    Z = Z,
    prior_beta_fam = priors$beta$family_int,
    beta_loc       = priors$beta$loc,
    beta_scale     = priors$beta$scale,
    beta_df        = priors$beta$df,
    prior_phi_fam  = priors$phi$family_int,
    phi_loc        = priors$phi$loc,
    phi_scale      = priors$phi$scale,
    phi_df         = priors$phi$df,
    prior_sd_fam   = priors$sd_u$family_int,
    sd_loc         = priors$sd_u$loc,
    sd_scale       = priors$sd_u$scale,
    sd_df          = priors$sd_u$df
  )

  ## attach FE interaction metadata (always set, even when empty)
  out$fe_interaction_terms <- shared$fe_interaction_terms
  if (length(shared$fe_interaction_terms) > 0L) {
    out$fe_interaction_colnames <- shared$fe_interaction_colnames
  }

  ## prediction metadata (not passed to Stan)
  out$id_levels <- shared$id_levels
  out$time_obs <- shared$time_obs
  out$index_obs <- shared$index_obs
  out$x_center_means <- shared$x_center_means
  out$row_map <- shared$row_map
  out$n_rows_data <- shared$n_rows_data
  out$design_spec <- shared$design_spec

  if (family == "ordinal") {
    out$C              <- C
    out$prior_kappa_fam <- priors$kappa$family_int
    out$kappa_loc      <- priors$kappa$loc
    out$kappa_scale    <- priors$kappa$scale
    out$kappa_df       <- priors$kappa$df
  }

  if (family == "gaussian") {
    out$prior_sigma_fam <- priors$sigma$family_int
    out$sigma_loc       <- priors$sigma$loc
    out$sigma_scale     <- priors$sigma$scale
    out$sigma_df        <- priors$sigma$df
  }

  if (family %in% c("gaussian", "bernoulli")) {
    out$prior_intercept_fam <- priors$intercept$family_int
    out$intercept_loc       <- priors$intercept$loc
    out$intercept_scale     <- priors$intercept$scale
    out$intercept_df        <- priors$intercept$df
  }

  # An attribute, not a list element: bvar() filters standata by name before
  # $sample(), but .bvar_nodewise() passes its node data list through
  # unfiltered, and cmdstanr serialises names() only.
  attr(out, "priors_effective") <- priors

  return(out)
}


# ---- .to_stan_data_shared() — family-agnostic design matrices ----

#' Build shared design matrices (family-agnostic)
#'
#' Always includes the Intercept column in X. Does NOT pack priors or do
#' family-specific type casting. Used by both \code{to_stan_data()} (joint
#' path) and \code{.bvar_nodewise()} (mixed path).
#'
#' The occasion index (time_col alone, or beep_col within day/session labels)
#' must be integer-valued. Every required predecessor must exist in the same
#' subject and segment. Incomplete windows are always excluded as responses.
#'
#' @return A list with p, J, K, n_obs, n_fe, n_re, id, Y, X, B, Z,
#'   id_levels, x_center_means, row_map, n_rows_data, design_spec,
#'   fe_interaction_terms, fe_interaction_colnames.
#' @keywords internal
#' @noRd
.to_stan_data_shared <- function(data, id_col, time_col, y_cols, x_cols,
                                  center_x = FALSE, fe_interactions = NULL,
                                  re_interactions = NULL, re_cols = character(0),
                                  re_temporal = FALSE, K, na_action = "listwise",
                                  beep_col = NULL) {

  if (!is.numeric(K) || length(K) != 1L || !is.finite(K) ||
      K < 1 || K != floor(K) || K > .Machine$integer.max)
    stop("K must be a positive integer.", call. = FALSE)
  K <- as.integer(K)
  if (!is.data.frame(data)) stop("data must be a data frame.", call. = FALSE)
  keys <- c(id_col, time_col, beep_col)
  for (name in c("id_col", "time_col", "beep_col")) {
    value <- get(name)
    if (name == "beep_col" && is.null(value)) next
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value))
      stop(name, " must name one column in data.", call. = FALSE)
  }
  if (anyDuplicated(keys)) stop("Subject, time and beep columns must be distinct.", call. = FALSE)
  required <- unique(c(keys, y_cols, x_cols))
  absent <- setdiff(required, names(data))
  if (length(absent)) stop("Columns not found in data: ", paste(absent, collapse = ", "), call. = FALSE)
  if (!length(y_cols)) stop("y_cols must select at least one outcome.", call. = FALSE)
  for (key in keys) {
    value <- data[[key]]
    if (!(is.numeric(value) || is.character(value) || is.factor(value) || inherits(value, "Date")) ||
        !is.null(dim(value)) || anyNA(value) ||
        ((is.numeric(value) || inherits(value, "Date")) && any(!is.finite(value))))
      stop("Index column '", key, "' must contain nonmissing finite numeric, Date, or character/factor labels.", call. = FALSE)
  }
  occasion_col <- if (is.null(beep_col)) time_col else beep_col
  tt <- data[[occasion_col]]
  if (!is.numeric(tt) || is.object(tt) || any(abs(tt - round(tt)) > 1e-6))
    stop("`", occasion_col, "` must be integer-valued: one time unit = one lag step. ",
         "Rescale to the scheduled grid and preserve missing occasion numbers.", call. = FALSE)
  data[[occasion_col]] <- round(tt)
  data <- data[do.call(order, unname(data[keys])), , drop = FALSE]
  n_input <- nrow(data)

  na_action <- match.arg(na_action, "listwise")
  if (na_action == "listwise") {
    check_cols <- c(y_cols, x_cols)
    complete <- complete.cases(data[, check_cols, drop = FALSE])
    data <- data[complete, , drop = FALSE]
  }

  # Validate retained keys, preserving the listwise-deletion duplicate contract.
  dup <- duplicated(data[keys]) | duplicated(data[keys], fromLast = TRUE)
  if (any(dup)) {
    stop("Duplicate (", paste(keys, collapse = ", "), ") rows found. ",
         "Each retained subject/occasion key must be unique; deduplicate or aggregate before fitting.", call. = FALSE)
  }

  ids_unique <- unique(data[[id_col]])
  J <- length(ids_unique)
  p <- length(y_cols)
  q <- length(x_cols)
  PK <- p * K

  n_rows_data <- nrow(data)
  id_all <- match(data[[id_col]], ids_unique)
  segment <- if (is.null(beep_col)) rep(1L, n_rows_data) else
    match(data[[time_col]], unique(data[[time_col]]))
  groups <- split(seq_len(n_rows_data), interaction(id_all, segment, drop = TRUE))
  predecessors <- matrix(NA_integer_, n_rows_data, K)
  occasions <- data[[occasion_col]]
  for (group in groups) {
    for (lag in seq_len(K)) {
      predecessors[group, lag] <- group[match(occasions[group] - lag, occasions[group])]
    }
  }
  row_map <- which(rowSums(is.na(predecessors)) == 0L)
  n_obs <- length(row_map)
  if (n_obs == 0L) {
    stop("All observations removed after missing-data handling. ",
         "No complete lag window exists within a subject/segment. ",
         "Check missingness, K and the scheduled index spacing.", call. = FALSE)
  }
  Y <- as.matrix(data[row_map, y_cols, drop = FALSE])
  X <- as.matrix(data[row_map, x_cols, drop = FALSE])
  B <- do.call(cbind, lapply(seq_len(K), function(lag)
    as.matrix(data[predecessors[row_map, lag], y_cols, drop = FALSE])))
  rownames(Y) <- rownames(X) <- rownames(B) <- NULL
  id_out <- id_all[row_map]
  time_obs <- data[[time_col]][row_map]
  index_obs <- data[row_map, keys, drop = FALSE]
  rownames(index_obs) <- NULL
  message(sprintf(paste0("bvarnet: %d incomplete occasion(s) and %d unavailable lag window(s) ",
                         "excluded; %d responses remain."),
                  n_input - n_rows_data, n_rows_data - n_obs, n_obs))
  empty_subjects <- J - length(unique(id_out))
  if (empty_subjects > 0L)
    message(sprintf("bvarnet: %d subject(s) have no contributing responses.", empty_subjects))

  x_center_means <- NULL
  if (isTRUE(center_x)) {
    means_X <- colMeans(X)
    X <- sweep(X, 2, means_X, "-")
    x_center_means <- means_X
  }

  # Always add intercept (shared convention)
  X <- cbind(Intercept = 1, X)
  colnames(X) <- c("Intercept", x_cols)

  colnames(Y) <- y_cols
  b_names <- unlist(lapply(1:K, function(lag) paste0("lag", lag, "_", y_cols)))
  colnames(B) <- b_names

  ## build FE interactions
  tmp <- add_terms_to_X(X, B, fe_interactions)
  X <- tmp$X
  fe_interaction_terms    <- normalize_terms(fe_interactions)
  fe_interaction_colnames <- tmp$new_names

  Z <- build_Z(X, B, re_cols = re_cols, re_temporal = re_temporal)
  Z <- add_re_interactions_from_X(Z, X, B, re_interactions)

  list(
    p = p, J = J, K = K, n_obs = n_obs,
    n_fe = ncol(X), n_re = ncol(Z),
    id = id_out, time_obs = time_obs, index_obs = index_obs, Y = Y, X = X, B = B, Z = Z,
    id_levels = as.character(ids_unique),
    x_center_means = x_center_means,
    row_map = row_map, n_rows_data = n_rows_data,
    design_spec = list(
      id_col = id_col, time_col = time_col, beep_col = beep_col,
      y_cols = y_cols, x_cols = x_cols,
      center_x = center_x,
      fe_interactions = fe_interactions,
      re_interactions = re_interactions,
      re_cols = re_cols, re_temporal = re_temporal,
      K = K, na_action = na_action
    ),
    fe_interaction_terms = fe_interaction_terms,
    fe_interaction_colnames = fe_interaction_colnames
  )
}


# ---- .to_stan_data_node() — per-node Stan data for node-wise fitting ----

#' Build Stan data for a single node from shared matrices
#'
#' @param shared List from \code{.to_stan_data_shared()}.
#' @param node Integer, which node (1..p).
#' @param family Character scalar, family for this node.
#' @param priors A \code{bvarnet_priors} object (original, unmodified).
#' @return Named list ready for \code{CmdStanModel$sample()}.
#' @keywords internal
#' @noRd
.to_stan_data_node <- function(shared, node, family, priors) {
  Y_node <- shared$Y[, node]

  if (family == "ordinal") {
    if (any(Y_node != floor(Y_node), na.rm = TRUE)) {
      stop("Ordinal Y must contain integer values. Found non-integer entries.")
    }
    Y_out <- matrix(as.integer(Y_node), ncol = 1L)
    X_out <- shared$X[, -1, drop = FALSE]  # strip intercept
    C     <- max(Y_out)
    n_fe  <- ncol(X_out)
  } else if (family == "bernoulli") {
    Y_out <- matrix(as.integer(Y_node), ncol = 1L)
    X_out <- shared$X
    n_fe  <- ncol(X_out)
  } else {
    Y_out <- matrix(Y_node, ncol = 1L)
    X_out <- shared$X
    n_fe  <- ncol(X_out)
  }

  # Rebuild Z from the (possibly stripped) X for this node
  Z_out <- build_Z(X_out, shared$B, re_cols = shared$design_spec$re_cols,
                    re_temporal = shared$design_spec$re_temporal)
  Z_out <- add_re_interactions_from_X(Z_out, X_out, shared$B,
                                       shared$design_spec$re_interactions)

  # Per-node Gaussian prior scaling (D3)
  node_priors <- .scale_default_priors(priors, family, s_y = sd(Y_node, na.rm = TRUE))

  # Abused K: Stan sees p=1, so K_node = p_full * K_orig to keep B correct
  K_node <- shared$p * shared$K

  out <- list(
    p      = 1L,
    J      = shared$J,
    K      = K_node,
    n_obs  = shared$n_obs,
    n_fe   = n_fe,
    n_re   = ncol(Z_out),
    id     = shared$id,
    Y      = Y_out,
    X      = X_out,
    B      = shared$B,
    Z      = Z_out,
    prior_beta_fam = node_priors$beta$family_int,
    beta_loc  = node_priors$beta$loc,
    beta_scale = node_priors$beta$scale,
    beta_df   = node_priors$beta$df,
    prior_phi_fam = node_priors$phi$family_int,
    phi_loc   = node_priors$phi$loc,
    phi_scale = node_priors$phi$scale,
    phi_df    = node_priors$phi$df,
    prior_sd_fam = node_priors$sd_u$family_int,
    sd_loc    = node_priors$sd_u$loc,
    sd_scale  = node_priors$sd_u$scale,
    sd_df     = node_priors$sd_u$df
  )

  if (family == "gaussian") {
    out$prior_sigma_fam <- node_priors$sigma$family_int
    out$sigma_loc  <- node_priors$sigma$loc
    out$sigma_scale <- node_priors$sigma$scale
    out$sigma_df   <- node_priors$sigma$df
  }

  if (family %in% c("gaussian", "bernoulli")) {
    out$prior_intercept_fam <- node_priors$intercept$family_int
    out$intercept_loc       <- node_priors$intercept$loc
    out$intercept_scale     <- node_priors$intercept$scale
    out$intercept_df        <- node_priors$intercept$df
  }

  if (family == "ordinal") {
    out$C <- C
    out$prior_kappa_fam <- node_priors$kappa$family_int
    out$kappa_loc  <- node_priors$kappa$loc
    out$kappa_scale <- node_priors$kappa$scale
    out$kappa_df   <- node_priors$kappa$df
  }

  # See to_stan_data(): carried as an attribute so it never reaches Stan's JSON.
  attr(out, "priors_effective") <- node_priors

  out
}

build_Z <- function(X, B, re_cols = character(0), re_temporal = FALSE) {
  # X: n_obs x n_fe (incl intercept for bernoulli/gaussian; NO intercept for ordinal)
  # B: n_obs x (p*K)
  stopifnot(is.matrix(X), is.matrix(B))

  Z_list <- list()

  # random slopes on selected fixed-effect columns
  if (length(re_cols) > 0) {
    # Handle "Intercept" specially: for ordinal, X has no intercept column
    # (intercept is absorbed by the kappa thresholds), so create a column of 1s
    if ("Intercept" %in% re_cols && !"Intercept" %in% colnames(X)) {
      intercept_col <- matrix(1, nrow = nrow(X), ncol = 1)
      colnames(intercept_col) <- "Intercept"
      X <- cbind(X, intercept_col)
    }

    missing <- setdiff(re_cols, colnames(X))
    if (length(missing) > 0) {
      stop("re_cols not found in X: ", paste(missing, collapse = ", "))
    }
    Z_list[["X"]] <- X[, re_cols, drop = FALSE]
  }

  # random slopes on lag predictors (temporal structure)
  if (isTRUE(re_temporal)) {
    Z_list[["B"]] <- B  # all lag columns
  }

  if (length(Z_list) == 0) {
    # no random effects
    Z <- matrix(0.0, nrow(X), 0)
  } else {
    Z <- do.call(cbind, Z_list)
  }

  # ensure matrix + colnames
  Z <- as.matrix(Z)
  if (ncol(Z) > 0 && is.null(colnames(Z))) colnames(Z) <- paste0("re", seq_len(ncol(Z)))

  return(Z)
}

# Todo: Interaction effects!!!
# Todo: validate design matrix building
# Todo: validate missing data mechanism


### For debugging
# K <- 1
# data <- df
# id_col <- "id"
# time_col <- "t"
# y_cols <- c("y_1", "y_2", "y_3", "y_4")
# x_cols <- c("x_1", "x_2")
#
#
#
#
# fe_terms = list(
#   c("x_1","x_2"),                 # 2-way
#   c("x_1","x_2","lag"),          # 3-way
#   c("lag","x_1")                 # lag moderation (whole structure)
# )
#
# re_terms = list(
#   c("x_1","x_2"),                # random slope on interaction
#   c("lag","x_1")                  # random slopes on whole lag structure moderated by x_1
# )

normalize_terms <- function(terms) {
  if (is.null(terms) || length(terms) == 0) return(list())

  if (!is.list(terms)) stop("terms must be a list of character vectors, e.g. list(c('x1','x2'), c('lag','x1'))")

  #if ("lag" %in% t) t <- c("lag", t[t != "lag"]) # fix, put lags in front always

  out <- lapply(terms, function(t) {
    if (!is.character(t)) stop("Each term must be a character vector.")
    t <- trimws(t)
    if (length(t) < 2) stop("Each term must have at least 2 factors (e.g. c('x1','x2') or c('lag','x1')).")

    if (any(t == "")) stop("Empty factor name in term.")
    if (any(duplicated(t[t != "lag"]))) stop("Duplicate non-lag factor in term: ", paste(t, collapse=":"))

    t
  })

  out
}



add_terms_to_X <- function(X, B, terms) {
  terms <- normalize_terms(terms)
  if (length(terms) == 0) return(list(X = X, new_names = character(0)))

  blocks <- lapply(terms, function(f) make_term_matrix(X, B, f))
  W <- do.call(cbind, blocks)

  dup <- intersect(colnames(W), colnames(X))
  if (length(dup) > 0) stop("Fixed-effect interaction columns already exist in X: ", paste(dup, collapse=", "))

  X2 <- cbind(X, W)
  list(X = X2, new_names = colnames(W))
}

make_term_matrix <- function(Xc, B, factors) {
  has_lag <- any(factors == "lag")
  others  <- factors[factors != "lag"]

  # build multiplicative modifier from Xc columns
  if (length(others) == 0) {
    mod <- rep(1.0, nrow(Xc))
    suffix <- ""
  } else {
    missing <- setdiff(others, colnames(Xc))
    if (length(missing) > 0) stop("Unknown factor(s): ", paste(missing, collapse=", "))
    mod <- rep(1.0, nrow(Xc))
    for (v in others) mod <- mod * as.numeric(Xc[, v])
    suffix <- paste(others, collapse=":")
  }

  if (!has_lag) {
    M <- matrix(mod, nrow(Xc), 1)
    colnames(M) <- suffix
    return(M)
  }

  # lag expansion = whole temporal structure
  M <- B * mod
  colnames(M) <- if (suffix == "") colnames(B) else paste0(colnames(B), ":", suffix)
  M
}



required_fe_names_for_term <- function(B, factors) {
  has_lag <- any(factors == "lag")
  others  <- factors[factors != "lag"]
  suffix  <- if (length(others) == 0) "" else paste(others, collapse=":")

  if (!has_lag) {
    return(suffix)
  } else {
    if (is.null(colnames(B))) stop("B must have colnames for lag expansion.")
    if (suffix == "") return(colnames(B))
    return(paste0(colnames(B), ":", suffix))
  }
}

required_fe_names_for_re_terms <- function(B, re_terms) {
  re_terms <- normalize_terms(re_terms)
  unique(unlist(lapply(re_terms, function(f) required_fe_names_for_term(B, f))))
}

add_re_interactions_from_X <- function(Z, Xc, B, re_terms) {
  if (is.null(re_terms) || length(re_terms) == 0) return(Z)

  want <- required_fe_names_for_re_terms(B, re_terms)
  missing <- setdiff(want, colnames(Xc))
  if (length(missing) > 0) {
    stop(
      "Random-effect interaction(s) requested but not present in fixed effects (X). Missing: ",
      paste(missing, collapse = ", "),
      "\nAdd them to `interactions` / `lag_interactions` first."
    )
  }

  cbind(Z, Xc[, want, drop = FALSE])
}
