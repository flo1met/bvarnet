# ──────────────────────────────────────────────────────────────────────────────
# simulate.bvarnet  &  predict.bvarnet  — S3 methods
# ──────────────────────────────────────────────────────────────────────────────


# ═══════════════════════════════════════════════════════════════════════════════
#                          INTERNAL HELPERS
# ═══════════════════════════════════════════════════════════════════════════════

# --- extract draw-level parameter blocks -----------------------------------

#' Extract a single posterior draw (or posterior mean) for a parameter block
#'
#' @param object A `bvarnet` object.
#' @param parameter Character — one of "beta", "phi", "sigma", "kappa", "sd_u".
#' @param draw_index `NULL` for posterior mean, integer for a specific draw
#'   (1-based, over all chains flattened).
#' @return Named numeric vector (flat, CmdStan ordering).
#' @noRd
.extract_param_draw <- function(object, parameter, draw_index = NULL) {
  draws <- object$draws  # [iter, chain, param]
  idx   <- grep(paste0("^", parameter, "\\["), dimnames(draws)[[3]])
  if (length(idx) == 0L)
    stop(sprintf("No '%s[...]' parameters found in draws.", parameter),
         call. = FALSE)

  chunk <- draws[, , idx, drop = FALSE]
  S     <- prod(dim(chunk)[1:2])
  dim(chunk) <- c(S, length(idx))
  colnames(chunk) <- dimnames(draws)[[3]][idx]

  if (is.null(draw_index)) {
    # posterior mean
    colMeans(chunk)
  } else {
    stopifnot(is.numeric(draw_index), length(draw_index) == 1L,
              draw_index >= 1L, draw_index <= S)
    chunk[draw_index, ]
  }
}


#' Reshape flat beta vector into matrix[n_fe, p]
#'
#' CmdStan column-major: beta[1,1], beta[2,1], ..., beta[n_fe,1], beta[1,2], ...
#' @noRd
.reshape_beta <- function(beta_vec, standata) {
  n_fe <- standata$n_fe
  p    <- standata$p
  matrix(beta_vec, nrow = n_fe, ncol = p)  # column-major matches CmdStan
}


#' Reshape flat phi vector into matrix[p*K, p]
#' @noRd
.reshape_phi <- function(phi_vec, standata) {
  p  <- standata$p
  K  <- standata$K
  PK <- p * K
  matrix(phi_vec, nrow = PK, ncol = p)
}


#' Reshape flat sigma vector into numeric vector of length p
#' @noRd
.reshape_sigma <- function(sigma_vec, standata) {
  as.numeric(sigma_vec)
}


#' Reshape flat kappa vector into list of p ordered vectors of length C-1
#'
#' CmdStan for `array[p] ordered[C-1] kappa`:
#'   kappa[1,1], kappa[1,2], ..., kappa[1,C-1], kappa[2,1], ...
#' i.e. array index (node) varies slowest, ordered index varies fastest.
#' @noRd
.reshape_kappa <- function(kappa_vec, standata) {
  p <- standata$p
  Cm1 <- standata$C - 1L
  # CmdStan flat ordering: node varies slowest, cutpoint varies fastest
  m <- matrix(kappa_vec, nrow = Cm1, ncol = p)
  lapply(seq_len(p), function(node) m[, node])
}


#' Reshape flat sd_u vector into matrix[p, n_re]
#' @noRd
.reshape_sd_u <- function(sd_u_vec, standata) {
  p    <- standata$p
  n_re <- standata$n_re
  matrix(sd_u_vec, nrow = p, ncol = n_re)
}


# --- build prediction standata from newdata --------------------------------

#' Build design matrices from newdata using the same recipe as fit-time
#'
#' Re-runs the same data-wrangling logic that `to_stan_data()` uses:
#' ordering, NA removal, lagging, centering, interactions, Z construction.
#' Returns a list with X, B, Z, id_char (character IDs for subject matching),
#' row_map (indices into original newdata), and n_rows (nrow(newdata)).
#'
#' @param newdata Data frame in long format.
#' @param object A `bvarnet` object.
#' @return Named list with X, B, Z, id_char, row_map, n_rows.
#' @noRd
.build_pred_standata <- function(newdata, object) {
  sd   <- object$standata
  spec <- sd$design_spec

  id_col   <- spec$id_col
  time_col <- spec$time_col
  y_cols   <- spec$y_cols
  x_cols   <- spec$x_cols
  K        <- spec$K
  skip_lag <- spec$skip_lag
  family   <- object$family
  p        <- sd$p

  # --- validate required columns -----------------------------------------
  required <- c(id_col, time_col, y_cols, x_cols)
  missing_cols <- setdiff(required, colnames(newdata))
  if (length(missing_cols) > 0L)
    stop("Missing columns in newdata: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  # Keep original row indices for mapping back
  newdata$.orig_row <- seq_len(nrow(newdata))

  # --- order by id, time -------------------------------------------------
  newdata <- newdata[order(newdata[[id_col]], newdata[[time_col]]), ]

  # --- split by subject and build matrices --------------------------------
  df_split <- split(newdata, newdata[[id_col]])
  ids_unique <- unique(newdata[[id_col]])
  J_new <- length(ids_unique)

  # Count modeled rows
  n_obs <- sum(vapply(df_split, function(sub) max(0L, nrow(sub) - K), integer(1)))
  if (n_obs == 0L) stop("No predictable rows after dropping first K per subject.",
                         call. = FALSE)

  PK <- p * K
  q  <- length(x_cols)

  X  <- matrix(NA_real_, n_obs, q)
  B  <- matrix(0, n_obs, PK)
  id_char  <- character(n_obs)
  row_map  <- integer(n_obs)

  row <- 0L
  for (jj in seq_len(J_new)) {
    this_id <- ids_unique[jj]
    df_sub  <- df_split[[as.character(this_id)]]
    df_sub  <- df_sub[order(df_sub[[time_col]]), , drop = FALSE]
    Ti      <- nrow(df_sub)
    if (Ti <= K) next

    times <- df_sub[[time_col]]
    Ymat  <- as.matrix(df_sub[, y_cols, drop = FALSE])
    Xmat  <- as.matrix(df_sub[, x_cols, drop = FALSE])

    for (t in (K + 1L):Ti) {
      row <- row + 1L
      id_char[row] <- as.character(this_id)
      row_map[row] <- df_sub$.orig_row[t]
      X[row, ] <- Xmat[t, ]

      valid <- TRUE
      for (lag in 1:K) {
        if ((times[t] - times[t - lag]) != lag) {
          valid <- FALSE
          break
        }
      }
      if (valid) {
        for (lag in 1:K) {
          B[row, ((lag - 1L) * p + 1L):(lag * p)] <- Ymat[t - lag, ]
        }
      } else if (!skip_lag) {
        row <- row - 1L
        next
      }
    }
  }

  # Trim
  if (row < n_obs) {
    X       <- X[seq_len(row), , drop = FALSE]
    B       <- B[seq_len(row), , drop = FALSE]
    id_char <- id_char[seq_len(row)]
    row_map <- row_map[seq_len(row)]
    n_obs   <- row
  }

  # --- centering ---------------------------------------------------------
  if (isTRUE(spec$center_x) && !is.null(sd$x_center_means)) {
    X <- sweep(X, 2, sd$x_center_means, "-")
  }

  colnames(X) <- x_cols

  # --- intercept (family-specific) ----------------------------------------
  if (family %in% c("bernoulli", "gaussian")) {
    X <- cbind(Intercept = 1, X)
    colnames(X) <- c("Intercept", x_cols)
  }

  # --- lag column names ---------------------------------------------------
  b_names <- unlist(lapply(1:K, function(lag) paste0("lag", lag, "_", y_cols)))
  colnames(B) <- b_names

  # --- FE interactions ----------------------------------------------------
  tmp <- add_terms_to_X(X, B, spec$fe_interactions)
  X <- tmp$X

  # --- RE design matrix Z -------------------------------------------------
  Z <- build_Z(X, B, re_cols = spec$re_cols, re_temporal = spec$re_temporal)
  Z <- add_re_interactions_from_X(Z, X, B, spec$re_interactions)

  list(
    X        = X,
    B        = B,
    Z        = Z,
    id_char  = id_char,
    row_map  = row_map,
    n_rows   = nrow(newdata)
  )
}


# --- compute linear predictor per node ------------------------------------

#' Compute eta (linear predictor) for all nodes
#'
#' @param X  Matrix [n_obs, n_fe].
#' @param B  Matrix [n_obs, p*K].
#' @param beta Matrix [n_fe, p].
#' @param phi  Matrix [p*K, p].
#' @param Z    Matrix [n_obs, n_re] (or NULL / 0-col).
#' @param u_rows Matrix [n_obs, n_re] — per-row RE values already resolved
#'        for each node call; or NULL for no RE.
#' @param node  Integer, which node/column (1..p).
#' @return Numeric vector of length n_obs.
#' @noRd
.predict_eta_node <- function(X, B, beta, phi, Z, u_rows, node) {
  # Fixed part:  X %*% beta[, node] + B %*% phi[, node]
  eta <- as.numeric(X %*% beta[, node]) +
         as.numeric(B %*% phi[, node])

  # Random effect part
  if (!is.null(u_rows) && ncol(Z) > 0L) {
    eta <- eta + rowSums(Z * u_rows)
  }
  eta
}


# --- resolve RE for each modeled row --------------------------------------

#' Resolve u values for each row based on subject_re policy
#'
#' @param object A `bvarnet` object.
#' @param id_char Character vector of subject IDs (length n_modeled_rows).
#' @param subject_re Character: "zero", "posterior-mean", "sample".
#' @param new_subject Character: "zero", "sample".
#' @param draw_index NULL or integer (for sample mode).
#' @param node Integer node index.
#' @return Matrix [n_rows, n_re].
#' @noRd
.get_re_for_rows <- function(object, id_char, subject_re, new_subject,
                             draw_index, node) {
  sd   <- object$standata
  n_re <- sd$n_re
  n    <- length(id_char)

  if (n_re == 0L || subject_re == "zero") {
    return(matrix(0, n, max(n_re, 1L)))
  }

  id_levels <- sd$id_levels
  j_idx     <- match(id_char, id_levels)  # NA if unseen

  out <- matrix(0, n, n_re)

  if (subject_re == "posterior-mean") {
    u_mean <- .posterior_mean_u(object)  # [p, J, n_re]
    for (r in seq_len(n)) {
      if (!is.na(j_idx[r])) {
        out[r, ] <- u_mean[node, j_idx[r], ]
      } else {
        # unseen subject
        if (new_subject == "sample") {
          sd_u_mean <- colMeans(.reshape_sd_u(
            .extract_param_draw(object, "sd_u", NULL), sd))
          out[r, ] <- rnorm(n_re, 0, sd_u_mean[node, ])
        }
        # else "zero" → leave as 0
      }
    }
  } else if (subject_re == "sample") {
    stopifnot(!is.null(draw_index))
    u_draws <- .extract_u_draws(object)  # [S, p, J, n_re]
    sd_u_draw <- .reshape_sd_u(
      .extract_param_draw(object, "sd_u", draw_index), sd)
    for (r in seq_len(n)) {
      if (!is.na(j_idx[r])) {
        out[r, ] <- u_draws[draw_index, node, j_idx[r], ]
      } else {
        if (new_subject == "sample") {
          out[r, ] <- rnorm(n_re, 0, sd_u_draw[node, ])
        }
      }
    }
  }

  out
}


# --- eta to output conversion ---------------------------------------------

#' Convert linear predictor matrix to the requested output type
#'
#' @param eta_mat Matrix [n_obs, p] of linear predictor values.
#' @param family Character.
#' @param type Character: "link", "response", "probabilities".
#' @param sigma Numeric vector length p (gaussian only).
#' @param kappa List of p ordered vectors (ordinal only).
#' @return For "link"/"response": matrix [n_obs, p].
#'   For "probabilities": list of p matrices.
#' @noRd
.eta_to_output <- function(eta_mat, family, type, sigma = NULL, kappa = NULL) {
  n   <- nrow(eta_mat)
  p   <- ncol(eta_mat)

  if (type == "link") {
    return(eta_mat)
  }

  if (type == "response") {
    if (family == "bernoulli") {
      # inverse logit
      return(1 / (1 + exp(-eta_mat)))
    } else if (family == "gaussian") {
      return(eta_mat)  # identity link
    } else if (family == "ordinal") {
      # expected value = sum( c * P(Y=c) )
      C <- length(kappa[[1]]) + 1L
      out <- matrix(NA_real_, n, p)
      for (node in seq_len(p)) {
        probs <- .ordinal_probs(eta_mat[, node], kappa[[node]], C)
        out[, node] <- probs %*% seq_len(C)
      }
      return(out)
    }
  }

  if (type == "probabilities") {
    if (family == "bernoulli") {
      prob1 <- 1 / (1 + exp(-eta_mat))
      out <- vector("list", p)
      for (node in seq_len(p)) {
        m <- matrix(prob1[, node], ncol = 1)
        colnames(m) <- "p1"
        out[[node]] <- m
      }
      return(out)
    } else if (family == "gaussian") {
      out <- vector("list", p)
      for (node in seq_len(p)) {
        m <- cbind(mean = eta_mat[, node], sd = rep(sigma[node], n))
        out[[node]] <- m
      }
      return(out)
    } else if (family == "ordinal") {
      C <- length(kappa[[1]]) + 1L
      out <- vector("list", p)
      for (node in seq_len(p)) {
        probs <- .ordinal_probs(eta_mat[, node], kappa[[node]], C)
        colnames(probs) <- paste0("cat_", seq_len(C))
        out[[node]] <- probs
      }
      return(out)
    }
  }

  stop("Unknown type: ", type, call. = FALSE)
}


#' Compute adjacent-category probabilities for a vector of etas
#'
#' Matches Stan parameterisation:
#'   lambda[c] = (c-1)*eta - cumsum(kappa)[c];  P(Y=c) = softmax(lambda)
#' @param eta Numeric vector of length n.
#' @param kappa_node Ordered numeric vector of length C-1.
#' @param C Integer, number of categories.
#' @return Matrix [n, C] of probabilities.
#' @noRd
.ordinal_probs <- function(eta, kappa_node, C) {
  n <- length(eta)
  kappa_cumsum <- c(0, cumsum(kappa_node))
  probs <- matrix(NA_real_, n, C)

  for (i in seq_len(n)) {
    lambda <- (seq_len(C) - 1L) * eta[i] - kappa_cumsum
    lambda <- lambda - max(lambda)  # log-sum-exp stability
    probs[i, ] <- exp(lambda) / sum(exp(lambda))
  }
  probs
}


# --- expand compact results to full newdata rows ---------------------------

#' Expand modeled-row results back to nrow(newdata) with NA for dropped rows
#'
#' @param compact Matrix or list.
#' @param row_map Integer vector mapping modeled rows to original row indices.
#' @param n_rows Integer, nrow(newdata).
#' @return Object of the same type with NA-padded rows.
#' @noRd
.expand_to_full_rows <- function(compact, row_map, n_rows) {
  if (is.matrix(compact)) {
    out <- matrix(NA_real_, n_rows, ncol(compact))
    colnames(out) <- colnames(compact)
    out[row_map, ] <- compact
    return(out)
  }
  if (is.list(compact)) {
    # list of matrices (probabilities output)
    return(lapply(compact, function(m) {
      out <- matrix(NA_real_, n_rows, ncol(m))
      colnames(out) <- colnames(m)
      out[row_map, ] <- m
      out
    }))
  }
  stop("Unexpected compact type", call. = FALSE)
}


# --- draw index helpers ----------------------------------------------------

#' Get total number of posterior draws (S = iter * chains)
#' @noRd
.n_draws_total <- function(object) {
  prod(dim(object$draws)[1:2])
}

#' Sample draw indices (1..S)
#' @noRd
.sample_draw_indices <- function(object, ndraws, seed = NULL) {
  S <- .n_draws_total(object)
  ndraws <- min(ndraws, S)
  if (!is.null(seed)) set.seed(seed)
  sort(sample.int(S, ndraws))
}


# ═══════════════════════════════════════════════════════════════════════════════
#                       predict.bvarnet
# ═══════════════════════════════════════════════════════════════════════════════

#' Predict from a fitted bvarnet model
#'
#' Computes one-step-ahead predictions for long-format time-series data.
#' Supports population-level (\code{subject_re = "zero"}) and
#' subject-specific (\code{subject_re = "posterior-mean"}) predictions.
#' Also serves as the out-of-sample engine: fit on training data, call
#' \code{predict(fit, newdata = test_data)}.
#'
#' @param object A \code{bvarnet} object from \code{bvar()}.
#' @param newdata Data frame in long format. If \code{NULL}, the original
#'   training data design matrices (stored in \code{object$standata}) are used
#'   for in-sample fitted values.
#' @param type Character. Output type: \code{"link"} (linear predictor),
#'   \code{"response"} (mean on the outcome scale), or
#'   \code{"probabilities"} (category-level probabilities/mean+sd).
#' @param method Character. \code{"posterior-mean"} uses posterior means
#'   only (deterministic). \code{"posterior-sample"} averages over
#'   \code{ndraws} draws and returns \code{attr(,"sd")} with across-draw SD.
#' @param ndraws Integer. Number of posterior draws to use when
#'   \code{method = "posterior-sample"}. Defaults to the smaller of 100 and the
#'   total available draws.
#' @param seed Integer or NULL. RNG seed for draw selection and new-subject
#'   sampling.
#' @param subject_re Character. How to handle random effects:
#'   \code{"zero"} (population-level, u = 0),
#'   \code{"posterior-mean"} (posterior mean of u for seen subjects), or
#'   \code{"sample"} (draw-specific u for seen subjects).
#' @param new_subject Character. Fallback for unseen subjects:
#'   \code{"zero"} (u = 0) or
#'   \code{"sample"} (draw from RE distribution).
#' @param ... Ignored.
#'
#' @return For \code{type = "link"} or \code{"response"}: a numeric matrix with
#'   \code{nrow(newdata)} rows and \code{p} columns, \code{NA} for the first K
#'   rows per subject. For \code{type = "probabilities"}: a list of \code{p}
#'   matrices. When \code{method = "posterior-sample"}, the output carries
#'   \code{attr(,"sd")} and \code{attr(,"ndraws")}.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   sim <- sim_var(N = 5, T_obs = 30, p = 2, family = "gaussian", seed = 1)
#'   fit <- bvar(id_col = "id", time_col = "t",
#'               y_cols = c("y_1", "y_2"), x_cols = character(0),
#'               data = sim$data, family = "gaussian",
#'               iter = 200, warmup = 100, chains = 2, seed = 1)
#'   preds <- predict(fit, type = "response")
#' }
#' @importFrom stats predict
#' @export
predict.bvarnet <- function(object,
                            newdata     = NULL,
                            type        = c("link", "response", "probabilities"),
                            method      = c("posterior-mean", "posterior-sample"),
                            ndraws      = NULL,
                            seed        = NULL,
                            subject_re  = c("zero", "posterior-mean", "sample"),
                            new_subject = c("zero", "sample"),
                            ...) {

  type        <- match.arg(type)
  method      <- match.arg(method)
  subject_re  <- match.arg(subject_re)
  new_subject <- match.arg(new_subject)

  sd     <- object$standata
  family <- object$family
  p      <- sd$p
  n_re   <- sd$n_re

  # Degrade subject_re to "zero" when no RE

  if (n_re == 0L) subject_re <- "zero"

  # --- validate type vs family ---
  if (type == "probabilities" && family == "gaussian") {
    # gaussian probabilities returns mean + sd columns — this is fine
  }

  # --- build design matrices / use stored ---
  if (is.null(newdata)) {
    # In-sample: use stored matrices
    X       <- sd$X
    B       <- sd$B
    Z       <- sd$Z
    id_char <- as.character(sd$id_levels[sd$id])
    row_map <- NULL
    n_rows  <- NULL
  } else {
    pred_sd <- .build_pred_standata(newdata, object)
    X       <- pred_sd$X
    B       <- pred_sd$B
    Z       <- pred_sd$Z
    id_char <- pred_sd$id_char
    row_map <- pred_sd$row_map
    n_rows  <- pred_sd$n_rows
  }

  n_obs <- nrow(X)

  # --- dispatch by method ------------------------------------------------
  if (method == "posterior-mean") {
    # Single deterministic prediction
    beta  <- .reshape_beta(.extract_param_draw(object, "beta", NULL), sd)
    phi   <- .reshape_phi(.extract_param_draw(object, "phi", NULL), sd)
    sigma <- if (family == "gaussian") .reshape_sigma(
      .extract_param_draw(object, "sigma", NULL), sd) else NULL
    kappa <- if (family == "ordinal") .reshape_kappa(
      .extract_param_draw(object, "kappa", NULL), sd) else NULL

    eta_mat <- matrix(NA_real_, n_obs, p)
    for (node in seq_len(p)) {
      u_rows <- .get_re_for_rows(object, id_char, subject_re,
                                  new_subject, NULL, node)
      eta_mat[, node] <- .predict_eta_node(X, B, beta, phi, Z, u_rows, node)
    }
    colnames(eta_mat) <- colnames(sd$Y)

    result <- .eta_to_output(eta_mat, family, type, sigma, kappa)

    # Name columns for matrix outputs
    if (is.matrix(result)) colnames(result) <- colnames(sd$Y)

  } else {
    # --- posterior-sample: average over draws
    S <- .n_draws_total(object)
    if (is.null(ndraws)) ndraws <- min(100L, S)
    draw_idx <- .sample_draw_indices(object, ndraws, seed)

    # Initialize accumulators
    if (type %in% c("link", "response")) {
      sum_out  <- matrix(0, n_obs, p)
      sum_sq   <- matrix(0, n_obs, p)
    } else {
      # probabilities: need per-node accumulators
      if (family == "bernoulli") {
        sum_out <- lapply(seq_len(p), function(.) matrix(0, n_obs, 1))
        sum_sq  <- lapply(seq_len(p), function(.) matrix(0, n_obs, 1))
      } else if (family == "gaussian") {
        sum_out <- lapply(seq_len(p), function(.) matrix(0, n_obs, 2))
        sum_sq  <- lapply(seq_len(p), function(.) matrix(0, n_obs, 2))
      } else {  # ordinal
        C <- sd$C
        sum_out <- lapply(seq_len(p), function(.) matrix(0, n_obs, C))
        sum_sq  <- lapply(seq_len(p), function(.) matrix(0, n_obs, C))
      }
    }

    for (s in draw_idx) {
      beta_s  <- .reshape_beta(.extract_param_draw(object, "beta", s), sd)
      phi_s   <- .reshape_phi(.extract_param_draw(object, "phi", s), sd)
      sigma_s <- if (family == "gaussian") .reshape_sigma(
        .extract_param_draw(object, "sigma", s), sd) else NULL
      kappa_s <- if (family == "ordinal") .reshape_kappa(
        .extract_param_draw(object, "kappa", s), sd) else NULL

      eta_s <- matrix(NA_real_, n_obs, p)
      for (node in seq_len(p)) {
        u_rows <- .get_re_for_rows(object, id_char, subject_re,
                                    new_subject, s, node)
        eta_s[, node] <- .predict_eta_node(X, B, beta_s, phi_s, Z, u_rows, node)
      }
      colnames(eta_s) <- colnames(sd$Y)

      out_s <- .eta_to_output(eta_s, family, type, sigma_s, kappa_s)

      if (type %in% c("link", "response")) {
        sum_out <- sum_out + out_s
        sum_sq  <- sum_sq  + out_s^2
      } else {
        for (node in seq_len(p)) {
          sum_out[[node]] <- sum_out[[node]] + out_s[[node]]
          sum_sq[[node]]  <- sum_sq[[node]]  + out_s[[node]]^2
        }
      }
    }

    nd <- length(draw_idx)

    if (type %in% c("link", "response")) {
      result    <- sum_out / nd
      result_sd <- sqrt(pmax(sum_sq / nd - result^2, 0))
      colnames(result)    <- colnames(sd$Y)
      colnames(result_sd) <- colnames(sd$Y)
      attr(result, "sd")     <- result_sd
      attr(result, "ndraws") <- nd
    } else {
      result    <- lapply(seq_len(p), function(node) sum_out[[node]] / nd)
      result_sd <- lapply(seq_len(p), function(node) {
        sqrt(pmax(sum_sq[[node]] / nd - (sum_out[[node]] / nd)^2, 0))
      })
      # Copy column names
      for (node in seq_len(p)) {
        colnames(result[[node]])    <- colnames(sum_out[[node]])
        colnames(result_sd[[node]]) <- colnames(sum_out[[node]])
      }
      attr(result, "sd")     <- result_sd
      attr(result, "ndraws") <- nd
    }
  }

  # --- expand to full newdata rows if needed ---
  if (!is.null(row_map)) {
    if (is.matrix(result)) {
      sd_attr <- attr(result, "sd")
      nd_attr <- attr(result, "ndraws")
      result <- .expand_to_full_rows(result, row_map, n_rows)
      if (!is.null(sd_attr))
        attr(result, "sd") <- .expand_to_full_rows(sd_attr, row_map, n_rows)
      if (!is.null(nd_attr)) attr(result, "ndraws") <- nd_attr
    } else {
      sd_attr <- attr(result, "sd")
      nd_attr <- attr(result, "ndraws")
      result <- .expand_to_full_rows(result, row_map, n_rows)
      if (!is.null(sd_attr))
        attr(result, "sd") <- .expand_to_full_rows(sd_attr, row_map, n_rows)
      if (!is.null(nd_attr)) attr(result, "ndraws") <- nd_attr
    }
  }

  result
}


# ═══════════════════════════════════════════════════════════════════════════════
#                      simulate.bvarnet
# ═══════════════════════════════════════════════════════════════════════════════

#' Simulate new trajectories from a fitted bvarnet model
#'
#' Generates new multilevel VAR trajectories using the posterior parameters.
#' Reuses the existing response kernels from \code{sim_var()}.
#'
#' @param object A \code{bvarnet} object from \code{bvar()}.
#' @param nsim Integer. Number of time points per subject. Default 20.
#' @param seed Integer or NULL. RNG seed.
#' @param method Character. \code{"posterior-mean"} returns one data frame.
#'   \code{"posterior-sample"} returns a list of data frames (one per draw).
#' @param ndraws Integer. Number of posterior draws when
#'   \code{method = "posterior-sample"}. Default 10.
#' @param N Integer. Number of subjects. Defaults to J from the fitted model.
#' @param burnin Integer. Burn-in time points to discard. Default 200.
#' @param x_gen Function \code{f(N, T)} returning an \code{N x T x q} array
#'   of covariates, or NULL for zero covariates.
#' @param subject_re Character. \code{"sample"} draws REs from
#'   \code{N(0, sd_u)}. \code{"zero"} sets all REs to 0.
#' @param ... Ignored.
#'
#' @return For \code{method = "posterior-mean"}: a data frame with columns
#'   \code{id}, \code{t}, \code{y_1}, ..., \code{y_p} (and optional
#'   \code{x_1}, ..., \code{x_q}). For \code{"posterior-sample"}: a list of
#'   such data frames, one per draw.
#'
#' @examples
#' if (instantiate::stan_cmdstan_exists()) {
#'   sim <- sim_var(N = 5, T_obs = 30, p = 2, family = "gaussian", seed = 1)
#'   fit <- bvar(id_col = "id", time_col = "t",
#'               y_cols = c("y_1", "y_2"), x_cols = character(0),
#'               data = sim$data, family = "gaussian",
#'               iter = 200, warmup = 100, chains = 2, seed = 1)
#'   sim_data <- simulate(fit, nsim = 20, seed = 42)
#' }
#' @importFrom stats simulate
#' @export
simulate.bvarnet <- function(object,
                             nsim       = 20L,
                             seed       = NULL,
                             method     = c("posterior-mean", "posterior-sample"),
                             ndraws     = 10L,
                             N          = NULL,
                             burnin     = 200L,
                             x_gen      = NULL,
                             subject_re = c("sample", "zero"),
                             ...) {

  method     <- match.arg(method)
  subject_re <- match.arg(subject_re)

  sd     <- object$standata
  family <- object$family
  p      <- sd$p
  K      <- sd$K
  n_re   <- sd$n_re
  n_fe   <- sd$n_fe
  PK     <- p * K

  if (n_re == 0L) subject_re <- "zero"
  if (is.null(N))    N <- sd$J
  if (!is.null(seed)) set.seed(seed)

  nsim   <- as.integer(nsim)
  burnin <- as.integer(burnin)
  N      <- as.integer(N)

  stopifnot(nsim >= K + 1L, N >= 1L, burnin >= 0L)

  T_total <- nsim + burnin

  # Determine q (number of covariates, excluding intercept for non-ordinal)
  x_cols <- sd$design_spec$x_cols
  q <- length(x_cols)

  # --- inner simulation function ----------------------------------------
  .simulate_one <- function(draw_index) {
    beta  <- .reshape_beta(.extract_param_draw(object, "beta", draw_index), sd)
    phi   <- .reshape_phi(.extract_param_draw(object, "phi", draw_index), sd)
    sigma <- if (family == "gaussian") .reshape_sigma(
      .extract_param_draw(object, "sigma", draw_index), sd) else NULL
    kappa <- if (family == "ordinal") .reshape_kappa(
      .extract_param_draw(object, "kappa", draw_index), sd) else NULL
    sd_u  <- if (n_re > 0L) .reshape_sd_u(
      .extract_param_draw(object, "sd_u", draw_index), sd) else NULL

    # Extract intercept and covariate effects from beta
    # beta is [n_fe, p]; for bernoulli/gaussian row 1 = intercept
    if (family %in% c("bernoulli", "gaussian")) {
      alpha <- beta[1, ]  # intercept
      gamma <- if (n_fe > 1L) beta[2:n_fe, , drop = FALSE] else matrix(0, 0, p)
    } else {
      # ordinal: no intercept in beta
      alpha <- rep(0, p)
      gamma <- if (n_fe > 0L) beta else matrix(0, 0, p)
    }

    # Generate covariates
    if (q > 0L) {
      if (!is.null(x_gen)) {
        X_cov <- x_gen(N, T_total)
        stopifnot(is.array(X_cov), identical(dim(X_cov), c(N, T_total, q)))
      } else {
        X_cov <- array(0, dim = c(N, T_total, q))
      }
    } else {
      X_cov <- array(0, dim = c(N, T_total, 0L))
    }

    # Generate subject-level random effects
    u_i <- matrix(0, N, n_re)
    if (subject_re == "sample" && n_re > 0L) {
      for (i in seq_len(N)) {
        for (re in seq_len(n_re)) {
          # sd_u is [p, n_re]; for simulation we need per-node sd
          # Use the mean across nodes for simplicity in generating u
          # Actually: u[node, subject, re] means each node has its own sd.
          # We generate u_i per (subject, re) then the node-level modulation
          # is handled by Z * u. But in the Stan model, u[node][subject, re]
          # is node-specific. For simulation, we need u per (node, subject, re).
        }
      }
    }

    # Storage: per-subject, per-node RE
    # u_sim[i, node, re]
    u_sim <- array(0, dim = c(N, p, n_re))
    if (subject_re == "sample" && n_re > 0L) {
      for (i in seq_len(N)) {
        for (node in seq_len(p)) {
          u_sim[i, node, ] <- rnorm(n_re, 0, sd_u[node, ])
        }
      }
    }

    # Forward simulation
    if (family == "gaussian") {
      Y_full <- array(NA_real_, dim = c(N, T_total, p))
    } else if (family == "ordinal") {
      Y_full <- array(NA_integer_, dim = c(N, T_total, p))
    } else {
      Y_full <- array(0L, dim = c(N, T_total, p))
    }

    for (i in seq_len(N)) {
      # Build alpha_i: intercept + RE contribution for the intercept RE
      # The RE design depends on which columns are in Z.
      # For simplicity, we use the direct linear predictor approach:
      # eta = alpha + gamma' * x + phi' * lag_y + u_contribution
      # where u_contribution per node = sum over re of (Z_row_value * u_sim[i, node, re])
      # But in simulation we don't have Z rows. The REs in the fitted model
      # correspond to specific columns. For intercept RE, Z has a column
      # "Intercept" if re_cols includes it. For lag REs, Z has lag columns.
      #
      # For simulation, we handle it differently: we fold REs into alpha_i and phi_i.
      # u_sim[i, node, re] * (corresponding coefficient from Z design)
      # Since Z is built from X and B columns, the REs modify those effects.

      # Determine which REs correspond to which predictors
      re_colnames <- colnames(sd$Z)
      fe_colnames <- colnames(sd$X)
      b_colnames  <- colnames(sd$B)

      alpha_i <- alpha
      phi_i   <- phi
      gamma_i <- gamma

      if (n_re > 0L) {
        for (re in seq_len(n_re)) {
          re_name <- re_colnames[re]
          # Check if this RE corresponds to an intercept
          if (re_name == "Intercept") {
            for (node in seq_len(p)) {
              alpha_i[node] <- alpha_i[node] + u_sim[i, node, re]
            }
          }
          # Check if it corresponds to a covariate (in X)
          else if (re_name %in% x_cols) {
            fe_idx_in_gamma <- match(re_name, x_cols)
            if (!is.na(fe_idx_in_gamma) && fe_idx_in_gamma <= nrow(gamma_i)) {
              for (node in seq_len(p)) {
                gamma_i[fe_idx_in_gamma, node] <-
                  gamma_i[fe_idx_in_gamma, node] + u_sim[i, node, re]
              }
            }
          }
          # Check if it corresponds to a lag (in B)
          else if (re_name %in% b_colnames) {
            lag_idx <- match(re_name, b_colnames)
            if (!is.na(lag_idx)) {
              for (node in seq_len(p)) {
                phi_i[lag_idx, node] <- phi_i[lag_idx, node] + u_sim[i, node, re]
              }
            }
          }
        }
      }

      # --- Initialize first K time points ---
      for (t in seq_len(K)) {
        eta <- alpha_i
        if (q > 0L) {
          g_mat <- matrix(gamma_i, nrow = q, ncol = p)
          x_vec <- matrix(X_cov[i, t, ], nrow = 1L)
          eta   <- eta + as.numeric(x_vec %*% g_mat)
        }
        Y_full[i, t, ] <- generate_response(eta, family, sigma, kappa,
                                             if (family == "ordinal") sd$C else NULL, p)
      }

      # --- Forward simulate t = K+1 .. T_total ---
      for (t in (K + 1L):T_total) {
        lag_y <- numeric(PK)
        for (lag in seq_len(K)) {
          idx <- ((lag - 1L) * p + 1L):(lag * p)
          lag_y[idx] <- Y_full[i, t - lag, ]
        }

        eta <- alpha_i + as.numeric(t(phi_i) %*% lag_y)

        if (q > 0L) {
          g_mat <- matrix(gamma_i, nrow = q, ncol = p)
          x_vec <- matrix(X_cov[i, t, ], nrow = 1L)
          eta   <- eta + as.numeric(x_vec %*% g_mat)
        }

        Y_full[i, t, ] <- generate_response(eta, family, sigma, kappa,
                                             if (family == "ordinal") sd$C else NULL, p)
      }
    }

    # Discard burnin
    keep_idx <- (burnin + 1L):T_total
    Y_keep <- Y_full[, keep_idx, , drop = FALSE]
    X_keep <- X_cov[, keep_idx, , drop = FALSE]

    # Assemble long-format data frame
    y_cols_out <- paste0("y_", seq_len(p))
    rows <- vector("list", N)
    for (i in seq_len(N)) {
      row_data <- data.frame(id = rep(i, nsim), t = seq_len(nsim))
      for (j in seq_len(p)) {
        row_data[[y_cols_out[j]]] <- Y_keep[i, , j]
      }
      if (q > 0L) {
        for (j in seq_len(q)) {
          row_data[[paste0("x_", j)]] <- X_keep[i, , j]
        }
      }
      rows[[i]] <- row_data
    }
    do.call(rbind, rows)
  }

  # --- dispatch by method ------------------------------------------------
  if (method == "posterior-mean") {
    return(.simulate_one(NULL))
  } else {
    S <- .n_draws_total(object)
    if (is.null(ndraws)) ndraws <- min(10L, S)
    draw_idx <- .sample_draw_indices(object, ndraws, seed)
    result <- lapply(draw_idx, .simulate_one)
    return(result)
  }
}
