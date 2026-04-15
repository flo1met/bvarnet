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
  Y_obs <- matrix(NA_real_, n_obs, p)
  id_char  <- character(n_obs)
  time_obs <- rep(NA_real_, n_obs)
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
      time_obs[row] <- times[t]
      row_map[row] <- df_sub$.orig_row[t]
      X[row, ] <- Xmat[t, ]
      Y_obs[row, ] <- Ymat[t, ]

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
    Y_obs   <- Y_obs[seq_len(row), , drop = FALSE]
    id_char <- id_char[seq_len(row)]
    time_obs <- time_obs[seq_len(row)]
    row_map <- row_map[seq_len(row)]
    n_obs   <- row
  }

  # --- centering ---------------------------------------------------------
  if (isTRUE(spec$center_x) && !is.null(sd$x_center_means)) {
    X <- sweep(X, 2, sd$x_center_means, "-")
  }

  colnames(X) <- x_cols

  # --- intercept (add only if stored standata has one) --------------------
  # Combined mixed standata (.to_stan_data_shared) always includes Intercept.
  # Homogeneous ordinal (to_stan_data) does not.
  if ("Intercept" %in% colnames(sd$X)) {
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
    Y_obs    = Y_obs,
    id_char  = id_char,
    time_obs = time_obs,
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
.predict_eta_node <- function(X, B, beta, phi, Z, u_rows, node,
                              family_node = "gaussian") {
  # Fixed part: skip intercept row/col for ordinal nodes (D4)
  # Only strip if X actually has an Intercept column (combined standata from
  # .to_stan_data_shared always includes it; homogeneous ordinal from
  # to_stan_data() already excludes it).
  has_intercept <- "Intercept" %in% colnames(X)
  if (family_node == "ordinal" && has_intercept) {
    icept_col <- which(colnames(X) == "Intercept")
    X_node    <- X[, -icept_col, drop = FALSE]
    beta_node <- beta[-1, node]
  } else {
    X_node    <- X
    beta_node <- beta[, node]
  }
  eta <- as.numeric(X_node %*% beta_node) +
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
          sd_u_mean <- .reshape_sd_u(
            .extract_param_draw(object, "sd_u", NULL), sd)
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
#' @param family_vec Character vector of length p (or scalar, recycled).
#' @param type Character: "link", "response", "probabilities".
#' @param sigma Numeric vector length p (NA for non-gaussian nodes).
#' @param kappa List of p elements (NULL for non-ordinal nodes).
#' @return For "link"/"response": matrix [n_obs, p].
#'   For "probabilities": list of p matrices.
#' @noRd
.eta_to_output <- function(eta_mat, family_vec, type, sigma = NULL, kappa = NULL) {
  n <- nrow(eta_mat)
  p <- ncol(eta_mat)

  if (length(family_vec) == 1L) family_vec <- rep(family_vec, p)

  if (type == "link") return(eta_mat)

  if (type == "response") {
    out <- matrix(NA_real_, n, p)
    for (j in seq_len(p)) {
      if (family_vec[j] == "bernoulli") {
        out[, j] <- 1 / (1 + exp(-eta_mat[, j]))
      } else if (family_vec[j] == "gaussian") {
        out[, j] <- eta_mat[, j]
      } else if (family_vec[j] == "ordinal") {
        C_j <- length(kappa[[j]]) + 1L
        probs <- .ordinal_probs(eta_mat[, j], kappa[[j]], C_j)
        out[, j] <- probs %*% seq_len(C_j)
      }
    }
    return(out)
  }

  if (type == "probabilities") {
    out <- vector("list", p)
    for (j in seq_len(p)) {
      if (family_vec[j] == "bernoulli") {
        prob1 <- 1 / (1 + exp(-eta_mat[, j]))
        m <- matrix(prob1, ncol = 1)
        colnames(m) <- "p1"
        out[[j]] <- m
      } else if (family_vec[j] == "gaussian") {
        out[[j]] <- cbind(mean = eta_mat[, j],
                          sd = rep(sigma[j], n))
      } else if (family_vec[j] == "ordinal") {
        C_j <- length(kappa[[j]]) + 1L
        probs <- .ordinal_probs(eta_mat[, j], kappa[[j]], C_j)
        colnames(probs) <- paste0("cat_", seq_len(C_j))
        out[[j]] <- probs
      }
    }
    return(out)
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


# --- recursive forecasting helpers -----------------------------------------

#' Validate and normalise conditioning_window into a per-subject named vector
#'
#' @param conditioning_window NULL, integer scalar, or named integer vector.
#' @param unique_ids Character vector of unique subject IDs (in row order).
#' @param K Lag order.
#' @return Named integer vector keyed by subject ID.
#' @noRd
.resolve_conditioning_window_by_subject <- function(conditioning_window,
                                                     unique_ids, K) {
  if (is.null(conditioning_window)) conditioning_window <- K

  if (length(conditioning_window) == 1L && is.null(names(conditioning_window))) {
    cw <- stats::setNames(rep(as.integer(conditioning_window), length(unique_ids)),
                          as.character(unique_ids))
  } else if (!is.null(names(conditioning_window))) {
    ids_chr <- as.character(unique_ids)
    cw_names <- as.character(names(conditioning_window))
    matched <- match(ids_chr, cw_names)
    if (any(is.na(matched)))
      stop("conditioning_window is missing entries for subjects: ",
           paste(ids_chr[is.na(matched)], collapse = ", "), call. = FALSE)
    cw <- stats::setNames(as.integer(conditioning_window[matched]), ids_chr)
  } else {
    stop("conditioning_window must be NULL, a scalar, or a named vector.",
         call. = FALSE)
  }

  if (any(cw < K))
    stop(sprintf("conditioning_window must be >= K (%d) for all subjects.", K),
         call. = FALSE)

  cw
}


#' Shift lag buffer and insert new observation
#'
#' @param lag_buffer Numeric vector of length p*K
#'   ordered [lag1_y1..lag1_yp, lag2_y1..lag2_yp, ...].
#' @param y_new Numeric vector of length p (new values to insert as lag 1).
#' @param K Integer lag order.
#' @param p Integer number of outcome variables.
#' @return Updated lag buffer.
#' @noRd
.update_lag_buffer <- function(lag_buffer, y_new, K, p) {
  if (K > 1L) {
    lag_buffer[(p + 1L):(p * K)] <- lag_buffer[1L:(p * (K - 1L))]
  }
  lag_buffer[1:p] <- y_new
  lag_buffer
}


#' Convert a row of eta to the value used for recursive lag updates
#'
#' @param eta_row Numeric vector of length p (linear predictor per node).
#' @param family_vec Character vector of length p (or scalar, recycled).
#' @param sigma Numeric vector (NA for non-gaussian) or NULL.
#' @param kappa List of p elements (NULL for non-ordinal) or NULL.
#' @return Numeric vector of length p.
#' @noRd
.recursive_lag_value <- function(eta_row, family_vec, sigma, kappa) {
  p <- length(eta_row)
  if (length(family_vec) == 1L) family_vec <- rep(family_vec, p)
  vals <- numeric(p)
  for (j in seq_len(p)) {
    if (family_vec[j] == "gaussian") {
      vals[j] <- eta_row[j]
    } else if (family_vec[j] == "bernoulli") {
      vals[j] <- 1 / (1 + exp(-eta_row[j]))
    } else if (family_vec[j] == "ordinal") {
      C_j <- length(kappa[[j]]) + 1L
      probs <- .ordinal_probs(eta_row[j], kappa[[j]], C_j)
      vals[j] <- sum(probs * seq_len(C_j))
    }
  }
  vals
}


#' Extract sigma (length-p vector, NA for non-gaussian) and kappa
#' (length-p list, NULL for non-ordinal) from a bvarnet object
#'
#' @param object A `bvarnet` object.
#' @param draw_index NULL for posterior mean, integer for a specific draw.
#' @return Named list with `sigma` and `kappa`.
#' @noRd
.extract_sigma_kappa <- function(object, draw_index = NULL) {
  p <- object$standata$p
  family_vec <- object$family

  sigma <- rep(NA_real_, p)
  if (.family_has(object, "gaussian")) {
    sigma_raw <- .extract_param_draw(object, "sigma", draw_index)
    sigma_names <- names(sigma_raw)
    sigma_idx <- as.integer(gsub("sigma\\[|\\]", "", sigma_names))
    sigma[sigma_idx] <- as.numeric(sigma_raw)
  }

  kappa <- vector("list", p)
  if (.family_has(object, "ordinal")) {
    kappa_raw <- .extract_param_draw(object, "kappa", draw_index)
    kappa_names <- names(kappa_raw)
    parts <- strsplit(gsub("kappa\\[|\\]", "", kappa_names), ",")
    j_vals <- as.integer(vapply(parts, `[[`, character(1L), 1L))
    c_vals <- as.integer(vapply(parts, `[[`, character(1L), 2L))
    for (j in unique(j_vals)) {
      mask <- j_vals == j
      kappa[[j]] <- as.numeric(kappa_raw[mask])[order(c_vals[mask])]
    }
  }

  list(sigma = sigma, kappa = kappa)
}


#' Generate a response vector for all p nodes (per-node family dispatch)
#'
#' Used by \code{simulate.bvarnet()} in the forward simulation loop.
#'
#' @param eta Numeric vector of length p. Linear predictor per node.
#' @param family_vec Character vector of length p.
#' @param sigma Numeric vector of length p (NA for non-gaussian).
#' @param kappa List of length p (NULL entries for non-ordinal).
#' @param p Integer. Number of nodes.
#' @return Numeric vector of length p.
#' @noRd
.generate_response_node_vec <- function(eta, family_vec, sigma, kappa, p) {
  y <- numeric(p)
  for (j in seq_len(p)) {
    y[j] <- switch(family_vec[j],
      bernoulli = rbinom(1L, 1L, 1 / (1 + exp(-eta[j]))),
      gaussian  = rnorm(1L, eta[j], sigma[j]),
      ordinal   = {
        kappa_cumsum <- c(0, cumsum(kappa[[j]]))
        C_j <- length(kappa[[j]]) + 1L
        lambda <- (seq_len(C_j) - 1L) * eta[j] - kappa_cumsum
        lambda <- lambda - max(lambda)
        probs  <- exp(lambda) / sum(exp(lambda))
        sample.int(C_j, 1L, prob = probs)
      },
      stop("Unknown family: ", family_vec[j])
    )
  }
  y
}


#' Compute eta matrix using recursive forecasting (subject-by-subject)
#'
#' For each subject: rows up to the conditioning window use observed lag
#' values (from B); subsequent rows use a rolling lag buffer that is updated
#' with predicted values after each step.
#'
#' @param X   Matrix [n_obs, n_fe].
#' @param B   Matrix [n_obs, p*K] (observed lags).
#' @param Z   Matrix [n_obs, n_re] (or 0-col).
#' @param Y_obs Matrix [n_obs, p] (observed outcomes for modeled rows).
#' @param beta Matrix [n_fe, p].
#' @param phi  Matrix [p*K, p].
#' @param sigma Numeric vector length p (gaussian) or NULL.
#' @param kappa List of p ordered vectors (ordinal) or NULL.
#' @param id_char Character vector length n_obs.
#' @param object bvarnet object (for RE extraction).
#' @param subject_re Character.
#' @param new_subject Character.
#' @param draw_index NULL or integer.
#' @param family_vec Character vector of length p (or scalar, recycled).
#' @param cw_by_subject Named integer vector (conditioning window per subject).
#' @param K Integer lag order.
#' @param p Integer number of nodes.
#' @return Matrix [n_obs, p] of linear predictor values.
#' @noRd
.compute_recursive_eta <- function(X, B, Z, Y_obs, beta, phi, sigma, kappa,
                                    id_char, object, subject_re, new_subject,
                                    draw_index, family_vec, cw_by_subject, K, p) {
  n_obs <- nrow(X)
  eta_mat <- matrix(NA_real_, n_obs, p)
  if (length(family_vec) == 1L) family_vec <- rep(family_vec, p)

  # Pre-compute RE values for each node (full vectorised extraction)
  has_re <- !is.null(Z) && ncol(Z) > 0L
  u_list <- vector("list", p)
  if (has_re) {
    for (node in seq_len(p)) {
      u_list[[node]] <- .get_re_for_rows(object, id_char, subject_re,
                                          new_subject, draw_index, node)
    }
  }

  unique_ids <- unique(id_char)
  for (subj_id in unique_ids) {
    rows    <- which(id_char == subj_id)
    T_mod   <- length(rows)
    cw      <- cw_by_subject[subj_id]
    n_cond  <- max(0L, as.integer(cw) - K)

    # Initialise lag buffer from the first modeled row's observed lags
    lag_buffer <- B[rows[1], ]

    for (ti in seq_len(T_mod)) {
      r <- rows[ti]  # global row index

      B_row <- if (ti <= n_cond) B[r, ] else lag_buffer

      for (node in seq_len(p)) {
        u_rows_node <- if (has_re) u_list[[node]][r, , drop = FALSE] else NULL
        eta_val <- .predict_eta_node(
          X[r, , drop = FALSE], matrix(B_row, nrow = 1L),
          beta, phi, Z[r, , drop = FALSE], u_rows_node, node,
          family_node = family_vec[node]
        )
        eta_mat[r, node] <- eta_val
      }

      # Update lag buffer
      if (ti <= n_cond) {
        lag_buffer <- .update_lag_buffer(lag_buffer, Y_obs[r, ], K, p)
      } else {
        y_pred <- .recursive_lag_value(eta_mat[r, ], family_vec, sigma, kappa)
        lag_buffer <- .update_lag_buffer(lag_buffer, y_pred, K, p)
      }
    }
  }

  colnames(eta_mat) <- colnames(object$standata$Y)
  eta_mat
}


# ═══════════════════════════════════════════════════════════════════════════════
#                       predict.bvarnet
# ═══════════════════════════════════════════════════════════════════════════════

#' Predict from a fitted bvarnet model
#'
#' Computes one-step-ahead or recursive forecasts for long-format time-series
#' data.
#' Supports population-level (\code{subject_re = "zero"}) and
#' subject-specific (\code{subject_re = "posterior-mean"}) predictions.
#' Also serves as the out-of-sample engine: fit on training data, call
#' \code{predict(fit, newdata = test_data)}.
#'
#' @param object A \code{bvarnet} object from \code{bvar()}.
#' @param newdata Data frame in long format. If \code{NULL}, the original
#'   training data design matrices (stored in \code{object$standata}) are used
#'   for in-sample fitted values.
#' @param forecast Character. \code{"one-step"} (default) uses observed lag
#'   values from \code{newdata} for every row. \code{"recursive"} feeds
#'   predicted values back into the lag buffer after the conditioning window.
#' @param conditioning_window Integer scalar, named integer vector, or
#'   \code{NULL}. Number of observed time points (per subject) used to
#'   initialise the lag buffer before recursive forecasting begins.
#'   Must be \code{>= K}. If \code{NULL}, defaults to \code{K} (minimum lag
#'   history). Ignored when \code{forecast = "one-step"}.
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
#' @return For \code{type = "link"} or \code{"response"}: a data.frame with
#'   columns named after \code{id_col} and \code{time_col}, plus
#'   \code{predicted_<y>} for each outcome variable. Only the modeled
#'   observations are included (no NA padding).
#'   For \code{type = "probabilities"}: a named list of data.frames
#'   (one per outcome), each with id, time, and probability columns.
#'   When \code{method = "posterior-sample"}, additional \code{_sd} columns
#'   are included and \code{attr(,"ndraws")} records the number of draws used.
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
                            forecast    = c("one-step", "recursive"),
                            conditioning_window = NULL,
                            type        = c("link", "response", "probabilities"),
                            method      = c("posterior-mean", "posterior-sample"),
                            ndraws      = NULL,
                            seed        = NULL,
                            subject_re  = c("zero", "posterior-mean", "sample"),
                            new_subject = c("zero", "sample"),
                            ...) {

  forecast    <- match.arg(forecast)
  type        <- match.arg(type)
  method      <- match.arg(method)
  subject_re  <- match.arg(subject_re)
  new_subject <- match.arg(new_subject)

  sd     <- object$standata
  family_vec <- object$family
  p      <- sd$p
  n_re   <- sd$n_re

  # Degrade subject_re to "zero" when no RE

  if (n_re == 0L) subject_re <- "zero"

  # --- build design matrices / use stored ---
  if (is.null(newdata)) {
    # In-sample: use stored matrices
    X        <- sd$X
    B        <- sd$B
    Z        <- sd$Z
    Y_obs    <- sd$Y
    id_char  <- as.character(sd$id_levels[sd$id])
    time_vec <- sd$time_obs
    orig_order <- NULL
  } else {
    pred_sd  <- .build_pred_standata(newdata, object)
    X        <- pred_sd$X
    B        <- pred_sd$B
    Z        <- pred_sd$Z
    Y_obs    <- pred_sd$Y_obs
    id_char  <- pred_sd$id_char
    time_vec <- pred_sd$time_obs
    # row_map holds original newdata row indices; use to restore input order
    orig_order <- order(pred_sd$row_map)
  }

  n_obs <- nrow(X)
  K     <- sd$K

  # --- resolve conditioning window for recursive mode ---------------------
  cw_by_subject <- NULL
  if (forecast == "recursive") {
    unique_ids    <- unique(id_char)
    cw_by_subject <- .resolve_conditioning_window_by_subject(
      conditioning_window, unique_ids, K
    )
  }

  # --- dispatch by method ------------------------------------------------
  if (method == "posterior-mean") {
    # Single deterministic prediction
    beta  <- if (sd$n_fe > 0L) .reshape_beta(.extract_param_draw(object, "beta", NULL), sd) else matrix(0, 0L, sd$p)
    phi   <- .reshape_phi(.extract_param_draw(object, "phi", NULL), sd)
    sk    <- .extract_sigma_kappa(object, NULL)
    sigma <- sk$sigma
    kappa <- sk$kappa

    eta_mat <- matrix(NA_real_, n_obs, p)
    if (forecast == "recursive") {
      eta_mat <- .compute_recursive_eta(
        X, B, Z, Y_obs, beta, phi, sigma, kappa,
        id_char, object, subject_re, new_subject,
        NULL, family_vec, cw_by_subject, K, p
      )
    } else {
      for (node in seq_len(p)) {
        u_rows <- .get_re_for_rows(object, id_char, subject_re,
                                    new_subject, NULL, node)
        eta_mat[, node] <- .predict_eta_node(X, B, beta, phi, Z, u_rows, node,
                                              family_node = family_vec[node])
      }
    }
    colnames(eta_mat) <- colnames(sd$Y)

    result <- .eta_to_output(eta_mat, family_vec, type, sigma, kappa)

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
      # probabilities: per-node accumulators with family-specific ncol
      sk_mean <- .extract_sigma_kappa(object, NULL)
      sum_out <- vector("list", p)
      sum_sq  <- vector("list", p)
      for (j in seq_len(p)) {
        nc <- switch(family_vec[j],
          bernoulli = 1L,
          gaussian  = 2L,
          ordinal   = length(sk_mean$kappa[[j]]) + 1L
        )
        sum_out[[j]] <- matrix(0, n_obs, nc)
        sum_sq[[j]]  <- matrix(0, n_obs, nc)
      }
    }

    for (s in draw_idx) {
      beta_s  <- if (sd$n_fe > 0L) .reshape_beta(.extract_param_draw(object, "beta", s), sd) else matrix(0, 0L, sd$p)
      phi_s   <- .reshape_phi(.extract_param_draw(object, "phi", s), sd)
      sk_s    <- .extract_sigma_kappa(object, s)
      sigma_s <- sk_s$sigma
      kappa_s <- sk_s$kappa

      eta_s <- matrix(NA_real_, n_obs, p)
      if (forecast == "recursive") {
        eta_s <- .compute_recursive_eta(
          X, B, Z, Y_obs, beta_s, phi_s, sigma_s, kappa_s,
          id_char, object, subject_re, new_subject,
          s, family_vec, cw_by_subject, K, p
        )
      } else {
        for (node in seq_len(p)) {
          u_rows <- .get_re_for_rows(object, id_char, subject_re,
                                      new_subject, s, node)
          eta_s[, node] <- .predict_eta_node(X, B, beta_s, phi_s, Z, u_rows, node,
                                              family_node = family_vec[node])
        }
      }
      colnames(eta_s) <- colnames(sd$Y)

      out_s <- .eta_to_output(eta_s, family_vec, type, sigma_s, kappa_s)

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

  # --- wrap into data.frame output ---
  id_col_name   <- sd$design_spec$id_col
  time_col_name <- sd$design_spec$time_col
  y_names       <- colnames(sd$Y)

  if (type %in% c("link", "response")) {
    df <- data.frame(id_char, time_vec, stringsAsFactors = FALSE)
    names(df) <- c(id_col_name, time_col_name)

    pred_cols <- as.data.frame(result)
    names(pred_cols) <- paste0("predicted_", y_names)
    df <- cbind(df, pred_cols)

    if (method == "posterior-sample") {
      sd_mat <- attr(result, "sd")
      sd_cols <- as.data.frame(sd_mat)
      names(sd_cols) <- paste0("predicted_", y_names, "_sd")
      df <- cbind(df, sd_cols)
      attr(df, "ndraws") <- attr(result, "ndraws")
    }

    # Restore original newdata row order when applicable
    if (!is.null(orig_order)) df <- df[orig_order, , drop = FALSE]
    rownames(df) <- NULL
    result <- df
  } else {
    # type == "probabilities": list of data.frames, one per node
    sd_attr <- attr(result, "sd")
    nd_attr <- attr(result, "ndraws")
    result_list <- vector("list", p)
    for (j in seq_len(p)) {
      dj <- data.frame(id_char, time_vec, stringsAsFactors = FALSE)
      names(dj) <- c(id_col_name, time_col_name)
      prob_df <- as.data.frame(result[[j]])
      dj <- cbind(dj, prob_df)
      if (method == "posterior-sample" && !is.null(sd_attr)) {
        sd_j <- as.data.frame(sd_attr[[j]])
        names(sd_j) <- paste0(names(prob_df), "_sd")
        dj <- cbind(dj, sd_j)
      }
      if (!is.null(orig_order)) dj <- dj[orig_order, , drop = FALSE]
      rownames(dj) <- NULL
      result_list[[j]] <- dj
    }
    names(result_list) <- y_names
    if (method == "posterior-sample" && !is.null(nd_attr))
      attr(result_list, "ndraws") <- nd_attr
    result <- result_list
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
  family_vec <- object$family
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
    beta  <- if (n_fe > 0L) .reshape_beta(.extract_param_draw(object, "beta", draw_index), sd) else matrix(0, 0L, p)
    phi   <- .reshape_phi(.extract_param_draw(object, "phi", draw_index), sd)
    sk    <- .extract_sigma_kappa(object, draw_index)
    sigma <- sk$sigma
    kappa <- sk$kappa
    sd_u  <- if (n_re > 0L) .reshape_sd_u(
      .extract_param_draw(object, "sd_u", draw_index), sd) else NULL

    # Extract intercept and covariate effects from beta per-node
    # For mixed families, beta always has intercept row (from combined standata).
    # For homogeneous ordinal, to_stan_data strips intercept → all rows are covariates.
    has_intercept_row <- any(family_vec != "ordinal") || n_fe > q
    alpha <- numeric(p)
    if (has_intercept_row) {
      for (j in seq_len(p)) {
        if (family_vec[j] != "ordinal") {
          alpha[j] <- beta[1, j]  # intercept row
        }
        # ordinal: alpha stays 0 (absorbed in kappa)
      }
      gamma <- if (n_fe > 1L) beta[2:n_fe, , drop = FALSE] else matrix(0, 0, p)
    } else {
      # Homogeneous ordinal: no intercept in beta; all rows are covariates
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
    Y_full <- array(NA_real_, dim = c(N, T_total, p))

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
        Y_full[i, t, ] <- .generate_response_node_vec(
          eta, family_vec, sigma, kappa, p)
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

        Y_full[i, t, ] <- .generate_response_node_vec(
          eta, family_vec, sigma, kappa, p)
      }
    }

    # Discard burnin
    keep_idx <- (burnin + 1L):T_total
    Y_keep <- Y_full[, keep_idx, , drop = FALSE]
    X_keep <- X_cov[, keep_idx, , drop = FALSE]

    # Assemble long-format data frame (preserve original y_cols names if available)
    y_cols_out <- if (!is.null(sd$design_spec$y_cols) &&
                      length(sd$design_spec$y_cols) == p) {
      sd$design_spec$y_cols
    } else {
      paste0("y_", seq_len(p))
    }
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
