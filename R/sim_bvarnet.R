# ──────────────────────────────────────────────────────────────────────────────
# sim_var() — simulate data from a multilevel VAR model
#
# Generates data matching the exact generative model estimated by
# bvar() for all three families: bernoulli, ordinal, gaussian.
# ──────────────────────────────────────────────────────────────────────────────

#' Simulate data from a multilevel VAR model
#'
#' Generates data from the generative model implied by each Stan model family.
#' Useful for testing parameter recovery and model validation.
#'
#' @param N Integer. Number of subjects (groups).
#' @param T_obs Integer. Number of time points per subject.
#' @param p Integer. Number of outcome nodes.
#' @param K Integer. AR order (default 1).
#' @param family Character. One of \code{"bernoulli"}, \code{"ordinal"},
#'   \code{"gaussian"}.
#' @param alpha Numeric vector of length \code{p}. Population intercepts
#'   (on logit scale for bernoulli, identity for gaussian). For ordinal,
#'   this is absorbed into kappa and should be left NULL. Generated if NULL.
#' @param gamma Matrix \code{q x p}. Population covariate effects. Generated
#'   if NULL and \code{q > 0}.
#' @param Phi Matrix \code{(p*K) x p}. Population lag coefficients. Generated
#'   if NULL.
#' @param sigma Numeric vector of length \code{p}. Residual SD per node
#'   (gaussian only). Generated if NULL.
#' @param kappa List of \code{p} ordered vectors, each of length \code{C-1}.
#'   Cutpoints per node (ordinal only). Generated if NULL.
#' @param q Integer. Number of covariates (default 0).
#' @param x_gen Function \code{f(N, T_obs)} returning an \code{N x T_obs x q}
#'   array of covariates. If NULL, default generation is used.
#' @param sd_alpha Numeric. SD of random intercepts (scalar or p-vector).
#'   Default 0.5. Set to 0 to simulate a fixed-effects-only model with no
#'   between-person variation in intercepts.
#' @param sd_phi Numeric. SD of random lag coefficients (scalar or matrix).
#'   Default 0.2.
#' @param sd_gamma Numeric or NULL. SD of random covariate slopes. NULL means
#'   no random slopes on covariates.
#' @param re_temporal Logical. Include random slopes on lag predictors?
#'   Default FALSE.
#' @param C Integer. Number of ordinal categories (ordinal only, default 5).
#' @param burnin Integer. Number of time points to discard as warmup before
#'   recording data (default 500). The VAR process is simulated for
#'   \code{burnin + T_obs} time points per subject, and the first
#'   \code{burnin} are discarded. This allows the process to reach its
#'   stationary distribution before data collection begins.
#' @param seed Integer or NULL. RNG seed.
#'
#' @details
#' To simulate a VAR without any random effects (i.e. all subjects share
#' identical parameters), set \code{sd_alpha = 0}, \code{re_temporal = FALSE}
#' (the default), and \code{sd_gamma = NULL} (the default).
#'
#' @return A list with two components:
#' \describe{
#'   \item{data}{A long-format data frame with columns \code{id}, \code{t},
#'     \code{y_1}, ..., \code{y_p}, and optionally \code{x_1}, ..., \code{x_q}.}
#'   \item{truth}{A list of true generating parameters.}
#' }
#'
#' @export
sim_var <- function(
    N,
    T_obs,
    p,
    K        = 1L,
    family   = c("bernoulli", "ordinal", "gaussian"),
    alpha    = NULL,
    gamma    = NULL,
    Phi      = NULL,
    sigma    = NULL,
    kappa    = NULL,
    q        = 0L,
    x_gen    = NULL,
    sd_alpha = 0.5,
    sd_phi   = 0.2,
    sd_gamma = NULL,
    re_temporal = FALSE,
    C        = 5L,
    burnin   = 500L,
    seed     = NULL
) {
  ## ── 0. Validation & setup ──────────────────────────────────────────────────
  family <- match.arg(family)
  N      <- as.integer(N)
  T_obs  <- as.integer(T_obs)
  p      <- as.integer(p)
  K      <- as.integer(K)
  q      <- as.integer(q)
  C      <- as.integer(C)

  burnin <- as.integer(burnin)

  stopifnot(N >= 1L, T_obs >= K + 1L, p >= 1L, K >= 1L, q >= 0L, burnin >= 0L)

  if (family == "ordinal") stopifnot(C >= 2L)
  if (!is.null(seed)) set.seed(seed)

  PK <- p * K
  T_total <- T_obs + burnin  # total simulated length (burnin then recorded)

  # sd_alpha: broadcast scalar to p-vector
  if (length(sd_alpha) == 1L) sd_alpha <- rep(sd_alpha, p)
  stopifnot(length(sd_alpha) == p)

  ## ── 1. Covariates ─────────────────────────────────────────────────────────
  # Generate for the full T_total length (burnin + recorded)
  if (q > 0L) {
    if (!is.null(x_gen)) {
      X_cov <- x_gen(N, T_total)
      stopifnot(is.array(X_cov), identical(dim(X_cov), c(N, T_total, q)))
    } else {
      X_cov <- generate_default_covariates(N, T_total, q)
    }
  } else {
    X_cov <- array(0, dim = c(N, T_total, 0L))
  }

  ## ── 2. Population parameters ──────────────────────────────────────────────
  # Intercepts (not used for ordinal — kappa absorbs intercept)
  if (is.null(alpha)) {
    if (family == "ordinal") {
      alpha <- rep(0, p)  # no intercept for ordinal
    } else {
      alpha <- runif(p, -1, 1)
    }
  }
  stopifnot(length(alpha) == p)

  # Covariate effects
  if (q > 0L) {
    if (is.null(gamma)) {
      gamma <- matrix(runif(q * p, -0.5, 0.5), nrow = q, ncol = p)
    }
    stopifnot(is.matrix(gamma), nrow(gamma) == q, ncol(gamma) == p)
  } else {
    gamma <- matrix(0, nrow = 0L, ncol = p)
  }

  # Lag coefficients
  if (is.null(Phi)) {
    Phi <- generate_stable_phi(p, K)
  }
  stopifnot(is.matrix(Phi), nrow(Phi) == PK, ncol(Phi) == p)
  if (!check_var_stability(Phi, p, K)) {
    warning("Provided Phi is not VAR-stable (eigenvalues >= 1). ",
            "Rescaling to ensure stability.")
    Phi <- rescale_to_stable(Phi, p, K)
  }

  # Gaussian: residual SD
  if (family == "gaussian") {
    if (is.null(sigma)) sigma <- runif(p, 0.5, 1.5)
    stopifnot(length(sigma) == p, all(sigma > 0))
  }

  # Ordinal: cutpoints
  if (family == "ordinal") {
    if (is.null(kappa)) {
      kappa <- replicate(p, generate_default_kappa(C), simplify = FALSE)
    }
    stopifnot(is.list(kappa), length(kappa) == p)
    for (node in seq_len(p)) {
      stopifnot(length(kappa[[node]]) == C - 1L)
      stopifnot(!is.unsorted(kappa[[node]], strictly = TRUE))
    }
  }

  ## ── 3. Person-level random effects ────────────────────────────────────────
  # Random intercepts (always present for bernoulli/gaussian, not for ordinal)
  alpha_i <- matrix(NA_real_, N, p)
  for (i in seq_len(N)) {
    alpha_i[i, ] <- rnorm(p, mean = alpha, sd = sd_alpha)
  }

  # Random temporal effects
  Phi_i <- array(NA_real_, dim = c(N, PK, p))
  if (isTRUE(re_temporal)) {
    sd_phi_mat <- if (is.matrix(sd_phi)) sd_phi else matrix(sd_phi, PK, p)
    for (i in seq_len(N)) {
      Phi_i[i, , ] <- Phi + matrix(rnorm(PK * p, 0, sd_phi_mat), PK, p)
    }
  } else {
    for (i in seq_len(N)) {
      Phi_i[i, , ] <- Phi
    }
  }

  # Random covariate slopes
  gamma_i <- array(NA_real_, dim = c(N, max(q, 0L), p))
  if (!is.null(sd_gamma) && q > 0L) {
    sd_gamma_mat <- if (is.matrix(sd_gamma)) sd_gamma else matrix(sd_gamma, q, p)
    for (i in seq_len(N)) {
      gamma_i[i, , ] <- gamma + matrix(rnorm(q * p, 0, sd_gamma_mat), q, p)
    }
  } else if (q > 0L) {
    for (i in seq_len(N)) {
      gamma_i[i, , ] <- gamma
    }
  }

  ## ── 4. Forward simulation ─────────────────────────────────────────────────
  # Simulate T_total = burnin + T_obs time points; discard first `burnin`
  if (family == "gaussian") {
    Y_full <- array(NA_real_, dim = c(N, T_total, p))
  } else if (family == "ordinal") {
    Y_full <- array(NA_integer_, dim = c(N, T_total, p))
  } else {
    Y_full <- array(0L, dim = c(N, T_total, p))
  }

  for (i in seq_len(N)) {
    # --- Initialize first K time points (intercept + covariates only) ---
    for (t in seq_len(K)) {
      eta <- alpha_i[i, ]
      if (q > 0L) {
        g_mat <- matrix(gamma_i[i, , ], nrow = q, ncol = p)
        x_vec <- matrix(X_cov[i, t, ], nrow = 1L)
        eta <- eta + as.numeric(x_vec %*% g_mat)
      }
      Y_full[i, t, ] <- generate_response(eta, family, sigma, kappa, C, p)
    }

    # --- Forward simulate t = K+1 .. T_total ---
    for (t in (K + 1L):T_total) {
      # Build lag vector: [y_{t-1}, y_{t-2}, ..., y_{t-K}]
      lag_y <- numeric(PK)
      for (lag in seq_len(K)) {
        idx <- ((lag - 1L) * p + 1L):(lag * p)
        lag_y[idx] <- Y_full[i, t - lag, ]
      }

      # Linear predictor
      eta <- alpha_i[i, ] +
        as.numeric(t(Phi_i[i, , ]) %*% lag_y)

      if (q > 0L) {
        g_mat <- matrix(gamma_i[i, , ], nrow = q, ncol = p)
        x_vec <- matrix(X_cov[i, t, ], nrow = 1L)
        eta <- eta + as.numeric(x_vec %*% g_mat)
      }

      Y_full[i, t, ] <- generate_response(eta, family, sigma, kappa, C, p)
    }
  }

  # Discard burnin period
  keep_idx <- (burnin + 1L):T_total
  Y     <- Y_full[, keep_idx, , drop = FALSE]
  X_cov <- X_cov[, keep_idx, , drop = FALSE]

  ## ── 5. Assemble long-format data frame ────────────────────────────────────
  df <- assemble_long_df(Y, X_cov, N, T_obs, p, q)

  ## ── 6. Build truth object ─────────────────────────────────────────────────
  # Build sd_u to match Stan layout: matrix[p, n_re]
  #   n_re columns depend on re_temporal and sd_gamma
  sd_u_list <- list()
  if (!is.null(sd_gamma) && q > 0L) {
    sd_gamma_mat <- if (is.matrix(sd_gamma)) sd_gamma else matrix(sd_gamma, q, p)
    sd_u_list[["gamma"]] <- t(sd_gamma_mat)  # p x q
  }
  if (isTRUE(re_temporal)) {
    sd_phi_mat <- if (is.matrix(sd_phi)) sd_phi else matrix(sd_phi, PK, p)
    sd_u_list[["phi"]] <- t(sd_phi_mat)  # p x PK
  }
  if (length(sd_u_list) > 0L) {
    sd_u <- do.call(cbind, sd_u_list)
  } else {
    sd_u <- matrix(0, p, 0L)
  }

  truth <- list(
    alpha      = alpha,
    gamma      = if (q > 0L) gamma else NULL,
    Phi        = Phi,
    sigma      = if (family == "gaussian") sigma else NULL,
    kappa      = if (family == "ordinal") kappa else NULL,
    sd_alpha   = sd_alpha,
    sd_u       = sd_u,
    alpha_i    = alpha_i,
    Phi_i      = Phi_i,
    gamma_i    = if (q > 0L) gamma_i else NULL,
    family     = family,
    N          = N,
    T_obs      = T_obs,
    p          = p,
    K          = K,
    q          = q,
    C          = if (family == "ordinal") C else NULL,
    burnin     = burnin
  )

  list(data = df, truth = truth)
}


# ──────────────────────────────────────────────────────────────────────────────
# Internal helper functions
# ──────────────────────────────────────────────────────────────────────────────

#' Generate default covariates (odd = continuous, even = binary)
#' @noRd
generate_default_covariates <- function(N, T_obs, q) {
  X <- array(NA_real_, dim = c(N, T_obs, q))
  for (j in seq_len(q)) {
    if (j %% 2L == 1L) {
      X[, , j] <- rnorm(N * T_obs)
    } else {
      X[, , j] <- rbinom(N * T_obs, 1L, 0.5)
    }
  }
  X
}


#' Generate a stable VAR(K) lag coefficient matrix
#'
#' Creates a Phi matrix with moderate auto-regressive effects on the diagonal
#' and sparse cross-lag effects, rescaled to ensure VAR stability.
#' @noRd
generate_stable_phi <- function(p, K) {
  PK <- p * K
  Phi <- matrix(0, PK, p)

  for (lag in seq_len(K)) {
    block_rows <- ((lag - 1L) * p + 1L):(lag * p)
    block <- matrix(0, p, p)

    # Autoregressive effects (decay with lag)
    diag(block) <- runif(p, 0.2, 0.6) / lag

    # Sparse cross-lag effects (about 30% density)
    for (row in seq_len(p)) {
      for (col in seq_len(p)) {
        if (row != col && runif(1) < 0.3) {
          block[row, col] <- runif(1, -0.2, 0.2) / lag
        }
      }
    }

    Phi[block_rows, ] <- block
  }

  # Ensure stability
  if (!check_var_stability(Phi, p, K)) {
    Phi <- rescale_to_stable(Phi, p, K)
  }

  Phi
}


#' Check VAR stability (all companion matrix eigenvalues inside unit circle)
#' @noRd
check_var_stability <- function(Phi, p, K) {
  companion <- build_companion(Phi, p, K)
  max_ev <- max(abs(eigen(companion, only.values = TRUE)$values))
  max_ev < 1.0
}


#' Build companion matrix for VAR(K)
#' @noRd
build_companion <- function(Phi, p, K) {
  PK <- p * K
  if (K == 1L) {
    return(t(Phi))  # p x p
  }

  companion <- matrix(0, PK, PK)
  for (lag in seq_len(K)) {
    rows <- ((lag - 1L) * p + 1L):(lag * p)
    companion[1:p, rows] <- t(Phi[rows, ])
  }
  if (K > 1L) {
    companion[(p + 1L):PK, 1:(p * (K - 1L))] <- diag(p * (K - 1L))
  }
  companion
}


#' Rescale Phi to ensure VAR stability (spectral radius < 0.95)
#' @noRd
rescale_to_stable <- function(Phi, p, K, target = 0.95) {
  companion <- build_companion(Phi, p, K)
  max_ev <- max(abs(eigen(companion, only.values = TRUE)$values))
  if (max_ev >= 1.0) {
    scale_factor <- target / max_ev
    Phi <- Phi * scale_factor
  }
  Phi
}


#' Generate default ordered cutpoints for ordinal model
#' @noRd
generate_default_kappa <- function(C) {
  # C-1 cutpoints evenly spaced from approximately -2 to 2
  seq(-2, 2, length.out = C - 1L)
}


#' Generate response for a single time point (all p nodes)
#'
#' @param eta Numeric vector of length p. Linear predictor per node.
#' @param family Character. Model family.
#' @param sigma Numeric vector (gaussian only).
#' @param kappa List of p vectors (ordinal only).
#' @param C Integer (ordinal only).
#' @param p Integer. Number of nodes.
#' @return Vector of length p with generated responses.
#' @noRd
generate_response <- function(eta, family, sigma = NULL, kappa = NULL,
                              C = NULL, p) {
  switch(family,
    bernoulli = generate_response_binary(eta),
    gaussian  = generate_response_gaussian(eta, sigma),
    ordinal   = generate_response_ordinal(eta, kappa, C, p),
    stop("Unknown family: ", family)
  )
}


#' Binary response: logistic link + Bernoulli draw
#' @noRd
generate_response_binary <- function(eta) {
  prob <- 1 / (1 + exp(-eta))
  rbinom(length(eta), 1L, prob)
}


#' Gaussian response: identity link + normal noise
#' @noRd
generate_response_gaussian <- function(eta, sigma) {
  rnorm(length(eta), mean = eta, sd = sigma)
}


#' Ordinal response: adjacent-category logit
#'
#' Matches the Stan parameterisation:
#'   lambda[c] = (c-1) * eta - kappa_cumsum[c]
#'   P(Y = c) = exp(lambda[c]) / sum(exp(lambda))
#' @noRd
generate_response_ordinal <- function(eta, kappa, C, p) {
  y <- integer(p)
  for (node in seq_len(p)) {
    # Cumulative sum of kappa for this node
    kappa_cumsum <- c(0, cumsum(kappa[[node]]))

    # Adjacent-category log-odds
    lambda <- (seq_len(C) - 1L) * eta[node] - kappa_cumsum

    # Softmax with log-sum-exp stability
    lambda <- lambda - max(lambda)
    probs <- exp(lambda) / sum(exp(lambda))

    y[node] <- sample.int(C, size = 1L, prob = probs)
  }
  y
}


#' Assemble Y and X arrays into a long-format data frame
#' @noRd
assemble_long_df <- function(Y, X_cov, N, T_obs, p, q) {
  rows <- vector("list", N)
  for (i in seq_len(N)) {
    row_data <- data.frame(
      id = rep(i, T_obs),
      t  = seq_len(T_obs)
    )
    # y columns
    for (j in seq_len(p)) {
      row_data[[paste0("y_", j)]] <- Y[i, , j]
    }
    # x columns
    if (q > 0L) {
      for (j in seq_len(q)) {
        row_data[[paste0("x_", j)]] <- X_cov[i, , j]
      }
    }
    rows[[i]] <- row_data
  }
  do.call(rbind, rows)
}


# ──────────────────────────────────────────────────────────────────────────────
# compare_to_truth() — compare fitted parameters to generating truth
# ──────────────────────────────────────────────────────────────────────────────

#' Compare fitted model parameters to simulation truth
#'
#' Extracts posterior summaries from a fitted \code{bvarnet} object and
#' compares them to the true parameter values used for data generation.
#'
#' @param fit A fitted \code{bvarnet} object (output from \code{bvar()}).
#' @param truth The \code{truth} component from \code{sim_var()} output.
#' @param ci_width Numeric. Width of the credible interval (default 0.90).
#'
#' @return A data frame with columns: parameter, node, index, true_value,
#'   post_mean, post_sd, ci_lower, ci_upper, covered (logical).
#'
#' @export
compare_to_truth <- function(fit, truth, ci_width = 0.90) {
  stopifnot(inherits(fit, "bvarnet"))

  alpha_lo <- (1 - ci_width) / 2
  alpha_hi <- 1 - alpha_lo
  p <- truth$p
  family <- truth$family

  results <- list()

  # ── beta (intercept + covariates) ────────────────────────────────────────
  draws_beta <- extract_draws(fit, "beta")
  # Stan layout: beta[fe_idx, node]
  n_fe <- fit$standata$n_fe
  for (node in seq_len(p)) {
    for (fe in seq_len(n_fe)) {
      par_name <- paste0("beta[", fe, ",", node, "]")
      d <- draws_beta[, par_name]

      # Map to truth:
      #   fe=1 for bernoulli/gaussian is intercept → alpha[node]
      #   fe=1 for ordinal is first covariate (no intercept)
      #   fe>1 is gamma[fe-1, node] (bernoulli/gaussian) or gamma[fe, node] (ordinal)
      if (family %in% c("bernoulli", "gaussian") && fe == 1L) {
        true_val <- truth$alpha[node]
        par_label <- "intercept"
      } else {
        gamma_idx <- if (family == "ordinal") fe else fe - 1L
        true_val <- if (!is.null(truth$gamma) && gamma_idx <= nrow(truth$gamma)) {
          truth$gamma[gamma_idx, node]
        } else {
          NA_real_
        }
        par_label <- paste0("gamma_", gamma_idx)
      }

      results[[length(results) + 1L]] <- data.frame(
        parameter = par_label,
        node      = node,
        true_value = true_val,
        post_mean  = mean(d),
        post_sd    = sd(d),
        ci_lower   = unname(quantile(d, alpha_lo)),
        ci_upper   = unname(quantile(d, alpha_hi)),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── phi (lag coefficients) ──────────────────────────────────────────────
  draws_phi <- extract_draws(fit, "phi")
  PK <- p * truth$K
  for (node in seq_len(p)) {
    for (lag_idx in seq_len(PK)) {
      par_name <- paste0("phi[", lag_idx, ",", node, "]")
      d <- draws_phi[, par_name]
      true_val <- truth$Phi[lag_idx, node]

      results[[length(results) + 1L]] <- data.frame(
        parameter  = "phi",
        node       = node,
        true_value = true_val,
        post_mean  = mean(d),
        post_sd    = sd(d),
        ci_lower   = unname(quantile(d, alpha_lo)),
        ci_upper   = unname(quantile(d, alpha_hi)),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── sigma (gaussian only) ──────────────────────────────────────────────
  if (family == "gaussian") {
    draws_sigma <- extract_draws(fit, "sigma")
    for (node in seq_len(p)) {
      par_name <- paste0("sigma[", node, "]")
      d <- draws_sigma[, par_name]

      results[[length(results) + 1L]] <- data.frame(
        parameter  = "sigma",
        node       = node,
        true_value = truth$sigma[node],
        post_mean  = mean(d),
        post_sd    = sd(d),
        ci_lower   = unname(quantile(d, alpha_lo)),
        ci_upper   = unname(quantile(d, alpha_hi)),
        stringsAsFactors = FALSE
      )
    }
  }

  # ── kappa (ordinal only) ───────────────────────────────────────────────
  if (family == "ordinal") {
    draws_kappa <- extract_draws(fit, "kappa")
    C <- truth$C
    for (node in seq_len(p)) {
      for (k in seq_len(C - 1L)) {
        par_name <- paste0("kappa[", node, ",", k, "]")
        d <- draws_kappa[, par_name]

        results[[length(results) + 1L]] <- data.frame(
          parameter  = "kappa",
          node       = node,
          true_value = truth$kappa[[node]][k],
          post_mean  = mean(d),
          post_sd    = sd(d),
          ci_lower   = unname(quantile(d, alpha_lo)),
          ci_upper   = unname(quantile(d, alpha_hi)),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  # ── Combine and add coverage indicator ─────────────────────────────────
  out <- do.call(rbind, results)
  out$covered <- out$true_value >= out$ci_lower & out$true_value <= out$ci_upper
  rownames(out) <- NULL
  out
}
