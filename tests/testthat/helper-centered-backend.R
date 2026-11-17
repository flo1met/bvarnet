# Independent natural-coordinate density, including RE and truncation constants.
# Keep this in terms of the reference model, not the centered Stan likelihood.
centered_fixture <- function(family, n_re = 2L, prior_family = 1L, empty_x = FALSE) {
  x <- c(8, 12, 9, 11, 10, 10)
  w <- c(-1, 1, -2, 2, 0, 3)
  X <- cbind(x, w, x * w)
  if (empty_x) X <- X[, FALSE, drop = FALSE]
  if (family != "ordinal") X <- cbind(1, X)
  B <- matrix(seq_len(24) / 10, 6, 4)
  if (family == "bernoulli") B <- matrix(as.numeric(seq_len(24) %% 3 == 0), 6, 4)
  if (family == "ordinal") B <- matrix(rep(1:4, 6), 6, 4)
  B[2, ] <- 0 # genuine observed zeros need the same centering as other rows
  Z <- if (n_re > 0) cbind(x, B[, 1]) else matrix(numeric(), 6, 0)
  Y <- switch(family, gaussian = matrix(seq_len(12) / 3, 6, 2),
              bernoulli = matrix(rep(c(0L, 1L), 6), 6, 2),
              ordinal = matrix(rep(1:4, 3), 6, 2))
  sd <- list(p = 2L, K = 2L, J = 2L, n_obs = 6L, n_fe = ncol(X),
             n_re = ncol(Z), id = c(1L, 1L, 1L, 2L, 2L, 2L),
             Y = Y, X = X, B = B, Z = Z)
  for (block in c("beta", "phi", "sd", if (family == "ordinal") "kappa" else "intercept",
                  if (family == "gaussian") "sigma")) {
    sd[[paste0("prior_", block, "_fam")]] <- prior_family
    sd[[paste0(block, "_loc")]] <- if (block %in% c("sd", "sigma")) .3 else -.2
    sd[[paste0(block, "_scale")]] <- if (block == "phi") .6 else 1.3
    sd[[paste0(block, "_df")]] <- 5
  }
  if (family == "ordinal") sd$C <- 4L
  pars <- list(beta = matrix(seq_len(sd$n_fe * 2) / 30, sd$n_fe, 2),
               phi = matrix(c(.12, -.14, .08, -.03, -.1, .05, .02, .09), 4, 2),
               sd_u = matrix(seq_len(2 * sd$n_re) / 10 + .2, 2, sd$n_re),
               z_u = array(seq_len(2 * 2 * sd$n_re) / 20 - .2, c(2, 2, sd$n_re)))
  if (family == "gaussian") pars$sigma <- c(.8, 1.3)
  if (family == "ordinal") pars$kappa <- rbind(c(1.6, .2, 1.1), c(.9, -.4, 1.2))
  list(data = sd, pars = pars)
}

reference_density <- function(pars, sd, family) {
  lp <- function(x, block, half = FALSE) {
    code <- sd[[paste0("prior_", block, "_fam")]]
    loc <- sd[[paste0(block, "_loc")]]
    scale <- sd[[paste0(block, "_scale")]]
    df <- sd[[paste0(block, "_df")]]
    value <- oracle_lprior(x, code, loc, scale, df)
    if (half) {
      log_mass <- switch(as.character(code),
        "1" = pnorm(0, loc, scale, lower.tail = FALSE, log.p = TRUE),
        "2" = pt(-loc / scale, df, lower.tail = FALSE, log.p = TRUE),
        "3" = pcauchy(0, loc, scale, lower.tail = FALSE, log.p = TRUE))
      value <- value - length(x) * log_mass
    }
    value
  }
  eta <- sd$X %*% pars$beta + sd$B %*% pars$phi
  if (sd$n_re > 0) for (node in seq_len(sd$p)) {
    z <- matrix(pars$z_u[node, , ], sd$J, sd$n_re)
    u <- sweep(z, 2, pars$sd_u[node, ], "*")
    eta[, node] <- eta[, node] + rowSums(sd$Z * u[sd$id, , drop = FALSE])
  }
  prior <- lp(pars$phi, "phi") + lp(pars$sd_u, "sd", half = TRUE) +
    sum(dnorm(pars$z_u, log = TRUE))
  if (family == "ordinal") {
    prior <- prior + lp(pars$beta, "beta") + lp(pars$kappa, "kappa")
    likelihood <- sum(vapply(seq_len(sd$p), function(j)
      oracle_ll_acat_node(eta[, j], pars$kappa[j, ], sd$Y[, j], sd$C), numeric(1)))
  } else {
    prior <- prior + lp(pars$beta[1, ], "intercept") + lp(pars$beta[-1, ], "beta")
    likelihood <- if (family == "bernoulli")
      sum(dbinom(sd$Y, 1, plogis(eta), log = TRUE)) else {
        prior <- prior + lp(pars$sigma, "sigma", half = TRUE)
        sum(dnorm(sd$Y, eta, rep(pars$sigma, each = sd$n_obs), log = TRUE))
      }
  }
  c(likelihood = likelihood, prior = prior, total = likelihood + prior)
}

# Used only for the gradient oracle. Stan supplies its *sampled* natural
# coordinates and standard positivity transforms; reconstruct reference values
# independently, without requesting Stan's transformed parameters.
natural_stan_parameters <- function(pars, sd) {
  for (block in c("u", "z_u")) if (!is.null(pars[[block]])) {
    if (is.list(pars[[block]])) {
      values <- array(numeric(sd$p * sd$J * sd$n_re), c(sd$p, sd$J, sd$n_re))
      for (j in seq_len(sd$p)) values[j, , ] <- pars[[block]][[j]]
      pars[[block]] <- values
    } else {
      # utils::relist(), used by CmdStanR 0.9, drops rank-three dimensions.
      pars[[block]] <- array(pars[[block]], c(sd$p, sd$J, sd$n_re))
    }
  }
  for (block in c("kappa", "kappa_centered"))
    if (is.list(pars[[block]])) pars[[block]] <- do.call(rbind, pars[[block]])
  if (sd$n_fe == 0L) pars$beta <- matrix(numeric(), 0L, sd$p)
  pars
}

reference_from_sampled <- function(sampled, sd, family) {
  sampled <- natural_stan_parameters(sampled, sd)
  pars <- sampled
  if (family != "ordinal") pars$beta <- sampled$beta_centered
  if (family == "ordinal") pars$kappa <- sampled$kappa_centered
  for (node in seq_len(sd$p)) {
    shift <- sum(colMeans(sd$B) * sampled$phi[, node])
    if (family == "ordinal") {
      shift <- shift + sum(colMeans(sd$X) * sampled$beta[, node])
      pars$kappa[node, ] <- sampled$kappa_centered[node, ] + shift
    } else {
      shift <- shift + sum(colMeans(sd$X)[-1] * sampled$beta_centered[-1, node])
      pars$beta[1, node] <- sampled$beta_centered[1, node] - shift
    }
  }
  pars
}

centered_fixed_fit <- function(sd, pars, family) {
  init <- .bvarnet_centered_init(pars, sd, family, 1)[[1]]
  model <- compile_stan_source(.bvarnet_model_for_family(family))
  fit <- model$sample(data = sd, init = list(init), fixed_param = TRUE,
               chains = 1, parallel_chains = 1, iter_warmup = 0, iter_sampling = 1,
               seed = 1, refresh = 0, sig_figs = 18, show_messages = FALSE)
  # Some CmdStanR versions require this before density/gradient transformations.
  fit$init_model_methods()
  fit
}
