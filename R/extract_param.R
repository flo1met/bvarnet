extract_param <- function(object) {
  stopifnot(inherits(object, "bvarnet"))



  sd <- object$standata
  nm <- get_param_names(sd)

  # ---------- Intercepts & fixed effects (beta) ----------
  draws_beta <- extract_draws(object, "beta")
  beta_tab   <- build_summary_table(draws_beta, nm$fe, nm$y, "placeholder")
  beta_tab$type <- ifelse(
    beta_tab$predictor == "Intercept",
    "Intercept", "Fixed Effect"
  )

  # ---------- Temporal effects (phi) ----------
  draws_phi <- extract_draws(object, "phi")
  phi_tab   <- build_summary_table(draws_phi, nm$b, nm$y, "Temporal")

  # ---------- Random-effect SDs (sd_u) ----------
  re_sd_tab <- if (sd$n_re > 0) {
    draws_sd <- extract_draws(object, "sd_u")
    tab <- build_summary_table(draws_sd, nm$y, nm$re, "Random Effect SD")
    colnames(tab) <- c("type", "outcome", "random_effect",
                       "mean", "median", "q5", "q95")
    tab
  } else NULL

  # ---------- Residual SD (sigma, gaussian only) ----------
  sigma_tab <- if (object$family == "gaussian") {
    draws_sigma <- extract_draws(object, "sigma")
    build_summary_table(draws_sigma, nm$y, "sigma", "Residual SD")
  } else NULL

  # ---------- Thresholds (kappa, ordinal only) ----------
  kappa_tab <- if (object$family == "ordinal") {
    draws_kappa <- extract_draws(object, "kappa")
    build_summary_table(draws_kappa, paste0("kappa[", seq_len(sd$C - 1), "]"), "kappa", "Threshold")
  } else NULL

  out <- list(
    beta     = beta_tab,
    phi      = phi_tab,
    re_sd    = re_sd_tab,
    standata = sd,
    fit      = object$fit
  )

  if (!is.null(sigma_tab)) out$sigma <- sigma_tab
  if (!is.null(kappa_tab)) out$kappa <- kappa_tab

  class(out) <- "bvarnet_params"
  out
}
