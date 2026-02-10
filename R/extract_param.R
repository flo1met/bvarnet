extract_param <- function(object) {
  stopifnot(inherits(object, "bvarnet"))

  sd <- object$standata
  nm <- get_param_names(sd)

  # ---------- Intercepts & fixed effects (beta) ----------
  draws_beta <- extract_draws(object, "beta")
  beta_tab   <- build_summary_table(draws_beta, nm$fe, nm$y, "placeholder")
  beta_tab$type <- ifelse(beta_tab$predictor == nm$fe[1],
                          "Intercept", "Fixed Effect")

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

  out <- list(
    beta     = beta_tab,
    phi      = phi_tab,
    re_sd    = re_sd_tab,
    standata = sd,
    fit      = object$fit
  )
  class(out) <- "bvarnet_params"
  out
}
