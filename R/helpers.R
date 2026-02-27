## ---- resolve variable names from standata ----
get_param_names <- function(sd) {
  p    <- sd$p
  K    <- sd$K
  n_fe <- sd$n_fe
  n_re <- sd$n_re

  y_names  <- colnames(sd$Y)  %||% paste0("y", seq_len(p))
  fe_names <- colnames(sd$X)  %||% paste0("fe", seq_len(n_fe))
  b_names  <- colnames(sd$B)  %||% paste0("lag", rep(seq_len(K), each = p), "_", rep(y_names, K))
  re_names <- if (n_re > 0) (colnames(sd$Z) %||% paste0("re", seq_len(n_re))) else character(0)

  list(y = y_names, fe = fe_names, b = b_names, re = re_names)
}


## ---- build a labelled data.frame from a draws matrix ----
## draws: iterations x parameters matrix (from extract_draws)
## row_names / col_names map to the [row, col] Stan indices
## column order in draws must follow row-major: (1,1), (2,1), ..., (nr,1), (1,2), ...
build_summary_table <- function(draws, row_names, col_names, type) {
  nr   <- length(row_names)
  nc   <- length(col_names)
  ncol_draws <- ncol(draws)
  stopifnot(ncol_draws == nr * nc)

  d_mean   <- colMeans(draws)
  d_median <- apply(draws, 2, stats::median)
  d_q5     <- apply(draws, 2, stats::quantile, probs = 0.05)
  d_q95    <- apply(draws, 2, stats::quantile, probs = 0.95)

  data.frame(
    type      = rep(type, nr * nc),
    predictor = rep(row_names, times = nc),
    outcome   = rep(col_names, each  = nr),
    mean      = as.numeric(d_mean),
    median    = as.numeric(d_median),
    q5        = as.numeric(d_q5),
    q95       = as.numeric(d_q95),
    stringsAsFactors = FALSE
  )
}


## ---- extract posterior draws with readable column names ----
extract_draws <- function(object, parameter = c("beta", "phi", "sd_u", "sigma", "kappa")) {
  stopifnot(inherits(object, "bvarnet") || inherits(object, "bvarnet_params"))
  parameter <- match.arg(parameter, c("beta", "phi", "sd_u", "sigma", "kappa"))

  if (parameter == "sigma" && object$family != "gaussian") {
    stop("Parameter 'sigma' only exists for gaussian models.")
  }
  if (parameter == "kappa" && object$family != "ordinal") {
    stop("Parameter 'kappa' only exists for ordinal models.")
  }

  fit <- object$fit
  sd  <- object$standata

  nm <- get_param_names(sd)

  if (parameter == "sd_u" && sd$n_re == 0) {
    stop("Parameter 'sd_u' not available — model has no random effects (n_re = 0).")
  }

  

  p    <- sd$p
  K    <- sd$K
  n_fe <- sd$n_fe
  n_re <- sd$n_re

  draws     <- fit$draws(variables = parameter, format = "matrix")
  old_names <- colnames(draws)
  new_names <- old_names

  if (parameter == "beta") {
    for (j in seq_len(p)) {
      for (i in seq_len(n_fe)) {
        old <- paste0("beta[", i, ",", j, "]")
        new_names[old_names == old] <- paste0(nm$fe[i], " -> ", nm$y[j])
      }
    }
  } else if (parameter == "phi") {
    for (j in seq_len(p)) {
      for (i in seq_len(p * K)) {
        old <- paste0("phi[", i, ",", j, "]")
        new_names[old_names == old] <- paste0(nm$b[i], " -> ", nm$y[j])
      }
    }
  } else if (parameter == "sd_u") {
    for (j in seq_len(n_re)) {
      for (i in seq_len(p)) {
        old <- paste0("sd_u[", i, ",", j, "]")
        new_names[old_names == old] <- paste0("sd(", nm$re[j], " | ", nm$y[i], ")")
      }
    }
  } else if (parameter == "sigma") {
    for (j in seq_len(p)) {
      old <- paste0("sigma[", j, "]")
      new_names[old_names == old] <- paste0("sigma(", nm$y[j], ")")
    }
  } else if (parameter == "kappa") {
    for (j in seq_len(p)) {
      for (c in seq_len(sd$C - 1)) {
        old <- paste0("kappa[", j, ",", c, "]")
        new_names[old_names == old] <- paste0("kappa(", nm$y[j], ", c", c, ")")
      }
    }
  } 

  colnames(draws) <- new_names
  draws
}


## ---- print method for bvarnet_params ----
print.bvarnet_params <- function(x, ...) {
  rule <- function(title) {
    cat("\n", title, "\n", strrep("-", nchar(title) + 4), "\n", sep = "")
  }

  fmt <- function(df) {
    num_cols <- c("mean", "median", "q5", "q95")
    for (col in num_cols) df[[col]] <- sprintf("% .3f", df[[col]])
    df
  }

  cat("BVAR Network \u2014 Parameter Summary\n")
  cat(strrep("=", 40), "\n")

  # Intercepts
  intercepts <- x$beta[x$beta$type == "Intercept", ]
  if (nrow(intercepts) > 0) {
    rule("Intercepts")
    print(fmt(intercepts), row.names = FALSE)
  }

  # Fixed effects (non-intercept)
  fes <- x$beta[x$beta$type == "Fixed Effect", ]
  if (nrow(fes) > 0) {
    rule("Fixed Effects")
    print(fmt(fes), row.names = FALSE)
  }

  # Temporal
  rule("Temporal Effects")
  print(fmt(x$phi), row.names = FALSE)

  # Random-effect SDs
  if (!is.null(x$re_sd) && nrow(x$re_sd) > 0) {
    rule("Random Effect SDs")
    print(fmt(x$re_sd), row.names = FALSE)
  }

  if (!is.null(x$kappa) && nrow(x$kappa) > 0) {
    rule("Thresholds")
    print(fmt(x$kappa), row.names = FALSE)
  }

  if (!is.null(x$sigma) && nrow(x$sigma) > 0) {
    rule("Residual SD")
    print(fmt(x$sigma), row.names = FALSE)
  }

  cat("\n")
  invisible(x)
}
