#how to handle interaction effects?
# TODO: Interaction effects
# TODO: set up to write it like formula
# TODO: make centering predictors option


# dont remove first row, but set lags to 0?

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
                         skip_lag = TRUE) { # skipping lag mechanism on/off
  ## Input
  # df: data in "long" format, cols: id, time, y, covariates

  ## Output
  # list with
  #(
  # p: no parameters
  # q: no covariates
  # w: no interactions
  # J: no people/gorups
  # T: no timepoints
  # K: AR(K)
  # n_obs: length Y, J*T-K
  # n_fe: no FEs
  # n_re: no RE
  # id: person identifier
  # Y: outcome matrix in long format
  # X: design matrix FE
  # B: design matrix lags
  # Z: design matrix RE
  #)

  ## Checks
  family <- match.arg(family, choices = c("bernoulli", "ordinal"))
  K <- as.integer(K) # ensure its an integer to not break fun
  stopifnot(is.logical(skip_lag), length(skip_lag) == 1L, !is.na(skip_lag))
  ##

   

  ## ensure the data is properly ordered by id and time to prevent data wrangling errors
  data <- data[order(data[[id_col]], data[[time_col]]), ]

  ## Listwise deletion
  na_action <- match.arg(na_action)

  if (na_action == "listwise") {
    check_cols <- c(y_cols, x_cols)
    complete <- complete.cases(data[, check_cols, drop = FALSE])
    n_na <- sum(!complete)
    data <- data[complete , , drop = FALSE]
  }

  ## "summary statistics"
  ids_unique <- unique(data[[id_col]])
  J <- length(ids_unique)
  p <- length(y_cols)
  q <- length(x_cols)
  PK <- p*K
  

  if (family == "ordinal") { # Check C ## Todo: add recoding??
  if (any(Y < 1L | Y > C, na.rm = TRUE)) {
    stop(sprintf("Ordinal Y values must be in 1:%d. Found values outside this range.", C))
  }
}

  df_split <- split(data, data[[id_col]])

  #n_obs <- J * (T_obs - K)
  # total number of modeled rows (one per observed row after first K within each subject)
  n_obs_initial <- n_obs <- sum(vapply(df_split, function(sub) max(0L, nrow(sub) - K), integer(1)))

  ## initialize design matrices and id
  Y <- matrix(NA_integer_, n_obs, p)
  X <- matrix(NA_real_, n_obs, q)
  B <- matrix(0, n_obs, PK)
  id_out <- integer(n_obs)


  ## begin creation of design matrices
  row <- 0L # bookkeeping row number outcome

  for(jj in seq_len(J)) { # safe with J = 0
    this_id <- ids_unique[jj]
    df_sub <- df_split[[as.character(this_id)]]
    df_sub <- df_sub[order(df_sub[[time_col]]),, drop = FALSE] # ensure ordering stays!, drop = F to ensure dimension stays the same

    Ti <- nrow(df_sub)
    if (Ti <= K) next # skips subjects that cant contribute to likelihood, atleast K+1 obs needed per subject

    times <- df_sub[[time_col]]
    Ymat <- as.matrix(df_sub[, y_cols, drop = FALSE])
    Xmat <- as.matrix(df_sub[, x_cols, drop = FALSE])

    ## begin onstruction of B
    for (t in (K+1L):Ti) {
      row <- row + 1L

      id_out[row] <- jj
      Y[row, ] <- Ymat[t, ]
      X[row, ] <- Xmat[t, ]

      # missing data handling: if we skip lags as mechanism this is what we use. TODO: make an option

      valid <- TRUE
      for(lag in 1:K) {
        if ((times[t] - times[t - lag]) != lag) {
          valid <- FALSE
          break
        }
      }

      if(valid) {
        for(lag in 1:K) {
          B[row, ((lag - 1L)*p+1L):(lag*p)] <- Ymat[t - lag, ]
        }
      } else if (!skip_lag) {
        row <- row - 1L # deletes created row
        next
      }
    }
  }

  if (row < n_obs) { # trim in case of listwise deletion without lag skipping
    Y <- Y[seq_len(row), , drop = FALSE]
    X <- X[seq_len(row), , drop = FALSE]
    B <- B[seq_len(row), , drop = FALSE]
    id_out <- id_out[seq_len(row)]
    n_obs <- row
  }

  if (n_obs == 0L) { # if all obs removed, stop function
    stop("All observations removed after missing-data handling. ",
         "Check your data for NAs and irregular time gaps.")
  }
  if(family == "ordinal") { 
    C <- max(Y, na.rm = TRUE)
    if (min(Y, na.rm = TRUE) < 1L) {
      stop("Ordinal Y must be coded as integers starting at 1. Found values < 1.")
    }
  }

  n_dropped <- n_obs_initial - n_obs

  if (n_dropped > 0) {
    message(sprintf("bvarnet: %d row(s) removed (na_action = '%s', skip_lag = %s). %d rows remain.",
                n_dropped, na_action, skip_lag, n_obs))
  }

  # center and add intercept ##make centering predictors an option
  if (center_x == TRUE) {
    means_X <- colMeans(X)
    X <- sweep(X, 2, means_X, "-")

  }

  if(family == "bernoulli") {
    X <- cbind(Intercept = 1, X) # add option to not add intercepts ( = constrain intercept to 0)
  }
  
  if (family == "bernoulli") { #naming of X
    colnames(X) <- c("Intercept", x_cols) # name X
  } else {
    colnames(X) <- x_cols
  }


  ## build FE interactions, subj to change
  tmp <- add_terms_to_X(X, B, fe_interactions)
  X <- tmp$X

  ##

  Z <- build_Z(X, B, re_cols = re_cols, re_temporal = re_temporal)

  ## build RE interactions, subj to change
  Z <- add_re_interactions_from_X(Z, X, B, re_interactions)

  ##
  ## name
  colnames(Y) <- y_cols # keep y names
  b_names <- unlist(lapply(1:K, function(lag) paste0("lag", lag, "_", y_cols)))
  colnames(B) <- b_names # name B
  
  
  out <- list(
  p = p, 
  q = q, 
  J = J,
  T = max(data[[time_col]], na.rm = TRUE),
  K = K, 
  n_obs = n_obs,
  n_fe = ncol(X), 
  n_re = ncol(Z),
  id = id_out, 
  Y = Y, 
  X = X, 
  B = B, 
  Z = Z,
  prior_beta_fam = 1, beta_loc = 0, beta_scale = 1, beta_df = 1000,
  prior_phi_fam = 1, phi_loc = 0, phi_scale = 1, phi_df = 1000,
  prior_sd_fam = 1, sd_loc = 0, sd_scale = 1, sd_df = 1000
)

  if (family == "ordinal") {
    out$C <- C
    out$prior_kappa_fam <- 1
    out$kappa_loc <- 0
    out$kappa_scale <- 1
    out$kappa_df <- 1000
  }

return(out)
  
}

build_Z <- function(X, B, re_cols = character(0), re_temporal = FALSE) {
  # X: n_obs x n_fe (incl intercept)
  # B: n_obs x (p*K)
  stopifnot(is.matrix(X), is.matrix(B))

  Z_list <- list()

  # random slopes on selected fixed-effect columns
  if (length(re_cols) > 0) {
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
