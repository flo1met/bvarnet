#' Fit a Bayesian multilevel VAR network model
#'
#' The `bvar` function estimates the posterior distribution of the specified Bayesian (Multilevel) Vector Autoregression.
#'
#' @param id_col Character. Name of the subject/group identifier column.
#' @param time_col Character scalar naming the time-index column in `data`.
#'   Its meaning depends on `beep_col`:
#'   \itemize{
#'     \item If `beep_col = NULL` (default), this column contains the scheduled
#'       occasion index for one sequence per subject. Values must be finite and
#'       integer-valued, with one unit representing one scheduled lag step.
#'       Indices must continue across the subject's sequence rather than restart
#'       each day. For lag k, the predecessor must have index `time - k`.
#'       No separate day or overnight boundary is inferred in this mode.
#'     \item If `beep_col` names a column, `time_col` identifies the day or session.
#'       `beep_col` then contains finite integer-valued beep indices within each day/session, 
#'        with one unit representing one scheduled lag step. Beep indices may restart in each day/session.
#'       For lag k, the predecessor must have index `beep - k` in the same subject and day/session.
#'   }
#'   Preserve the original scheduled labels when observations are missing.
#'   Rows without a complete, correctly spaced lag window are excluded as
#'   responses from the likelihood; their observed outcomes may still be used
#'   as predictors for later eligible responses.
#' @param beep_col NULL (default), or a character scalar naming the scheduled
#'   occasion index within each day/session, specified in the `time_col` column. 
#'   If NULL, `time_col` itself supplies the occasion index for one sequence per subject, 
#'   and no day/session splitting is performed. If supplied, `time_col` must identify
#'   the day/session and this column must contain finite integer-valued beep
#'   indices, with one unit representing one scheduled lag step. Beep indices
#'   may restart in each day/session. For lag k, the predecessor must have
#'   index `beep - k` in the same subject and day/session. 
#' @param y_cols Character vector. Names of the outcome columns.
#' @param x_cols Character vector or NULL. Names of the covariate columns.
#' @param center_x Logical. Grand-mean centre covariates before fitting?
#'   Default \code{FALSE}.
#' @param fe_interactions List or NULL. Fixed-effect interaction terms to add
#'   to the design matrix. Each element is a character vector of column names
#'   to interact, or \code{c("lag", "x")} to interact all lag columns with
#'   a covariate.
#' @param re_interactions List or NULL. Random-effect interaction terms.
#' @param re_cols Character vector. Columns from X and/or "Intercept" to include as random slopes.
#' @param re_temporal Logical. Include random slopes on lag predictors?
#'   Default \code{FALSE}.
#' @param K Integer. AR order. Default 1.
#' @param na_action Character. Missing-data strategy; currently only
#'   \code{"listwise"}.
#' @param data Data frame in long format.
#' @param family Character scalar or vector. Observation model per node.
#'   A scalar is recycled to all \code{y_cols}. A vector of length
#'   \code{length(y_cols)} (named or positional) specifies per-node families.
#'   Valid values: \code{"bernoulli"}, \code{"ordinal"}, \code{"gaussian"}.
#' @param priors A \code{bvarnet_priors} object from \code{set_priors()}.
#'   Defaults to \code{set_priors()} (package defaults).
#' @param iter Integer. Number of post-warmup iterations per chain. Default 4000.
#' @param warmup Integer. Number of warmup iterations per chain. Default 1000.
#' @param chains Integer. Number of MCMC chains. Default 4.
#' @param cores Integer. Number of chains to run in parallel. Default 1.
#' @param seed Integer or NULL. RNG seed.
#' @param adapt_delta Numeric in (0, 1). Target average proposal acceptance
#'   probability during warmup adaptation. Higher values (e.g., 0.95–0.99)
#'   reduce divergences at the cost of slower sampling. Default \code{NULL}
#'   (CmdStan default of 0.8).
#' @param max_treedepth Integer. Maximum depth of the NUTS binary tree.
#'   Increasing this allows the sampler to take more leapfrog steps per
#'   iteration, which can help with difficult posteriors (e.g., funnels in
#'   hierarchical logistic models) but increases computation. Default
#'   \code{NULL} (CmdStan default of 10).
#' @param save_data Logical. Keep the Stan data used for estimation in the
#'   returned object's \code{standata} slot. Default \code{TRUE}. Set to
#'   \code{FALSE} to remove these data from the fitted model before saving or
#'   sharing it.
#' @param ... Additional named arguments forwarded to the CmdStanR
#'   \code{$sample()} method, e.g. \code{init}, \code{refresh}, \code{thin},
#'   \code{step_size}, or \code{show_messages}. See
#'   \code{?cmdstanr::"model-method-sample"} for the full list. Arguments that
#'   \code{bvar()} sets itself (\code{data}, \code{seed}, \code{iter_warmup},
#'   \code{iter_sampling}, \code{chains}, \code{parallel_chains},
#'   \code{adapt_delta}, \code{max_treedepth}, and CmdStanR's deprecated
#'   aliases for them) are rejected with an error; use the corresponding
#'   \code{bvar()} argument instead.
#'
#' @section Scheduled time and missing observations:
#' Index columns identify eligible transitions; they do not automatically add
#' time or day effects as covariates. A daily series with one observation per
#' day can use `time_col = "day"` alone: consecutive days are ordinary lag steps.
#' With multiple observations per day, supply both day and beep columns to
#' exclude overnight transitions. Day/session labels may be numeric, character,
#' factor or Date; occasion indices must be numeric integers. Retained keys must
#' be unique within subject (and day/session when supplied), with no missing keys.
#'
#' Only complete lag windows enter the likelihood. For indices `1, 2, 4, 5`,
#' a VAR(1) uses `1 -> 2` and `4 -> 5`, excluding `2 -> 4`. With day/beep indices,
#' day 1 beep 5 never predicts day 2 beep 1, while day 2 beep 1 can predict day 2
#' beep 2. Subject random effects remain shared across that subject's days.
#' Occasion-level listwise deletion is applied to all selected outcomes and
#' covariates before constructing lags. This is not a general correction for
#' missingness under MAR or MNAR.
#'
#' The former `skip_lag` argument and zero-filled gap likelihood were removed
#' because they could bias temporal and non-temporal parameters. Passing the
#' removed argument is an error. Fits may use fewer responses and give different
#' estimates than the former default.
#'
#' @section Internal predictor centering:
#' Stan samples intercepts (Gaussian/Bernoulli) or thresholds (ordinal) using
#' centered completed fixed-effect and lag columns. Reported \code{beta},
#' \code{kappa}, priors and Bayes-factor hypotheses retain the reference of
#' the public design produced with the requested \code{center_x} setting.
#' Random-effect designs and priors are unchanged. This coordinate change can
#' improve sampling, but does not guarantee a speed improvement.
#'
#' @return A \code{bvarnet} object (a named list) with slots:
#'   \code{draws}, \code{convergence}, \code{diagnostics}, \code{timing},
#'   \code{metadata}, \code{return_codes}, \code{family}, \code{model_info},
#'   \code{priors}, and \code{priors_effective}. The \code{standata} slot
#'   contains the estimation data when \code{save_data = TRUE} (default) and
#'   is \code{NULL} when \code{save_data = FALSE}. The \code{model_info} slot
#'   retains only dimensions and variable/interaction labels needed for
#'   summaries and parameter extraction. No additional data frame is stored.
#'
#' @section Requested versus effective priors:
#' \code{priors} holds the priors as you specified them. \code{priors_effective}
#' holds the priors Stan actually sampled under: one \code{bvarnet_priors}
#' object per outcome. The two differ for Gaussian outcomes left at their
#' default priors, where \code{intercept}, \code{beta} and \code{sigma} scales
#' are multiplied by the outcome SD so that unit-scale defaults stay weakly
#' informative on the raw data scale. User-supplied priors are never rescaled.
#'
#' Bayes factors (\code{\link{bf_table}}) divide by \code{priors_effective},
#' as the Savage-Dickey density ratio requires. Note that the joint path
#' (\code{family} a single value) scales every outcome by the mean SD across
#' outcomes, whereas the nodewise path (mixed \code{family}) scales each outcome
#' by its own SD, so the same data can imply slightly different effective priors
#' depending on which path runs.
#'
#' @seealso \code{\link{bvarnet_setup_models}}, which must be run once before
#'   the first \code{bvar()} call to set up the required Stan models (either
#'   by downloading precompiled binaries or compiling them locally).
#'
#' @examples
#' # Scheduled indices retain missing occasions rather than compressing time.
#' single_index <- data.frame(id = 1, occasion = c(1, 2, 4, 5), y = c(1, 2, 3, 4))
#' # time_col = "occasion": transitions 1 -> 2 and 4 -> 5 are eligible.
#' day_beep <- data.frame(id = 1, day = c(1, 1, 1, 2, 2),
#'                        beep = c(1, 2, 4, 1, 2), y = c(1, 2, 3, 4, 5))
#' # time_col = "day", beep_col = "beep": only within-day beeps 1 -> 2.
#' \dontrun{
#' # Run bvar on studentlife data
#' data(studentlife, package = "bvarnet")
#' fit <- bvar(
#'   id_col = "id",
#'   time_col = "day",
#'   y_cols = c("anxious", "calm", "conventional", "critical", "dependable"),
#'   re_temporal = TRUE,
#'   K = 1,
#'   data = studentlife,
#'   family = "ordinal",
#'   priors = set_priors(),
#'   seed = 1337)
#'
#' summary(fit)
#' }

#' @export
bvar <- function(id_col,
                 time_col,
                 beep_col = NULL,
                 y_cols,
                 x_cols = NULL,
                 center_x = FALSE,
                 fe_interactions = NULL,
                 re_interactions = NULL,
                 re_cols = NULL,
                 re_temporal = FALSE,
                 K = 1,
                 na_action = c("listwise"),
                 data,
                 family = c("bernoulli", "ordinal", "gaussian"),
                 priors = set_priors(),
                 iter = 4000,
                 warmup = 1000,
                 chains = 4,
                 cores = 1,
                 seed = NULL,
                 adapt_delta = NULL,
                 max_treedepth = NULL,
                 save_data = TRUE,
                 ...

  ) {

  if (!is.logical(save_data) || length(save_data) != 1L || is.na(save_data))
    stop("save_data must be TRUE or FALSE.", call. = FALSE)

  dots <- .check_sampler_dots(...)

  family_vec <- .parse_family(family, y_cols)
  is_mixed   <- length(unique(family_vec)) > 1L

  # Validate and ensure all required prior slots exist
  if (!inherits(priors, "bvarnet_priors"))
    stop("'priors' must be a bvarnet_priors object from set_priors().",
         call. = FALSE)
  priors <- .ensure_prior_slots(priors, family_vec)

  if (is_mixed) {
    return(.bvar_nodewise(
      id_col = id_col, time_col = time_col, y_cols = y_cols,
      x_cols = x_cols, center_x = center_x,
      fe_interactions = fe_interactions,
      re_interactions = re_interactions,
      re_cols = re_cols, re_temporal = re_temporal,
      K = K, na_action = na_action, beep_col = beep_col,
      data = data, family_vec = family_vec, priors = priors,
      iter = iter, warmup = warmup, chains = chains, cores = cores,
      seed = seed, adapt_delta = adapt_delta, max_treedepth = max_treedepth,
      save_data = save_data, dots = dots
    ))
  }

  # --- existing joint path (homogeneous family) ---
  family <- unname(family_vec[1])  # scalar for backward compat
  model_name <- .bvarnet_model_for_family(family)

  standata <- to_stan_data(data = data,
                            family = family,
                            id_col = id_col,
                            time_col = time_col,
                            y_cols = y_cols,
                            x_cols = x_cols,
                            center_x = center_x,
                            fe_interactions = fe_interactions,
                            re_interactions = re_interactions,
                            re_cols = re_cols,
                            re_temporal = re_temporal,
                            K = K,
                            na_action = na_action,
                            beep_col = beep_col,
                            priors = priors
                           )

  # Prior warnings (uses actual n_re from built design)
  priors_needed <- .prior_warnings(priors, family_vec, standata$n_re)
  resolved <- .bvarnet_stan_model(model_name)
  stanmodel <- resolved$model

  if ("init" %in% names(dots))
    dots$init <- .bvarnet_centered_init(dots$init, standata, family, chains)

  sample_call <- function() {
    do.call(stanmodel$sample, c(
      list(data = standata[!names(standata) %in%
                                  c("fe_interaction_terms",
                                    "fe_interaction_colnames",
                                    "id_levels",
                                    "x_center_means",
                                    "time_obs",
                                    "index_obs",
                                    "row_map",
                                    "n_rows_data",
                                    "design_spec")],
           seed = seed,
           iter_warmup = warmup,
           iter_sampling = iter,
           chains = chains,
           parallel_chains = cores,
           adapt_delta = adapt_delta,
           max_treedepth = max_treedepth),
      dots
    ))
  }
  # Cache-resolved models need the bundled TBB findable at $sample() time
  # install-tree models keep cmdstanr's own resolution untouched.
  stanfit <- if (!is.null(resolved$cache_dir)) {
    .bvarnet_with_libpath(file.path(resolved$cache_dir, "lib", "tbb"), sample_call())
  } else {
    sample_call()
  }

  # Extract everything from CmdStanMCMC into plain base-R objects, then discard
  # the fit object (CSV refs, compiled binary, lazy draws) to keep memory lean.
  raw_draws   <- .bvarnet_public_draws(stanfit$draws(format = "array"))
  draws       <- unclass(raw_draws)         # strip draws_array class; dimnames preserved
  attr(draws, "class") <- NULL              # ensure it is a plain array

  # Compute convergence diagnostics from posterior draws
  conv_tbl    <- posterior::summarise_draws(raw_draws,
                   posterior::rhat, posterior::ess_bulk, posterior::ess_tail)
  convergence <- as.data.frame(conv_tbl)
  names(convergence) <- gsub("^posterior::", "", names(convergence))

  diagnostics  <- as.data.frame(stanfit$diagnostic_summary())
  timing       <- stanfit$time()
  metadata     <- stanfit$metadata()
  return_codes <- stanfit$return_codes()

  out <- structure(
    list(
      draws         = draws,
      convergence   = convergence,
      diagnostics   = diagnostics,
      timing        = timing,
      metadata      = metadata,
      return_codes  = return_codes,
      family        = family_vec,
      standata      = if (save_data) standata else NULL,
      model_info    = .bvarnet_model_info(standata),
      priors        = priors,
      # One s_y scales every outcome in the joint path, so all nodes share the
      # same effective prior. Stored per node anyway, to match the mixed path.
      priors_effective = stats::setNames(
        rep(list(attr(standata, "priors_effective")), length(family_vec)),
        names(family_vec)
      ),
      priors_needed = priors_needed
    ),
    class = "bvarnet"
  )
  out
}


# ---- .parse_family() — validate & normalise family to named character vector ----

#' Parse and validate the family argument
#'
#' Returns a named character vector of length \code{p} with names from
#' \code{y_cols}.
#'
#' @param family Character scalar or vector.
#' @param y_cols Character vector of outcome names.
#' @return Named character vector.
#' @keywords internal
#' @noRd
.parse_family <- function(family, y_cols) {
  valid <- c("bernoulli", "ordinal", "gaussian")
  p <- length(y_cols)

  if (length(family) == 1L && is.null(names(family))) {
    family <- match.arg(family, valid)
    out <- stats::setNames(rep(family, p), y_cols)
  } else if (length(family) == p && is.null(names(family))) {
    bad <- setdiff(family, valid)
    if (length(bad) > 0L)
      stop("Invalid family value(s): ", paste(bad, collapse = ", "),
           ". Must be one of: ", paste(valid, collapse = ", "), call. = FALSE)
    out <- stats::setNames(family, y_cols)
  } else if (!is.null(names(family))) {
    if (!setequal(names(family), y_cols))
      stop("Names of 'family' must match y_cols: ",
           paste(y_cols, collapse = ", "), call. = FALSE)
    bad <- setdiff(family, valid)
    if (length(bad) > 0L)
      stop("Invalid family value(s): ", paste(bad, collapse = ", "),
           ". Must be one of: ", paste(valid, collapse = ", "), call. = FALSE)
    out <- family[y_cols]  # reorder to match y_cols
  } else {
    stop("'family' must be a scalar, a length-p vector, or a named vector ",
         "with names matching y_cols.", call. = FALSE)
  }
  out
}


# ---- .bvar_nodewise() — fit each node independently ----

#' @keywords internal
.bvar_nodewise <- function(id_col, time_col, y_cols, x_cols,
                           center_x, fe_interactions, re_interactions,
                           re_cols, re_temporal, K, na_action, beep_col,
                           data, family_vec, priors,
                           iter, warmup, chains, cores,
                           seed, adapt_delta, max_treedepth,
                           save_data = TRUE, dots = list()) {
  p <- length(y_cols)

  # --- Shared matrices (D3) ---
  shared <- .to_stan_data_shared(
    data = data, id_col = id_col, time_col = time_col,
    y_cols = y_cols, x_cols = x_cols, center_x = center_x,
    fe_interactions = fe_interactions, re_interactions = re_interactions,
    re_cols = re_cols, re_temporal = re_temporal, K = K,
    na_action = na_action, beep_col = beep_col
  )

  # Prior warnings (uses actual n_re from built design)
  priors_needed <- .prior_warnings(priors, family_vec, shared$n_re)

  # --- Resolve each node's model up front (Section C1) ---
  resolved <- lapply(family_vec, function(fam_j) {
    .bvarnet_stan_model(.bvarnet_model_for_family(fam_j))
  })
  # All cache-resolved models share one cache dir (keyed by package version +
  # source hash), so at most one distinct lib/tbb needs wrapping around the
  # whole loop rather than per iteration (many $sample() calls happen here).
  cache_dir <- Find(Negate(is.null), lapply(resolved, `[[`, "cache_dir"))

  # Built up front rather than inside run_sampling(): each node scales the
  # default Gaussian priors by its own sd(Y_node), and those effective priors
  # have to reach the fitted object for savage_dickey() to use.
  sd_nodes <- lapply(seq_len(p), function(j)
    .to_stan_data_node(shared, j, family_vec[j], priors))

  # Resolve every initializer before starting any node, so a malformed later
  # node cannot leave an expensive partially fitted model behind.
  dots_nodes <- lapply(seq_len(p), function(j) {
    node_dots <- .bvarnet_node_dots(dots, j)
    if ("init" %in% names(node_dots))
      node_dots$init <- .bvarnet_centered_init(node_dots$init, sd_nodes[[j]],
                                             family_vec[j], chains, node_id = j)
    node_dots
  })

  run_sampling <- function() {
    out <- vector("list", p)
    for (j in seq_len(p)) {
      node_dots <- dots_nodes[[j]]
      out[[j]] <- do.call(resolved[[j]]$model$sample, c(
        list(
          data            = sd_nodes[[j]],
          seed            = seed,
          iter_warmup     = warmup,
          iter_sampling   = iter,
          chains          = chains,
          parallel_chains = min(cores, chains),
          adapt_delta     = adapt_delta,
          max_treedepth   = max_treedepth
        ),
        node_dots
      ))
    }
    out
  }
  fits <- if (!is.null(cache_dir)) {
    .bvarnet_with_libpath(file.path(cache_dir, "lib", "tbb"), run_sampling())
  } else {
    run_sampling()
  }

  # --- Combine into bvarnet object ---
  priors_effective <- stats::setNames(
    lapply(sd_nodes, attr, "priors_effective"), names(family_vec)
  )
  .combine_nodewise_fits(fits, family_vec, shared, priors, priors_effective,
                         priors_needed, iter, chains, save_data = save_data)
}


# ---- .combine_nodewise_fits() — merge per-node fits into one bvarnet object ----

#' @keywords internal
.combine_nodewise_fits <- function(fits, family_vec, shared, priors,
                                   priors_effective, priors_needed,
                                   iter, chains, save_data = TRUE) {
  p <- length(family_vec)
  n_fe_full <- shared$n_fe   # includes Intercept

  # --- a) Combined draws array ---
  draw_chunks <- vector("list", p)
  lp_chunks   <- vector("list", p)
  conv_chunks <- vector("list", p)

  for (j in seq_len(p)) {
    raw_draws <- .bvarnet_public_draws(fits[[j]]$draws(format = "array"))
    arr <- unclass(raw_draws)
    attr(arr, "class") <- NULL
    par_names <- dimnames(arr)[[3]]

    # Extract and remove lp__
    lp_idx <- which(par_names == "lp__")
    if (length(lp_idx) > 0L) {
      lp_chunks[[j]] <- arr[, , lp_idx, drop = TRUE]  # [iter, chains]
      arr <- arr[, , -lp_idx, drop = FALSE]
      par_names <- par_names[-lp_idx]
    }

    # Rename parameters from p=1 to full-model indexing
    fam_j <- family_vec[j]
    new_names <- character(length(par_names))
    for (k in seq_along(par_names)) {
      new_names[k] <- .rename_node_param(par_names[k], j, fam_j, n_fe_full)
    }
    dimnames(arr)[[3]] <- new_names

    # For ordinal: insert NA sentinel for beta[1, j]
    if (fam_j == "ordinal") {
      sentinel_name <- sprintf("beta[1,%d]", j)
      sentinel_arr <- array(NA_real_,
                            dim = c(dim(arr)[1], dim(arr)[2], 1L),
                            dimnames = list(NULL, NULL, sentinel_name))
      arr <- abind_simple(sentinel_arr, arr)
    }

    draw_chunks[[j]] <- arr

    # Convergence: rename variables
    conv_tbl <- posterior::summarise_draws(raw_draws,
                  posterior::rhat, posterior::ess_bulk, posterior::ess_tail)
    conv_df <- as.data.frame(conv_tbl)
    names(conv_df) <- gsub("^posterior::", "", names(conv_df))
    conv_df <- conv_df[conv_df$variable != "lp__", , drop = FALSE]
    conv_df$variable <- vapply(conv_df$variable, function(nm) {
      .rename_node_param(nm, j, fam_j, n_fe_full)
    }, character(1L), USE.NAMES = FALSE)
    if (fam_j == "ordinal") {
      sentinel_row <- data.frame(
        variable = sprintf("beta[1,%d]", j),
        rhat = NA_real_, ess_bulk = NA_real_, ess_tail = NA_real_,
        stringsAsFactors = FALSE
      )
      conv_df <- rbind(sentinel_row, conv_df)
    }
    conv_chunks[[j]] <- conv_df
  }

  # Bind draws across 3rd dimension
  combined_draws <- do.call(function(...) abind_simple(...), draw_chunks)

  # Add summed lp__
  if (length(lp_chunks) > 0L && !is.null(lp_chunks[[1]])) {
    lp_sum <- Reduce(`+`, lp_chunks)
    lp_arr <- array(lp_sum, dim = c(dim(lp_sum), 1L),
                    dimnames = list(NULL, NULL, "lp__"))
    combined_draws <- abind_simple(combined_draws, lp_arr)
  }

  # --- b) Convergence ---
  combined_convergence <- do.call(rbind, conv_chunks)
  rownames(combined_convergence) <- NULL

  # --- c) Diagnostics (D8) ---
  diag_list <- lapply(fits, function(f) as.data.frame(f$diagnostic_summary()))
  diagnostics <- data.frame(
    num_divergent     = Reduce(`+`, lapply(diag_list, `[[`, "num_divergent")),
    num_max_treedepth = Reduce(`+`, lapply(diag_list, `[[`, "num_max_treedepth")),
    ebfmi             = Reduce(pmin, lapply(diag_list, `[[`, "ebfmi"))
  )

  # --- d) Timing, metadata, return codes ---
  timing   <- fits[[1]]$time()
  metadata <- fits[[1]]$metadata()
  rc_list  <- lapply(fits, function(f) f$return_codes())
  return_codes <- Reduce(pmax, lapply(rc_list, abs))

  # --- e) Augment shared standata (D6) ---
  standata_full <- shared
  stopifnot(standata_full$p == p, standata_full$K == shared$K)

  ord_idx <- which(family_vec == "ordinal")
  if (length(ord_idx) > 0L) {
    standata_full$C_per_node <- vapply(ord_idx, function(j)
      as.integer(max(shared$Y[, j])), integer(1L))
    names(standata_full$C_per_node) <- as.character(ord_idx)
  }

  # --- f) Return object ---
  structure(
    list(
      draws            = combined_draws,
      convergence      = combined_convergence,
      diagnostics      = diagnostics,
      node_diagnostics = diag_list,
      timing           = timing,
      metadata         = metadata,
      return_codes     = return_codes,
      family           = family_vec,
      standata         = if (save_data) standata_full else NULL,
      model_info       = .bvarnet_model_info(standata_full),
      priors           = priors,
      # Per node: each Gaussian node scaled its defaults by its own sd(Y_node).
      priors_effective = priors_effective,
      priors_needed    = priors_needed
    ),
    class = "bvarnet"
  )
}


# ---- helpers for .combine_nodewise_fits() ----

#' Rename a per-node (p=1) Stan parameter name to full-model naming
#' @keywords internal
#' @noRd
.rename_node_param <- function(name, j, family, n_fe_full) {
  if (grepl("^beta\\[", name)) {
    parts <- regmatches(name, regexec("^beta\\[(\\d+),(\\d+)\\]$", name))[[1]]
    k <- as.integer(parts[2])
    if (family == "ordinal") {
      return(sprintf("beta[%d,%d]", k + 1L, j))
    }
    return(sprintf("beta[%d,%d]", k, j))
  }
  if (grepl("^phi\\[", name)) {
    parts <- regmatches(name, regexec("^phi\\[(\\d+),(\\d+)\\]$", name))[[1]]
    m <- as.integer(parts[2])
    return(sprintf("phi[%d,%d]", m, j))
  }
  if (grepl("^sigma\\[", name)) {
    return(sprintf("sigma[%d]", j))
  }
  if (grepl("^kappa\\[", name)) {
    parts <- regmatches(name, regexec("^kappa\\[(\\d+),(\\d+)\\]$", name))[[1]]
    c_val <- as.integer(parts[3])
    return(sprintf("kappa[%d,%d]", j, c_val))
  }
  if (grepl("^sd_u\\[", name)) {
    parts <- regmatches(name, regexec("^sd_u\\[(\\d+),(\\d+)\\]$", name))[[1]]
    r <- as.integer(parts[3])
    return(sprintf("sd_u[%d,%d]", j, r))
  }
  if (grepl("^u\\[", name)) {
    parts <- regmatches(name, regexec("^u\\[(\\d+),(\\d+),(\\d+)\\]$", name))[[1]]
    i_val <- as.integer(parts[3])
    r_val <- as.integer(parts[4])
    return(sprintf("u[%d,%d,%d]", j, i_val, r_val))
  }
  if (grepl("^z_u\\[", name)) {
    parts <- regmatches(name, regexec("^z_u\\[(\\d+),(\\d+),(\\d+)\\]$", name))[[1]]
    i_val <- as.integer(parts[3])
    r_val <- as.integer(parts[4])
    return(sprintf("z_u[%d,%d,%d]", j, i_val, r_val))
  }
  name
}


#' Bind arrays along the third dimension (simple abind replacement)
#' @keywords internal
#' @noRd
abind_simple <- function(...) {
  arrays <- list(...)
  if (length(arrays) == 1L && is.list(arrays[[1]]) && !is.array(arrays[[1]]))
    arrays <- arrays[[1]]
  arrays <- Filter(Negate(is.null), arrays)
  if (length(arrays) == 0L) return(NULL)
  if (length(arrays) == 1L) return(arrays[[1]])

  d1 <- dim(arrays[[1]])[1]
  d2 <- dim(arrays[[1]])[2]
  all_names <- unlist(lapply(arrays, function(a) dimnames(a)[[3]]))
  total_d3 <- sum(vapply(arrays, function(a) dim(a)[3], integer(1L)))

  out <- array(NA_real_, dim = c(d1, d2, total_d3))
  offset <- 0L
  for (a in arrays) {
    d3 <- dim(a)[3]
    out[, , (offset + 1L):(offset + d3)] <- a
    offset <- offset + d3
  }
  dimnames(out) <- list(dimnames(arrays[[1]])[[1]],
                        dimnames(arrays[[1]])[[2]],
                        all_names)
  out
}
