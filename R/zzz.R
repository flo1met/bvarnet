#' Check whether the precompiled Stan models are available
#'
#' Returns the names of the bundled Stan models whose compiled executables are
#' missing from the installed package. instantiate compiles the models at
#' install time, but only when CmdStan is available. On CRAN binaries
#' (built without CmdStan) the executables are absent and the models must be
#' compiled by installing from source on a machine with CmdStan set up.
#'
#' @return Character vector of model names lacking a compiled executable.
#' @noRd
.bvarnet_missing_models <- function() {
  models <- c("model_gaussian", "model_binary", "model_ordinal")
  ext <- if (.Platform$OS.type == "windows") ".exe" else ""
  missing <- character(0)
  for (name in models) {
    exe <- system.file(file.path("bin", "stan", paste0(name, ext)),
                       package = "bvarnet")
    if (!nzchar(exe) || !file.exists(exe)) {
      missing <- c(missing, name)
    }
  }
  missing
}

#' Detect first-time vs. post-upgrade "models not set up" states
#'
#' The model cache is keyed by `<installed version>-<source hash>` (see
#' `.bvarnet_cache_key()` in R/download_models.R), so a package upgrade
#' silently invalidates the cache: the new version looks in a new subdir. This
#' distinguishes "nothing has ever been set up" from "something was set up,
#' but for a different version" using only local filesystem checks -- never
#' the network, never erroring -- so it is safe to call from `.onAttach()`.
#'
#' @return One of `"ok"`, `"stale"`, or `"absent"`.
#' @noRd
.bvarnet_models_state <- function() {
  if (length(.bvarnet_missing_models()) == 0L) return("ok")
  cache_root <- .bvarnet_cache_root()
  cur_key <- tryCatch(.bvarnet_cache_key(), error = function(e) NA_character_)
  if (!is.na(cur_key) && .bvarnet_cache_complete(file.path(cache_root, cur_key))) {
    return("ok")
  }
  if (!dir.exists(cache_root)) return("absent")
  others <- setdiff(list.dirs(cache_root, recursive = FALSE, full.names = FALSE), cur_key)
  if (length(others)) "stale" else "absent"
}

#' Stop with an actionable, state-aware message when a model isn't available
#' @param model_name Character scalar naming the Stan model.
#' @return Does not return; always raises an error.
#' @noRd
.bvarnet_models_missing_error <- function(model_name) {
  state <- tryCatch(.bvarnet_models_state(), error = function(e) "absent")
  if (identical(state, "stale")) {
    stop(
      "bvarnet was updated to ", utils::packageVersion("bvarnet"), "; its Stan models ",
      "were set up for a previous version.\n",
      "Run bvarnet_setup_models() to set them up for this version, then retry.\n",
      "See ?bvarnet_setup_models for details.",
      call. = FALSE
    )
  }
  stop(
    "bvarnet: Stan models are not set up.\n",
    "Run bvarnet_setup_models() to either download precompiled binaries or ",
    "compile locally, then retry.\n",
    "See ?bvarnet_setup_models for details.",
    call. = FALSE
  )
}

.onAttach <- function(libname, pkgname) {
  tryCatch({
    state <- .bvarnet_models_state()
    if (identical(state, "ok")) return(invisible())
    msg <- if (identical(state, "stale")) {
      paste0(
        "bvarnet was updated to ", utils::packageVersion("bvarnet"), "; its Stan models ",
        "were set up for a previous version.\n",
        "Run bvarnet_setup_models() to set them up for this version."
      )
    } else {
      paste0(
        "bvarnet: the precompiled Stan models are not set up.\n",
        "Run bvarnet_setup_models() to either download precompiled binaries or ",
        "compile locally.\n",
        "Calls to bvar() will fail until this is done.\n",
        "See ?bvarnet_setup_models for details."
      )
    }
    packageStartupMessage(msg)
  }, error = function(e) invisible())
}
