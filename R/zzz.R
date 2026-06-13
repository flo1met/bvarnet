#' Check whether the precompiled Stan models are available
#'
#' Returns the names of the bundled Stan models whose compiled executables are
#' missing from the installed package. instantiate compiles the models at
#' install time, but only when CmdStan is available. On CRAN/r-universe binaries
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

#' Stop with an installation hint if a model's executable is missing
#' @param model_name Character scalar naming the Stan model.
#' @return Invisibly `TRUE` if the model is compiled; otherwise an error.
#' @noRd
.check_compiled_model <- function(model_name) {
  ext <- if (.Platform$OS.type == "windows") ".exe" else ""
  exe <- system.file(file.path("bin", "stan", paste0(model_name, ext)),
                     package = "bvarnet")
  if (nzchar(exe) && file.exists(exe)) {
    return(invisible(TRUE))
  }
  stop(
    "bvarnet: the precompiled Stan model \"", model_name, "\" is not available, ",
    "so bvar() cannot run.\n",
    "This happens when the package is installed from a CRAN/r-universe binary, ",
    "or from source without CmdStan set up, so the models were never compiled.\n",
    "See the installation instructions at ",
    "https://github.com/flo1met/bvarnet#installation ",
    "for how to install cmdstanr + CmdStan and reinstall bvarnet from source.",
    call. = FALSE
  )
}

.onAttach <- function(libname, pkgname) {
  missing <- .bvarnet_missing_models()
  if (length(missing) == 0L) {
    return(invisible())
  }
  packageStartupMessage(
    "bvarnet: the precompiled Stan models are not available ",
    "(missing: ", paste(missing, collapse = ", "), ").\n",
    "This happens when the package is installed from a CRAN/r-universe binary, ",
    "or from source without CmdStan set up, so the models were never compiled.\n",
    "Calls to bvar() will fail until the models are compiled.\n",
    "See the installation instructions at ",
    "https://github.com/flo1met/bvarnet#installation ",
    "for how to install cmdstanr + CmdStan and reinstall bvarnet from source."
  )
}
