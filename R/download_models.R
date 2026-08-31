## Precompiled Stan model download.
##
## Two TEST-ONLY seams (never user-facing knobs):
## - `BVARNET_MODELS_BASE_URL` lets tests/humans point the client at a local server
##   instead of the real GitHub release host. Only the HTTPS requirement is relaxed
##   when it is set, since the seam is only ever used for http://localhost.
## - `BVARNET_MODELS_TAG` overrides the release-tag path segment (default
##   `v<packageVersion>`), so CI can test against a throwaway pre-release tag like
##   `citest-<sha>` without ever creating a real `v<version>` release. Without this,
##   the "test on a pre-release tag" flow is impossible: GitHub serves assets at
##   releases/download/<tag>/..., and the client would otherwise always look under
##   the v<version> tag, which deliberately must not exist before the real release.
## Every safety check (SHA-256 of the asset, source_hash vs the installed .stan
## sources) applies unconditionally regardless of either seam.

## ---- constants ----

.bvarnet_stan_model_names <- function() c("model_gaussian", "model_binary", "model_ordinal")

.bvarnet_exe_ext <- function() if (.Platform$OS.type == "windows") ".exe" else ""

## ---- platform key ----

#' @noRd
.bvarnet_platform_key <- function(sysname = Sys.info()[["sysname"]],
                                  arch = R.version$arch) {
  os <- switch(sysname,
    Darwin  = "macos",
    Linux   = "linux",
    Windows = "windows",
    stop("bvarnet: unsupported platform: ", sysname, call. = FALSE)
  )
  arch_key <- switch(arch,
    aarch64 = "arm64",
    arm64   = "arm64",
    x86_64  = "x86_64",
    stop("bvarnet: unsupported CPU architecture: ", arch, call. = FALSE)
  )
  paste0(os, "-", arch_key)
}

## ---- manifest URL + fetch ----

#' @noRd
.bvarnet_models_base_url <- function() {
  Sys.getenv("BVARNET_MODELS_BASE_URL",
             unset = "https://github.com/flo1met/bvarnet/releases/download")
}

#' @noRd
.bvarnet_assert_https <- function(url) {
  overridden <- nzchar(Sys.getenv("BVARNET_MODELS_BASE_URL", unset = ""))
  if (!overridden && !startsWith(url, "https://")) {
    stop("bvarnet: refusing a non-HTTPS download URL: ", url, call. = FALSE)
  }
  invisible(TRUE)
}

#' Release-tag path segment of the download URL. Defaults to `v<version>` (the
#' documented release-tag rule); `BVARNET_MODELS_TAG` is the test-only override
#' that lets CI fetch from a throwaway pre-release tag (see file header).
#' @noRd
.bvarnet_release_tag <- function(version = as.character(utils::packageVersion("bvarnet"))) {
  Sys.getenv("BVARNET_MODELS_TAG", unset = paste0("v", version))
}

#' @noRd
.bvarnet_manifest_url <- function(version = as.character(utils::packageVersion("bvarnet"))) {
  paste0(.bvarnet_models_base_url(), "/", .bvarnet_release_tag(version), "/manifest.json")
}

#' Fetch the per-version manifest.
#'
#' Three distinct outcomes, per the plan's client resolution rules:
#'   - no release / no manifest for this version -> NULL (caller falls back silently)
#'   - manifest present but not valid JSON -> hard error (unexpected, worth surfacing)
#'   - manifest present and parses -> the parsed list
#' @noRd
.bvarnet_fetch_manifest <- function(version = as.character(utils::packageVersion("bvarnet")),
                                    quiet = TRUE) {
  url <- .bvarnet_manifest_url(version)
  .bvarnet_assert_https(url)
  dest <- tempfile(fileext = ".json")
  on.exit(unlink(dest), add = TRUE)
  ok <- tryCatch({
    utils::download.file(url, dest, method = "libcurl", mode = "wb", quiet = quiet)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok) || !file.exists(dest) || file.size(dest) == 0L) {
    return(NULL)  # no release/manifest for this version -> compile-from-source fallback
  }
  manifest <- tryCatch(jsonlite::fromJSON(dest, simplifyVector = TRUE),
                       error = function(e) NULL)
  if (is.null(manifest)) {
    stop(
      "bvarnet: the manifest at ", url, " could not be parsed as JSON. ",
      "This is unexpected for a published release and may indicate a problem ",
      "with the release itself -- please report this at ",
      "https://github.com/flo1met/bvarnet/issues.",
      call. = FALSE
    )
  }
  manifest
}

## ---- source hash ----

#' @noRd
.bvarnet_stan_source_files <- function(package = "bvarnet") {
  names <- sort(c("model_gaussian.stan", "model_binary.stan", "model_ordinal.stan"))
  vapply(names, function(f) {
    system.file(file.path("bin", "stan", f), package = package, mustWork = TRUE)
  }, character(1))
}

#' Read a .stan file as raw bytes with CR stripped (CRLF/LF normalization).
#' @noRd
.bvarnet_read_stan_bytes_lf <- function(path) {
  raw <- readBin(path, what = "raw", n = file.size(path))
  raw[raw != as.raw(0x0d)]
}

#' SHA-256 over the concatenated, line-ending-normalized model_*.stan sources.
#' Must match the build's hasher byte-for-byte.
#' @noRd
.bvarnet_source_hash <- function(package = "bvarnet") {
  files <- .bvarnet_stan_source_files(package)
  bytes <- do.call(c, lapply(files, .bvarnet_read_stan_bytes_lf))
  digest::digest(bytes, algo = "sha256", serialize = FALSE)
}

## ---- cache dir / key ----

#' @noRd
.bvarnet_cache_root <- function() tools::R_user_dir("bvarnet", which = "cache")

#' @noRd
.bvarnet_cache_key <- function() {
  hash8 <- substr(.bvarnet_source_hash(), 1, 8)
  paste0(as.character(utils::packageVersion("bvarnet")), "-", hash8)
}

#' @noRd
.bvarnet_models_cache_dir <- function() {
  file.path(.bvarnet_cache_root(), .bvarnet_cache_key())
}

#' @noRd
.bvarnet_runtime_dir <- function() {
  file.path(.bvarnet_models_cache_dir(), "runtime")
}

## ---- cache integrity ----

#' Is the cache dir complete and verified? Missing/mismatched checksum ==
#' absent, never run. The completion marker is written last at staging time,
#' so its presence (with matching checksums) is the whole signal.
#' @noRd
.bvarnet_cache_complete <- function(dir) {
  marker <- file.path(dir, "_complete.json")
  if (!file.exists(marker)) return(FALSE)
  info <- tryCatch(jsonlite::fromJSON(marker, simplifyVector = TRUE), error = function(e) NULL)
  if (is.null(info) || is.null(info$sha256)) return(FALSE)
  ext <- .bvarnet_exe_ext()
  for (name in names(info$sha256)) {
    exe <- file.path(dir, paste0(name, ext))
    if (!file.exists(exe)) return(FALSE)
    got <- tryCatch(unname(digest::digest(exe, algo = "sha256", file = TRUE)),
                    error = function(e) NA_character_)
    if (!identical(got, unname(info$sha256[[name]]))) return(FALSE)
  }
  TRUE
}

#' Resolver-facing accessor: a verified cache entry for `model_name`, or NULL.
#' @noRd
.bvarnet_cache_model <- function(model_name) {
  dir <- .bvarnet_models_cache_dir()
  if (!.bvarnet_cache_complete(dir)) return(NULL)
  ext  <- .bvarnet_exe_ext()
  exe  <- file.path(dir, paste0(model_name, ext))
  stan <- file.path(dir, paste0(model_name, ".stan"))
  if (!file.exists(exe) || !file.exists(stan)) return(NULL)
  list(exe = exe, stan = stan, dir = dir)
}

## ---- cmdstan_valid guard (mirrors instantiate:::cmdstan_valid) ----

#' @noRd
.bvarnet_cmdstan_valid <- function(path) {
  length(path) > 0L && nzchar(path) && dir.exists(path) &&
    file.exists(file.path(path, "makefile")) &&
    any(grepl("^CMDSTAN_VERSION",
              tryCatch(readLines(file.path(path, "makefile")), error = function(e) character(0))))
}

## ---- loader-path helper ----

#' Platform env var the dynamic loader reads for extra search paths.
#' NA on Windows: TBB is resolved by co-locating the DLL beside the exe instead.
#' @noRd
.bvarnet_libpath_var <- function() {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Darwin") "DYLD_FALLBACK_LIBRARY_PATH"
  else if (sys == "Windows") NA_character_
  else "LD_LIBRARY_PATH"
}

#' Transiently PREPEND `dir` to the platform loader-path env var around `code`,
#' restoring the exact prior state on exit (including truly-unset, via
#' `Sys.getenv(unset = NA)` -- restoring to `""` would be a subtly different,
#' wrong state). Prepend, not replace: R itself sets DYLD_FALLBACK_LIBRARY_PATH
#' on macOS, and replacing it would clobber R's own fallback paths for the
#' duration of the call.
#' @noRd
.bvarnet_with_libpath <- function(dir, code) {
  var <- .bvarnet_libpath_var()
  if (is.na(var)) return(force(code))
  old <- Sys.getenv(var, unset = NA)
  new <- if (is.na(old) || !nzchar(old)) dir else paste(dir, old, sep = .Platform$path.sep)
  on.exit(
    if (is.na(old)) Sys.unsetenv(var) else do.call(Sys.setenv, stats::setNames(list(old), var)),
    add = TRUE
  )
  do.call(Sys.setenv, stats::setNames(list(new), var))
  force(code)
}

## ---- unified resolver ----

#' Which bundled Stan model backs a given outcome family.
#' Shared by both `bvar()` call sites.
#' @noRd
.bvarnet_model_for_family <- function(family) {
  switch(family,
    bernoulli = "model_binary",
    ordinal   = "model_ordinal",
    gaussian  = "model_gaussian",
    stop("Unknown family: ", family, call. = FALSE)
  )
}

#' Resolve a ready-to-sample CmdStanModel for `model_name`.
#'
#' Resolution order: install-tree exe (source / CRAN-with-CmdStan
#' installs) -> verified cache exe (downloaded or compiled-to-cache) -> error.
#'
#' Returns `list(model = <CmdStanModel>, cache_dir = <path or NULL>)`.
#' `cache_dir` is NULL when resolved from the install tree (no loader-env
#' wrapping needed -- cmdstanr's own rpath-based resolution applies); it is the
#' cache directory when resolved from the cache, so the caller can wrap the
#' `$sample()` call in `.bvarnet_with_libpath(file.path(cache_dir, "lib",
#' "tbb"), ...)` to make the bundled TBB findable. This wrapping is
#' deliberately done by the CALLER around `$sample()`, not here: the loader env
#' must be live when `$sample()` spawns the child, which is after this
#' function has already returned.
#' @noRd
.bvarnet_stan_model <- function(model_name) {
  ext <- .bvarnet_exe_ext()
  exe <- system.file(file.path("bin", "stan", paste0(model_name, ext)), package = "bvarnet")
  if (nzchar(exe) && file.exists(exe)) {
    return(list(
      model = instantiate::stan_package_model(name = model_name, package = "bvarnet"),
      cache_dir = NULL
    ))
  }

  cached <- .bvarnet_cache_model(model_name)
  if (!is.null(cached)) {
    if (!requireNamespace("cmdstanr", quietly = TRUE)) {
      .bvarnet_models_missing_error(model_name)
    }
    # Mirror instantiate::stan_package_model()'s own idiom exactly: capture the
    # prior path and restore it on return, but ONLY if it was already a valid
    # CmdStan. This guard is what makes the restore-on-return safe for BOTH
    # cases: a download-only user has no valid prior path, so nothing is
    # restored and the stub stays live through $sample(); a user with a real
    # CmdStan gets it back, with no leaked state.
    path_old <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
    if (.bvarnet_cmdstan_valid(path_old)) {
      on.exit(suppressMessages(cmdstanr::set_cmdstan_path(path_old)), add = TRUE)
    }
    suppressMessages(cmdstanr::set_cmdstan_path(.bvarnet_runtime_dir()))
    model <- cmdstanr::cmdstan_model(stan_file = cached$stan, exe_file = cached$exe,
                                     compile = FALSE)  # never $compile() -- OQ1 invariant
    return(list(model = model, cache_dir = cached$dir))
  }

  .bvarnet_models_missing_error(model_name)
}

## ---- cmdstanr install-on-demand + toolchain probe ----

#' @noRd
.bvarnet_ensure_cmdstanr <- function() {
  if (requireNamespace("cmdstanr", quietly = TRUE)) return(invisible(TRUE))
  utils::install.packages(
    "cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
  )
  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "bvarnet: failed to install the 'cmdstanr' package, which is required to ",
      "run Stan models. Please install it manually:\n",
      "install.packages(\"cmdstanr\", repos = c(\"https://mc-stan.org/r-packages/\", ",
      "getOption(\"repos\")))",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' @noRd
.bvarnet_toolchain_ok <- function() {
  if (!requireNamespace("cmdstanr", quietly = TRUE)) return(FALSE)
  tryCatch({
    cmdstanr::check_cmdstan_toolchain(quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)
}

## ---- consent ----

#' @noRd
.bvarnet_consent_download <- function(asset_name, ask = interactive()) {
  if (!isTRUE(ask)) {
    return(
      isTRUE(getOption("bvarnet.allow_download")) ||
        tolower(Sys.getenv("BVARNET_ALLOW_DOWNLOAD", "")) %in% c("1", "true", "yes")
    )
  }
  if (!interactive()) return(FALSE)
  cat(sprintf(
    "bvarnet will download a precompiled Stan executable (%s) from\n%s\nand run it on your machine.\n",
    asset_name, .bvarnet_models_base_url()
  ))
  identical(utils::menu(c("Yes", "No"), title = "Proceed with download?"), 1L)
}

## ---- staging a downloaded asset ----

#' Verify the inner SHA256SUMS (if present), set exec bits, strip macOS
#' quarantine (belt-and-suspenders), co-locate tbb.dll on Windows,
#' and write the completion marker LAST. `dir` is a private per-process temp
#' dir under the cache root; the caller atomically renames it into place.
#' @noRd
.bvarnet_stage_extracted_asset <- function(dir) {
  names <- .bvarnet_stan_model_names()
  ext   <- .bvarnet_exe_ext()

  sums_file <- file.path(dir, "SHA256SUMS")
  if (file.exists(sums_file)) {
    for (line in readLines(sums_file)) {
      parts <- strsplit(trimws(line), "\\s+")[[1]]
      if (length(parts) < 2L) next
      sha <- parts[1]
      rel <- parts[length(parts)]
      f <- file.path(dir, rel)
      if (!file.exists(f)) {
        stop("bvarnet: asset missing file listed in SHA256SUMS: ", rel, call. = FALSE)
      }
      got <- unname(digest::digest(f, algo = "sha256", file = TRUE))
      if (!identical(got, sha)) {
        stop("bvarnet: asset file failed inner SHA-256 verification: ", rel, call. = FALSE)
      }
    }
  }

  exe_paths <- file.path(dir, paste0(names, ext))
  for (exe in exe_paths) {
    if (!file.exists(exe)) {
      stop("bvarnet: asset missing expected model executable: ", basename(exe), call. = FALSE)
    }
    Sys.chmod(exe, "0755")
  }
  stanc <- file.path(dir, "runtime", "bin", paste0("stanc", ext))
  if (file.exists(stanc)) Sys.chmod(stanc, "0755")

  sys <- Sys.info()[["sysname"]]
  if (sys == "Darwin") {
    for (exe in exe_paths) {
      suppressWarnings(system2("xattr", c("-d", "com.apple.quarantine", shQuote(exe)),
                               stdout = FALSE, stderr = FALSE))
    }
  }
  if (sys == "Windows") {
    tbb_dlls <- Sys.glob(file.path(dir, "lib", "tbb", "*.dll"))
    for (exe in exe_paths) file.copy(tbb_dlls, dirname(exe), overwrite = TRUE)
  }

  .bvarnet_write_completion_marker(dir, exe_paths, names,
                                   extra = list(platform = .bvarnet_platform_key(),
                                               downloaded_at = format(Sys.time(), tz = "UTC", usetz = TRUE)))
  dir
}

#' @noRd
.bvarnet_write_completion_marker <- function(dir, exe_paths, names, extra = list()) {
  sha <- stats::setNames(
    as.list(vapply(exe_paths, function(f) unname(digest::digest(f, algo = "sha256", file = TRUE)),
                   character(1))),
    names
  )
  marker <- c(list(sha256 = sha), extra)
  jsonlite::write_json(marker, file.path(dir, "_complete.json"), auto_unbox = TRUE)
  invisible(TRUE)
}

## ---- download path ----

#' @noRd
.bvarnet_download_models <- function(ask = interactive(), quiet = FALSE) {
  .bvarnet_ensure_cmdstanr()
  version  <- as.character(utils::packageVersion("bvarnet"))
  manifest <- .bvarnet_fetch_manifest(version)
  if (is.null(manifest)) {
    message(
      "bvarnet: no precompiled-model release found for version ", version, ". ",
      "Falling back to compile-from-source; see ?bvarnet_setup_models."
    )
    return(invisible(FALSE))
  }
  if (isTRUE(manifest$revoked)) {
    message(
      "bvarnet: the precompiled binaries for version ", version, " were revoked (",
      manifest$revoked_reason %||% "no reason given", "). ",
      "Falling back to compile-from-source; see ?bvarnet_setup_models."
    )
    return(invisible(FALSE))
  }

  local_hash    <- .bvarnet_source_hash()
  manifest_hash <- sub("^sha256:", "", manifest$source_hash %||% "")
  if (!identical(local_hash, manifest_hash)) {
    stop(
      "bvarnet: the installed package's Stan sources do not match the published ",
      "manifest for version ", version, " (source hash mismatch). This can happen ",
      "with a development install between releases. Refusing to download a ",
      "possibly mismatched binary; use bvarnet_setup_models(method = \"compile\") instead.",
      call. = FALSE
    )
  }

  plat  <- .bvarnet_platform_key()
  entry <- manifest$platforms[[plat]]
  if (is.null(entry)) {
    message(
      "bvarnet: no precompiled binary published for platform \"", plat, "\" (version ",
      version, "). Falling back to compile-from-source."
    )
    return(invisible(FALSE))
  }
  asset_name <- entry$asset
  asset_sha  <- entry$sha256
  asset_url  <- paste0(.bvarnet_models_base_url(), "/", .bvarnet_release_tag(version), "/", asset_name)
  .bvarnet_assert_https(asset_url)

  if (!.bvarnet_consent_download(asset_name, ask = ask)) {
    message("bvarnet: download declined; models were not set up. See ?bvarnet_setup_models.")
    return(invisible(FALSE))
  }

  cache_root <- .bvarnet_cache_root()
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  work <- tempfile("bvarnet_dl_", tmpdir = cache_root)
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

  tar_path <- file.path(work, asset_name)
  ok <- tryCatch({
    utils::download.file(asset_url, tar_path, method = "libcurl", mode = "wb", quiet = quiet)
    TRUE
  }, error = function(e) FALSE, warning = function(w) FALSE)
  if (!isTRUE(ok) || !file.exists(tar_path)) {
    stop("bvarnet: failed to download ", asset_url, call. = FALSE)
  }

  actual_sha <- unname(digest::digest(tar_path, algo = "sha256", file = TRUE))
  if (!identical(actual_sha, asset_sha)) {
    unlink(tar_path)
    stop(
      "bvarnet: downloaded asset failed SHA-256 verification (expected ", asset_sha,
      ", got ", actual_sha, "). Aborting -- the file was deleted.",
      call. = FALSE
    )
  }

  extract_dir <- file.path(work, "extracted")
  dir.create(extract_dir)
  utils::untar(tar_path, exdir = extract_dir)

  staged <- .bvarnet_stage_extracted_asset(extract_dir)

  final_dir <- .bvarnet_models_cache_dir()
  if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE, force = TRUE)
  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  file.rename(staged, final_dir)  # atomic on the same filesystem (both under cache_root)

  if (!quiet) message("bvarnet: Stan models set up at ", final_dir)
  invisible(TRUE)
}

## ---- compile-to-cache path ----

#' @noRd
.bvarnet_compile_models <- function(ask = interactive(), quiet = FALSE) {
  .bvarnet_ensure_cmdstanr()
  if (!.bvarnet_toolchain_ok()) {
    stop(
      "bvarnet: no working C++ toolchain was detected. Install one and ",
      "re-run bvarnet_setup_models(method = \"compile\"). See ?bvarnet_setup_models.",
      call. = FALSE
    )
  }
  if (isTRUE(ask) && interactive()) {
    cat("Compiling the Stan models can take a couple of minutes.\n")
    if (!identical(utils::menu(c("Yes", "No"), title = "Proceed?"), 1L)) {
      message("bvarnet: compilation declined; models were not set up.")
      return(invisible(FALSE))
    }
  }

  cs_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!.bvarnet_cmdstan_valid(cs_path)) {
    cmdstanr::install_cmdstan(quiet = quiet)
    cs_path <- cmdstanr::cmdstan_path()
  }

  cache_root <- .bvarnet_cache_root()
  dir.create(cache_root, recursive = TRUE, showWarnings = FALSE)
  work <- tempfile("bvarnet_compile_", tmpdir = cache_root)
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)

  names <- .bvarnet_stan_model_names()
  ext   <- .bvarnet_exe_ext()
  for (name in names) {
    stan_file <- system.file(file.path("bin", "stan", paste0(name, ".stan")), package = "bvarnet")
    exe_file  <- file.path(work, paste0(name, ext))
    cmdstanr::cmdstan_model(stan_file = stan_file, exe_file = exe_file, compile = TRUE,
                            quiet = quiet)
    file.copy(stan_file, file.path(work, paste0(name, ".stan")))
  }

  runtime_bin <- file.path(work, "runtime", "bin")
  dir.create(runtime_bin, recursive = TRUE, showWarnings = FALSE)
  file.copy(file.path(cs_path, "makefile"), file.path(work, "runtime", "makefile"), overwrite = TRUE)
  file.copy(file.path(cs_path, "bin", paste0("stanc", ext)),
           file.path(runtime_bin, paste0("stanc", ext)), overwrite = TRUE)
  Sys.chmod(file.path(runtime_bin, paste0("stanc", ext)), "0755")

  exe_paths <- file.path(work, paste0(names, ext))
  .bvarnet_write_completion_marker(work, exe_paths, names,
                                   extra = list(platform = .bvarnet_platform_key(),
                                               compiled_at = format(Sys.time(), tz = "UTC", usetz = TRUE)))

  final_dir <- .bvarnet_models_cache_dir()
  if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE, force = TRUE)
  dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
  file.rename(work, final_dir)

  if (!quiet) message("bvarnet: Stan models compiled and cached at ", final_dir)
  invisible(TRUE)
}

## ---- interactive menu ----

#' @noRd
.bvarnet_prompt_method <- function() {
  manifest <- tryCatch(.bvarnet_fetch_manifest(), error = function(e) NULL)
  plat     <- tryCatch(.bvarnet_platform_key(), error = function(e) NA_character_)
  download_available <-
    !is.null(manifest) && !isTRUE(manifest$revoked) && !is.na(plat) &&
    !is.null(manifest$platforms[[plat]]) &&
    identical(sub("^sha256:", "", manifest$source_hash %||% ""), .bvarnet_source_hash())

  choices <- character(0)
  actions <- character(0)
  if (download_available) {
    choices <- c(choices, "Download precompiled binaries (toolchain-free)")
    actions <- c(actions, "download")
  }
  choices <- c(choices, "Compile locally (requires a C++ toolchain)", "Cancel")
  actions <- c(actions, "compile", "cancel")

  sel <- utils::menu(choices, title = "bvarnet: set up Stan models")
  if (sel == 0L || identical(actions[sel], "cancel")) return(NULL)
  actions[sel]
}

## ---- exported entry points ----

#' Set up bvarnet's precompiled Stan models
#'
#' `bvar()` requires compiled Stan model executables. When installing from a
#' CRAN/r-universe binary (no CmdStan at install time) they are not built
#' automatically, so this function sets them up: either by downloading
#' precompiled, toolchain-free binaries for your platform, or by compiling
#' them locally with your own CmdStan installation. Either way the result is
#' cached in a per-user cache directory (see [bvarnet_model_cache_dir()]) and
#' persists across sessions.
#'
#' @param method `"download"` or `"compile"`. If omitted, and the session is
#'   interactive, you are prompted to choose.
#' @param force Logical. If `TRUE`, clear the current cache entry first and
#'   set up from scratch.
#' @param ask Logical. Whether to prompt for confirmation before downloading
#'   or compiling. Defaults to `interactive()`. When `FALSE`, a download is
#'   only performed if pre-authorised via `options(bvarnet.allow_download =
#'   TRUE)` or the `BVARNET_ALLOW_DOWNLOAD` environment variable.
#'
#' @details
#' Downloaded binaries are only ever fetched over HTTPS from this package's
#' official GitHub repository, and are verified against a published SHA-256
#' checksum before being made executable; a mismatch aborts and deletes the
#' download. The download is also pinned to a hash of your installed
#' package's Stan sources, so a binary can never silently be paired with the
#' wrong model version.
#'
#' `bvar()` itself never downloads or prompts -- all network activity happens
#' here, in this explicitly user-invoked function, so scripts and
#' non-interactive sessions never trigger a download by surprise.
#'
#' @seealso [bvarnet_model_cache_dir()] to inspect where models are cached,
#'   [bvarnet_clear_model_cache()] to remove a cached entry, and [bvar()] for
#'   the function that uses the resulting models.
#'
#' @return Invisibly, `TRUE` on success, `FALSE` if declined or unavailable.
#' @examples
#' \dontrun{
#' # Interactive: prompts you to choose download vs. compile
#' bvarnet_setup_models()
#'
#' # Non-interactive / scripted, pre-authorised for download:
#' options(bvarnet.allow_download = TRUE)
#' bvarnet_setup_models(method = "download", ask = FALSE)
#'
#' # Force a clean re-setup (e.g. after a corrupted cache)
#' bvarnet_setup_models(method = "download", force = TRUE)
#' }
#' @export
bvarnet_setup_models <- function(method = c("download", "compile"),
                                 force = FALSE, ask = interactive()) {
  explicit_method <- !missing(method)
  method <- if (explicit_method) match.arg(method, c("download", "compile")) else NULL

  if (isTRUE(force)) bvarnet_clear_model_cache(quiet = TRUE)

  if (is.null(method)) {
    if (!isTRUE(ask)) {
      stop(
        "bvarnet: bvarnet_setup_models() needs `method` to be specified ",
        "(\"download\" or \"compile\") when `ask = FALSE` or in a non-interactive session.",
        call. = FALSE
      )
    }
    method <- .bvarnet_prompt_method()
    if (is.null(method)) return(invisible(FALSE))
  }

  if (identical(method, "download")) {
    .bvarnet_download_models(ask = ask)
  } else {
    .bvarnet_compile_models(ask = ask)
  }
}

#' Remove bvarnet's cached Stan model binaries
#'
#' @param all_versions Logical. Remove the entire cache (all versions), not
#'   just the entry for the currently installed version.
#' @param quiet Logical. Suppress the confirmation message.
#' @return Invisibly, the path removed.
#' @seealso [bvarnet_setup_models()] to set models up again afterward,
#'   [bvarnet_model_cache_dir()] to inspect the cache path.
#' @export
bvarnet_clear_model_cache <- function(all_versions = FALSE, quiet = FALSE) {
  root   <- .bvarnet_cache_root()
  target <- if (isTRUE(all_versions)) root else .bvarnet_models_cache_dir()
  if (dir.exists(target)) {
    unlink(target, recursive = TRUE, force = TRUE)
    if (!quiet) message("bvarnet: removed model cache at ", target)
  }
  invisible(target)
}

#' Path to bvarnet's Stan model cache directory
#'
#' @return Character scalar: the cache directory for the currently installed
#'   package version (may not exist yet).
#' @seealso [bvarnet_setup_models()] to populate this directory,
#'   [bvarnet_clear_model_cache()] to remove it.
#' @export
bvarnet_model_cache_dir <- function() {
  .bvarnet_models_cache_dir()
}
