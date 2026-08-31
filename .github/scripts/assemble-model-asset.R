#!/usr/bin/env Rscript
# assemble-model-asset.R
#
# Build-side counterpart to R/download_models.R. Run AFTER `bvarnet` has been installed 
# from source on a runner with CmdStan set up (so bin/stan/model_*[.exe] exist), and 
# after Windows static-linking has already been applied to CmdStan's make/local, if applicable.
#
# Produces, under `out_dir`:
#   staging/                      -- the exact asset layout (tarball root)
#     model_gaussian[.exe]        model_gaussian.stan
#     model_binary[.exe]          model_binary.stan
#     model_ordinal[.exe]         model_ordinal.stan
#     runtime/makefile            runtime/bin/stanc[.exe]
#     lib/tbb/<tbb libs>
#     SHA256SUMS
#   bvarnet-models-<version>-<platform>.tar.gz   -- the actual release asset
#   asset-info.json                              -- side-car for the manifest job:
#                                                    {platform, asset, sha256, source_hash,
#                                                     cmdstan_version, version}
#
# Deliberately reuses bvarnet:::.bvarnet_source_hash() (the SAME function the client calls)
# so the build and the client can never disagree about what "the source hash" means.

assemble_model_asset <- function(platform, cmdstan_version, out_dir = "asset-out") {
  stopifnot(requireNamespace("bvarnet", quietly = TRUE))
  stopifnot(requireNamespace("digest", quietly = TRUE))
  stopifnot(requireNamespace("jsonlite", quietly = TRUE))

  is_windows <- .Platform$OS.type == "windows"
  ext        <- if (is_windows) ".exe" else ""
  version    <- as.character(utils::packageVersion("bvarnet"))
  model_names <- c("model_gaussian", "model_binary", "model_ordinal")

  bin_stan <- system.file("bin", "stan", package = "bvarnet")
  stopifnot(dir.exists(bin_stan))

  staging <- file.path(out_dir, "staging")
  unlink(staging, recursive = TRUE, force = TRUE)
  dir.create(file.path(staging, "runtime", "bin"), recursive = TRUE)
  dir.create(file.path(staging, "lib", "tbb"), recursive = TRUE)

  # -- model executables + their .stan sources ---------------------------------
  for (name in model_names) {
    exe  <- file.path(bin_stan, paste0(name, ext))
    stan <- file.path(bin_stan, paste0(name, ".stan"))
    stopifnot(file.exists(exe), file.exists(stan))
    file.copy(exe, file.path(staging, paste0(name, ext)))
    file.copy(stan, file.path(staging, paste0(name, ".stan")))
  }

  # -- runtime stub: one-line makefile + standalone stanc ----------------
  # bin/stansummary and bin/diagnose are deliberately NOT bundled -- bvar() only
  # ever calls $sample()/$draws()/$time()/$metadata(), never those binaries.
  cmdstan_dir <- cmdstanr::cmdstan_path()
  writeLines(paste0("CMDSTAN_VERSION := ", cmdstan_version),
             file.path(staging, "runtime", "makefile"))
  stanc_src <- file.path(cmdstan_dir, "bin", paste0("stanc", ext))
  stopifnot(file.exists(stanc_src))
  file.copy(stanc_src, file.path(staging, "runtime", "bin", paste0("stanc", ext)))

  # -- bundled TBB runtime libs -- uniform layout on every platform ------
  tbb_dir <- file.path(cmdstan_dir, "stan", "lib", "stan_math", "lib", "tbb")
  stopifnot(dir.exists(tbb_dir))
  tbb_pattern <- if (is_windows) "\\.dll$" else if (Sys.info()[["sysname"]] == "Darwin") "\\.dylib$" else "\\.so(\\.[0-9]+)*$"
  tbb_libs <- list.files(tbb_dir, pattern = tbb_pattern, full.names = TRUE)
  stopifnot(length(tbb_libs) > 0L)
  file.copy(tbb_libs, file.path(staging, "lib", "tbb"))

  # -- SHA256SUMS over every staged file (client verifies before exec'ing) -----
  all_files <- list.files(staging, recursive = TRUE, full.names = FALSE)
  sums <- vapply(all_files, function(f) {
    unname(digest::digest(file.path(staging, f), algo = "sha256", file = TRUE))
  }, character(1))
  writeLines(paste(sums, all_files), file.path(staging, "SHA256SUMS"))

  # -- tarball ------------------------------------------------------------------
  asset_name <- sprintf("bvarnet-models-%s-%s.tar.gz", version, platform)
  asset_path <- file.path(out_dir, asset_name)
  old_wd <- setwd(staging)
  on.exit(setwd(old_wd), add = TRUE)
  utils::tar(file.path("..", basename(asset_name)), files = ".", compression = "gzip",
             tar = "internal")
  setwd(old_wd)

  # -- side-car info for the manifest job ---------------------------------------
  info <- list(
    platform        = platform,
    version         = version,
    asset           = asset_name,
    sha256          = unname(digest::digest(asset_path, algo = "sha256", file = TRUE)),
    source_hash     = paste0("sha256:", bvarnet:::.bvarnet_source_hash()),
    cmdstan_version = cmdstan_version
  )
  jsonlite::write_json(info, file.path(out_dir, "asset-info.json"), auto_unbox = TRUE)
  invisible(info)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  assemble_model_asset(
    platform        = args[1],
    cmdstan_version = args[2],
    out_dir         = if (length(args) >= 3L) args[3] else "asset-out"
  )
}
