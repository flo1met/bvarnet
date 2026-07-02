#!/usr/bin/env Rscript
# build-manifest.R
#
# Aggregates each platform's asset-info.json (written by assemble-model-asset.R)
# into the single per-version manifest.json release asset (dev/precompiled_model_download_plan.md,
# Section A manifest decision + schema).
#
# Asserts all platforms agree on `version`, `source_hash`, and `cmdstan_version` --
# this IS the CI-side "verify end-to-end" check from Section A step 4: since every
# platform computes source_hash via the SAME bvarnet:::.bvarnet_source_hash() (which
# strips CR before hashing), a Windows runner's CRLF checkout must still produce an
# identical hash to macOS/Linux. If it doesn't, something is wrong and the build
# must fail loudly rather than publish a manifest with a silently platform-specific hash.

build_manifest <- function(artifacts_dir, out_file = "manifest.json") {
  stopifnot(requireNamespace("jsonlite", quietly = TRUE))

  info_files <- list.files(artifacts_dir, pattern = "^asset-info\\.json$",
                           recursive = TRUE, full.names = TRUE)
  if (length(info_files) == 0L) {
    stop("build_manifest: no asset-info.json files found under ", artifacts_dir, call. = FALSE)
  }
  infos <- lapply(info_files, jsonlite::fromJSON)

  field <- function(name) vapply(infos, `[[`, character(1), name)

  versions <- unique(field("version"))
  if (length(versions) != 1L) {
    stop("build_manifest: platforms disagree on package version: ",
         paste(versions, collapse = ", "), call. = FALSE)
  }
  hashes <- unique(field("source_hash"))
  if (length(hashes) != 1L) {
    stop("build_manifest: platforms disagree on source_hash (CRLF/LF normalization ",
         "regression?): ", paste(hashes, collapse = ", "), call. = FALSE)
  }
  cmdstan_versions <- unique(field("cmdstan_version"))
  if (length(cmdstan_versions) != 1L) {
    stop("build_manifest: platforms disagree on cmdstan_version: ",
         paste(cmdstan_versions, collapse = ", "), call. = FALSE)
  }

  platform_keys <- field("platform")
  if (anyDuplicated(platform_keys)) {
    stop("build_manifest: duplicate platform key(s): ",
         paste(platform_keys[duplicated(platform_keys)], collapse = ", "), call. = FALSE)
  }

  platforms <- stats::setNames(
    lapply(infos, function(i) list(asset = i$asset, sha256 = i$sha256)),
    platform_keys
  )

  manifest <- list(
    version         = versions,
    source_hash     = hashes,
    cmdstan_version = cmdstan_versions,
    revoked         = FALSE,
    revoked_reason  = NA,  # reserved for the recall escape valve; NA -> JSON null
    platforms       = platforms
  )
  jsonlite::write_json(manifest, out_file, auto_unbox = TRUE, null = "null", pretty = TRUE,
                      na = "null")
  message("build_manifest: wrote ", out_file, " covering platforms: ",
          paste(platform_keys, collapse = ", "))
  invisible(manifest)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  build_manifest(
    artifacts_dir = args[1],
    out_file      = if (length(args) >= 2L) args[2] else "manifest.json"
  )
}
