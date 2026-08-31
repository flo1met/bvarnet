# test-download-models.R — unit tests for the precompiled-model download
# fallback (R/download_models.R, R/zzz.R). All mocked/local-filesystem; no
# network access, safe on CRAN.

# ═══════════════════════════════════════════════════════════════════════════
# §1 platform key
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_platform_key() maps sysname+arch to platform strings", {
  expect_equal(.bvarnet_platform_key("Darwin", "aarch64"), "macos-arm64")
  expect_equal(.bvarnet_platform_key("Darwin", "arm64"), "macos-arm64")
  expect_equal(.bvarnet_platform_key("Darwin", "x86_64"), "macos-x86_64") # incl. Rosetta
  expect_equal(.bvarnet_platform_key("Linux", "x86_64"), "linux-x86_64")
  expect_equal(.bvarnet_platform_key("Windows", "x86_64"), "windows-x86_64")
})

test_that(".bvarnet_platform_key() errors on unsupported platform/arch", {
  expect_error(.bvarnet_platform_key("SunOS", "x86_64"), "unsupported platform")
  expect_error(.bvarnet_platform_key("Linux", "mips"), "unsupported CPU architecture")
})

# ═══════════════════════════════════════════════════════════════════════════
# §2 manifest URL / base URL / HTTPS assertion
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_manifest_url() uses the installed version, never \"latest\"", {
  url <- .bvarnet_manifest_url(version = "1.0.1")
  expect_equal(url, "https://github.com/flo1met/bvarnet/releases/download/v1.0.1/manifest.json")
  expect_false(grepl("latest", url))
})

test_that(".bvarnet_manifest_url() defaults to utils::packageVersion", {
  url <- .bvarnet_manifest_url()
  expect_equal(url, .bvarnet_manifest_url(as.character(utils::packageVersion("bvarnet"))))
})

test_that("BVARNET_MODELS_TAG overrides the release-tag path segment (CI pre-release seam)", {
  old <- Sys.getenv("BVARNET_MODELS_TAG", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("BVARNET_MODELS_TAG")
          else Sys.setenv(BVARNET_MODELS_TAG = old), add = TRUE)

  Sys.unsetenv("BVARNET_MODELS_TAG")
  expect_equal(.bvarnet_release_tag("1.0.1"), "v1.0.1")
  expect_match(.bvarnet_manifest_url("1.0.1"), "/v1\\.0\\.1/manifest\\.json$")

  Sys.setenv(BVARNET_MODELS_TAG = "citest-abc1234")
  expect_equal(.bvarnet_release_tag("1.0.1"), "citest-abc1234")
  expect_match(.bvarnet_manifest_url("1.0.1"), "/citest-abc1234/manifest\\.json$")
  expect_false(grepl("v1\\.0\\.1", .bvarnet_manifest_url("1.0.1")))
})

test_that(".bvarnet_models_base_url() honours the BVARNET_MODELS_BASE_URL test seam", {
  old <- Sys.getenv("BVARNET_MODELS_BASE_URL", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("BVARNET_MODELS_BASE_URL")
          else Sys.setenv(BVARNET_MODELS_BASE_URL = old), add = TRUE)

  Sys.unsetenv("BVARNET_MODELS_BASE_URL")
  expect_equal(.bvarnet_models_base_url(), "https://github.com/flo1met/bvarnet/releases/download")

  Sys.setenv(BVARNET_MODELS_BASE_URL = "http://localhost:8000")
  expect_equal(.bvarnet_models_base_url(), "http://localhost:8000")
})

test_that(".bvarnet_assert_https() refuses non-HTTPS unless the test seam is set", {
  old <- Sys.getenv("BVARNET_MODELS_BASE_URL", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("BVARNET_MODELS_BASE_URL")
          else Sys.setenv(BVARNET_MODELS_BASE_URL = old), add = TRUE)

  Sys.unsetenv("BVARNET_MODELS_BASE_URL")
  expect_error(.bvarnet_assert_https("http://example.com/x"), "non-HTTPS")
  expect_true(.bvarnet_assert_https("https://example.com/x"))

  Sys.setenv(BVARNET_MODELS_BASE_URL = "http://localhost:8000")
  expect_true(.bvarnet_assert_https("http://localhost:8000/x"))
})

# ═══════════════════════════════════════════════════════════════════════════
# §3 source hash — reproducibility and CRLF/LF normalization
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_source_hash() is reproducible for the installed package", {
  h1 <- .bvarnet_source_hash()
  h2 <- .bvarnet_source_hash()
  expect_equal(h1, h2)
  expect_match(h1, "^[0-9a-f]{64}$")
})

test_that(".bvarnet_read_stan_bytes_lf() strips CR so CRLF and LF hash identically", {
  dir <- bvarnet_test_tempdir()
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)

  lf_path   <- file.path(dir, "lf.stan")
  crlf_path <- file.path(dir, "crlf.stan")
  content <- "parameters { real y; }\nmodel { y ~ std_normal(); }\n"
  writeBin(charToRaw(content), lf_path)
  writeBin(charToRaw(gsub("\n", "\r\n", content)), crlf_path)

  bytes_lf   <- .bvarnet_read_stan_bytes_lf(lf_path)
  bytes_crlf <- .bvarnet_read_stan_bytes_lf(crlf_path)
  expect_identical(bytes_lf, bytes_crlf)
  expect_equal(digest::digest(bytes_lf, algo = "sha256", serialize = FALSE),
              digest::digest(bytes_crlf, algo = "sha256", serialize = FALSE))
})

test_that(".bvarnet_stan_source_files() pins an explicit sorted list, never a glob", {
  files <- .bvarnet_stan_source_files()
  expect_equal(basename(files), c("model_binary.stan", "model_gaussian.stan", "model_ordinal.stan"))
  expect_false(any(grepl("functions\\.stan|\\.DS_Store", files)))
})

# ═══════════════════════════════════════════════════════════════════════════
# §4 cmdstan_valid guard (mirrors instantiate:::cmdstan_valid)
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_cmdstan_valid() requires dir + makefile + CMDSTAN_VERSION", {
  dir <- bvarnet_test_tempdir()
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)

  expect_false(.bvarnet_cmdstan_valid(""))
  expect_false(.bvarnet_cmdstan_valid(file.path(dir, "nonexistent")))
  expect_false(.bvarnet_cmdstan_valid(dir))  # dir exists but no makefile

  writeLines("SOME_OTHER_VAR := 1", file.path(dir, "makefile"))
  expect_false(.bvarnet_cmdstan_valid(dir))  # makefile without CMDSTAN_VERSION

  writeLines("CMDSTAN_VERSION := 2.39.0", file.path(dir, "makefile"))
  expect_true(.bvarnet_cmdstan_valid(dir))
})

# ═══════════════════════════════════════════════════════════════════════════
# §5 loader-path helper — prepend not replace, unset-vs-""
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_libpath_var() returns the right var per platform, NA on Windows", {
  local_mocked_bindings(
    .bvarnet_libpath_var = function() {
      sys <- "Darwin"
      if (sys == "Darwin") "DYLD_FALLBACK_LIBRARY_PATH" else NA_character_
    },
    .package = "bvarnet"
  )
  expect_equal(.bvarnet_libpath_var(), "DYLD_FALLBACK_LIBRARY_PATH")
})

test_that(".bvarnet_with_libpath() PREPENDS to an existing value and restores it exactly", {
  local_mocked_bindings(.bvarnet_libpath_var = function() "BVARNET_TEST_LIBPATH", .package = "bvarnet")
  old <- Sys.getenv("BVARNET_TEST_LIBPATH", unset = NA)
  on.exit(if (is.na(old)) Sys.unsetenv("BVARNET_TEST_LIBPATH")
          else Sys.setenv(BVARNET_TEST_LIBPATH = old), add = TRUE)

  Sys.setenv(BVARNET_TEST_LIBPATH = "/existing/path")
  seen <- .bvarnet_with_libpath("/cache/lib/tbb", Sys.getenv("BVARNET_TEST_LIBPATH"))
  expect_equal(seen, paste0("/cache/lib/tbb", .Platform$path.sep, "/existing/path"))
  expect_equal(Sys.getenv("BVARNET_TEST_LIBPATH"), "/existing/path")  # restored
})

test_that(".bvarnet_with_libpath() restores a truly-unset var to unset, not \"\"", {
  local_mocked_bindings(.bvarnet_libpath_var = function() "BVARNET_TEST_LIBPATH2", .package = "bvarnet")
  Sys.unsetenv("BVARNET_TEST_LIBPATH2")
  expect_true(is.na(Sys.getenv("BVARNET_TEST_LIBPATH2", unset = NA)))

  seen <- .bvarnet_with_libpath("/cache/lib/tbb", Sys.getenv("BVARNET_TEST_LIBPATH2"))
  expect_equal(seen, "/cache/lib/tbb")
  expect_true(is.na(Sys.getenv("BVARNET_TEST_LIBPATH2", unset = NA)))  # unset again, not ""
})

test_that(".bvarnet_with_libpath() is a no-op when the platform var is NA (Windows)", {
  local_mocked_bindings(.bvarnet_libpath_var = function() NA_character_, .package = "bvarnet")
  expect_equal(.bvarnet_with_libpath("/cache/lib/tbb", 41 + 1), 42)
})

# ═══════════════════════════════════════════════════════════════════════════
# §6 cache key / cache dir derivation
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_cache_key() combines version and an 8-char hash prefix", {
  local_mocked_bindings(
    .bvarnet_source_hash = function(...) paste0(rep("ab", 32), collapse = ""),
    .package = "bvarnet"
  )
  key <- .bvarnet_cache_key()
  expect_equal(key, paste0(as.character(utils::packageVersion("bvarnet")), "-abababab"))
})

# ═══════════════════════════════════════════════════════════════════════════
# §7 cache completeness / cache_model — verify-before-use
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_cache_complete() is FALSE when the marker is absent", {
  dir <- bvarnet_test_tempdir()
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  expect_false(.bvarnet_cache_complete(dir))
})

test_that(".bvarnet_cache_complete() verifies checksums and rejects a tampered exe", {
  dir <- bvarnet_test_tempdir()
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)

  exe <- file.path(dir, paste0("model_gaussian", .bvarnet_exe_ext()))
  writeBin(as.raw(1:10), exe)
  sha <- unname(digest::digest(exe, algo = "sha256", file = TRUE))
  jsonlite::write_json(list(sha256 = list(model_gaussian = sha)),
                       file.path(dir, "_complete.json"), auto_unbox = TRUE)
  expect_true(.bvarnet_cache_complete(dir))

  # Tamper with the exe after the marker was written -> must be treated as absent
  writeBin(as.raw(11:20), exe)
  expect_false(.bvarnet_cache_complete(dir))
})

test_that(".bvarnet_cache_model() returns NULL for an incomplete/unverified cache", {
  dir <- bvarnet_test_tempdir()
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  local_mocked_bindings(.bvarnet_models_cache_dir = function() dir, .package = "bvarnet")
  expect_null(.bvarnet_cache_model("model_gaussian"))
})

test_that(".bvarnet_cache_model() returns exe+stan paths for a verified cache", {
  dir <- bvarnet_test_tempdir()
  on.exit(unlink(dir, recursive = TRUE, force = TRUE), add = TRUE)
  local_mocked_bindings(.bvarnet_models_cache_dir = function() dir, .package = "bvarnet")

  ext <- .bvarnet_exe_ext()
  exe  <- file.path(dir, paste0("model_gaussian", ext))
  stan <- file.path(dir, "model_gaussian.stan")
  writeBin(as.raw(1:10), exe)
  writeLines("// stub", stan)
  sha <- unname(digest::digest(exe, algo = "sha256", file = TRUE))
  jsonlite::write_json(list(sha256 = list(model_gaussian = sha)),
                       file.path(dir, "_complete.json"), auto_unbox = TRUE)

  res <- .bvarnet_cache_model("model_gaussian")
  expect_equal(res$exe, exe)
  expect_equal(res$stan, stan)
  expect_equal(res$dir, dir)
})

# ═══════════════════════════════════════════════════════════════════════════
# §8 manifest fetch — three distinct failure modes
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_fetch_manifest() returns NULL on 404 / download failure", {
  local_mocked_bindings(
    download.file = function(...) stop("HTTP error 404"),
    .package = "utils"
  )
  expect_null(.bvarnet_fetch_manifest(version = "9.9.9"))
})

test_that(".bvarnet_fetch_manifest() errors clearly on malformed JSON", {
  local_mocked_bindings(
    download.file = function(url, destfile, ...) { writeLines("not json {{{", destfile); 0L },
    .package = "utils"
  )
  expect_error(.bvarnet_fetch_manifest(version = "1.0.1"), "could not be parsed as JSON")
})

test_that(".bvarnet_fetch_manifest() parses a well-formed manifest", {
  manifest <- list(
    version = "1.0.1", source_hash = "sha256:abc", cmdstan_version = "2.39.0",
    revoked = FALSE, revoked_reason = NULL,
    platforms = list(`macos-arm64` = list(asset = "bvarnet-models-1.0.1-macos-arm64.tar.gz",
                                          sha256 = "deadbeef"))
  )
  local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      jsonlite::write_json(manifest, destfile, auto_unbox = TRUE, null = "null")
      0L
    },
    .package = "utils"
  )
  m <- .bvarnet_fetch_manifest(version = "1.0.1")
  expect_equal(m$version, "1.0.1")
  expect_equal(m$source_hash, "sha256:abc")
  expect_false(isTRUE(m$revoked))
  expect_equal(m$platforms[["macos-arm64"]]$asset, "bvarnet-models-1.0.1-macos-arm64.tar.gz")
})

# ═══════════════════════════════════════════════════════════════════════════
# §9 download path — source_hash mismatch, revoked, missing platform, consent
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_download_models() aborts on a source_hash mismatch, never runs the binary", {
  local_mocked_bindings(
    .bvarnet_ensure_cmdstanr = function() invisible(TRUE),
    .bvarnet_fetch_manifest = function(...) list(
      version = "1.0.1", revoked = FALSE,
      source_hash = "sha256:deadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeefdeadbeef",
      platforms = list()
    ),
    .bvarnet_source_hash = function(...) "0000000000000000000000000000000000000000000000000000000000000000",
    .package = "bvarnet"
  )
  expect_error(.bvarnet_download_models(ask = FALSE), "source hash mismatch")
})

test_that(".bvarnet_download_models() falls back to compile-from-source when revoked", {
  local_mocked_bindings(
    .bvarnet_ensure_cmdstanr = function() invisible(TRUE),
    .bvarnet_fetch_manifest = function(...) list(
      version = "1.0.1", revoked = TRUE, revoked_reason = "bad TBB build"
    ),
    .package = "bvarnet"
  )
  expect_message(res <- .bvarnet_download_models(ask = FALSE), "revoked")
  expect_false(res)
})

test_that(".bvarnet_download_models() falls back when no manifest exists for this version", {
  local_mocked_bindings(
    .bvarnet_ensure_cmdstanr = function() invisible(TRUE),
    .bvarnet_fetch_manifest = function(...) NULL,
    .package = "bvarnet"
  )
  expect_message(res <- .bvarnet_download_models(ask = FALSE), "no precompiled-model release")
  expect_false(res)
})

test_that(".bvarnet_download_models() falls back when the platform key is absent from platforms", {
  h <- .bvarnet_source_hash()
  local_mocked_bindings(
    .bvarnet_ensure_cmdstanr = function() invisible(TRUE),
    .bvarnet_fetch_manifest = function(...) list(
      version = "1.0.1", revoked = FALSE, source_hash = paste0("sha256:", h),
      platforms = list()
    ),
    .bvarnet_platform_key = function(...) "made-up-platform",
    .package = "bvarnet"
  )
  expect_message(res <- .bvarnet_download_models(ask = FALSE), "no precompiled binary published")
  expect_false(res)
})

test_that("consent gating: non-interactive with no pre-authorisation never downloads", {
  old_opt <- getOption("bvarnet.allow_download")
  old_env <- Sys.getenv("BVARNET_ALLOW_DOWNLOAD", unset = NA)
  on.exit({
    options(bvarnet.allow_download = old_opt)
    if (is.na(old_env)) Sys.unsetenv("BVARNET_ALLOW_DOWNLOAD")
    else Sys.setenv(BVARNET_ALLOW_DOWNLOAD = old_env)
  }, add = TRUE)

  options(bvarnet.allow_download = NULL)
  Sys.unsetenv("BVARNET_ALLOW_DOWNLOAD")
  expect_false(.bvarnet_consent_download("asset.tar.gz", ask = FALSE))

  options(bvarnet.allow_download = TRUE)
  expect_true(.bvarnet_consent_download("asset.tar.gz", ask = FALSE))

  options(bvarnet.allow_download = NULL)
  Sys.setenv(BVARNET_ALLOW_DOWNLOAD = "true")
  expect_true(.bvarnet_consent_download("asset.tar.gz", ask = FALSE))
})

# ═══════════════════════════════════════════════════════════════════════════
# §10 models_state() — ok / stale / absent
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_models_state() is \"ok\" when install-tree exes are present", {
  local_mocked_bindings(.bvarnet_missing_models = function() character(0), .package = "bvarnet")
  expect_equal(.bvarnet_models_state(), "ok")
})

test_that(".bvarnet_models_state() is \"absent\" when nothing is cached anywhere", {
  root <- bvarnet_test_tempdir()  # exists but empty
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  local_mocked_bindings(
    .bvarnet_missing_models = function() c("model_gaussian"),
    .bvarnet_cache_root = function() root,
    .bvarnet_cache_key = function() "1.0.1-abcd1234",
    .package = "bvarnet"
  )
  expect_equal(.bvarnet_models_state(), "absent")
})

test_that(".bvarnet_models_state() is \"ok\" when the current key's cache is complete", {
  root <- bvarnet_test_tempdir()
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  cur  <- file.path(root, "1.0.1-abcd1234")
  dir.create(cur)
  local_mocked_bindings(
    .bvarnet_missing_models = function() c("model_gaussian"),
    .bvarnet_cache_root = function() root,
    .bvarnet_cache_key = function() "1.0.1-abcd1234",
    .bvarnet_cache_complete = function(dir) identical(dir, cur),
    .package = "bvarnet"
  )
  expect_equal(.bvarnet_models_state(), "ok")
})

test_that(".bvarnet_models_state() is \"stale\" when only a different-key cache dir exists", {
  root <- bvarnet_test_tempdir()
  on.exit(unlink(root, recursive = TRUE, force = TRUE), add = TRUE)
  dir.create(file.path(root, "1.0.0-deadbeef"))  # a prior version's cache
  local_mocked_bindings(
    .bvarnet_missing_models = function() c("model_gaussian"),
    .bvarnet_cache_root = function() root,
    .bvarnet_cache_key = function() "1.0.1-abcd1234",
    .bvarnet_cache_complete = function(dir) FALSE,
    .package = "bvarnet"
  )
  expect_equal(.bvarnet_models_state(), "stale")
})

test_that(".onAttach() never errors, even with a malformed cache root", {
  local_mocked_bindings(
    .bvarnet_models_state = function() stop("simulated failure reading cache"),
    .package = "bvarnet"
  )
  expect_silent(bvarnet:::.onAttach(NULL, "bvarnet"))
})

# ═══════════════════════════════════════════════════════════════════════════
# §11 resolver — errors when neither install-tree nor cache has the model
# ═══════════════════════════════════════════════════════════════════════════

test_that(".bvarnet_stan_model() raises the missing-models error when nothing is available", {
  local_mocked_bindings(
    .bvarnet_cache_model = function(model_name) NULL,
    .bvarnet_models_state = function() "absent",
    .package = "bvarnet"
  )
  expect_error(.bvarnet_stan_model("not_a_real_model"), "Stan models are not set up")
})

test_that(".bvarnet_model_for_family() maps families to bundled model names", {
  expect_equal(.bvarnet_model_for_family("bernoulli"), "model_binary")
  expect_equal(.bvarnet_model_for_family("ordinal"), "model_ordinal")
  expect_equal(.bvarnet_model_for_family("gaussian"), "model_gaussian")
  expect_error(.bvarnet_model_for_family("nope"), "Unknown family")
})
