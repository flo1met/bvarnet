#!/usr/bin/env Rscript
# check-windows-static-link.R
#
# Asserts that a compiled Stan model .exe carries NO dependency on the MinGW C++
# runtime DLLs -- i.e. that the build statically linked them. A regression here means a
# toolchain-free Windows user's download would fail to load with a missing-DLL
# error, so this is ship-blocking and belongs in CI, not a manual checklist.
#
# What it checks, per exe:
#   HARD FAIL  if any of the forbidden MinGW runtime DLLs is imported
#              (libstdc++-6, libgcc_s_seh-1, libgcc_s_dw2-1, libwinpthread-1).
#   NOTE       whether tbb.dll is imported. It is expected to be (TBB stays
#              dynamic and is bundled + co-located), but if a future build ever
#              static-links TBB too that is strictly better, so its ABSENCE is
#              only informational -- never a failure.
#
# It works off either `objdump -p` (from Rtools/mingw, usually already on PATH
# after r-lib/actions/setup-r) or `dumpbin /dependents` (from the preinstalled
# Visual Studio) -- whichever is available. Same script a human can run locally.
#
# Usage:
#   Rscript check-windows-static-link.R path/to/model_binary.exe [more.exe ...]
# Exit code 0 = all clean, 1 = a forbidden DLL was found (or no exe/tool found).

FORBIDDEN <- c(
  "libstdc++-6.dll",
  "libgcc_s_seh-1.dll",
  "libgcc_s_dw2-1.dll",   # 32-bit variant, belt-and-suspenders
  "libwinpthread-1.dll"
)

`%||%` <- function(a, b) if (length(a) == 0L || is.null(a)) b else a

# Return a lowercase, de-duplicated character vector of imported DLL names for
# `exe`, or NULL if no tool could be run. Parses whichever tool is present.
imported_dlls <- function(exe) {
  parse_objdump <- function(txt) {
    # objdump -p prints lines like:  \tDLL Name: tbb.dll
    m <- regmatches(txt, regexpr("(?i)DLL Name:\\s*\\S+", txt, perl = TRUE))
    sub("(?i)DLL Name:\\s*", "", m, perl = TRUE)
  }
  parse_dumpbin <- function(txt) {
    # dumpbin /dependents lists each dependency indented on its own line,
    # under "Image has the following dependencies:". Grab bare *.dll tokens.
    tok <- trimws(txt)
    tok[grepl("(?i)^\\S+\\.dll$", tok, perl = TRUE)]
  }

  run <- function(cmd, args) {
    out <- suppressWarnings(tryCatch(
      system2(cmd, args, stdout = TRUE, stderr = TRUE),
      error = function(e) NULL
    ))
    if (is.null(out) || !is.null(attr(out, "status")) && attr(out, "status") != 0L && !length(out))
      return(NULL)
    out
  }

  # Prefer objdump (reliably on PATH via Rtools); fall back to dumpbin.
  if (nzchar(Sys.which("objdump"))) {
    out <- run("objdump", c("-p", shQuote(exe)))
    if (!is.null(out)) return(unique(tolower(parse_objdump(out))))
  }
  if (nzchar(Sys.which("dumpbin"))) {
    out <- run("dumpbin", c("/dependents", shQuote(exe)))
    if (!is.null(out)) return(unique(tolower(parse_dumpbin(out))))
  }
  NULL
}

main <- function(args) {
  if (!length(args)) {
    message("usage: Rscript check-windows-static-link.R <model.exe> [more.exe ...]")
    quit(status = 1L)
  }
  forbidden <- tolower(FORBIDDEN)
  failed <- FALSE

  for (exe in args) {
    cat(sprintf("\n== %s ==\n", exe))
    if (!file.exists(exe)) {
      cat("  ERROR: file not found\n"); failed <- TRUE; next
    }
    dlls <- imported_dlls(exe)
    if (is.null(dlls)) {
      cat("  ERROR: neither `objdump` nor `dumpbin` is available on PATH\n")
      failed <- TRUE; next
    }
    cat("  imports:", paste(dlls, collapse = ", ") %||% "(none)", "\n")

    hits <- intersect(dlls, forbidden)
    if (length(hits)) {
      cat("  FAIL: MinGW runtime DLL(s) present (static-linking regressed):",
          paste(hits, collapse = ", "), "\n")
      failed <- TRUE
    } else {
      cat("  OK: no MinGW runtime DLLs\n")
    }
    if ("tbb.dll" %in% dlls) {
      cat("  note: links tbb.dll (expected; must be co-located beside the exe at staging)\n")
    } else {
      cat("  note: no tbb.dll import (TBB may be static-linked -- fine, just verify co-location isn't relied on)\n")
    }
  }

  if (failed) {
    cat("\nRESULT: FAIL -- see above\n"); quit(status = 1L)
  }
  cat("\nRESULT: PASS -- no forbidden MinGW DLLs in any exe\n"); quit(status = 0L)
}

main(commandArgs(trailingOnly = TRUE))
