# Contributing to bvarnet

To be added...


## Release checklist (precompiled model binaries)

`bvarnet_setup_models(method = "download")` resolves a GitHub Release asset from the *installed* package version (`releases/download/v<version>/...`). This makes the release process itself part of the package's correctness, not just packaging hygiene. From this, two rules follow:

### 1. Never delete or mutate a published release or its assets

Once a release's assets are published, they must stay exactly as they are, forever (with the one narrow exception below). Deleting or re-uploading an asset breaks
`bvarnet_setup_models(method = "download")` for every user still pinned to that version. An old script, a reproducible analysis, or a citation of the package would silently lose its precompiled binaries with no warning! This matters more than usual for an academic package that others may cite and expect to keep working.

- **To fix a bad build**, publish a **new** patch version with corrected assets. Do not overwrite the old release.
- **Deleting a release only ever *downgrades* affected users to compile-from-source** (the   `.stan` sources always ship in the package, so that path still works) but it is a silent regression for the toolchain-free users this feature exists for, so don't do it.
- **The one exception:** a shipped-but-broken binary (e.g. a bad build, a TBB regression) can be disabled by setting `"revoked": true` (and a short `"revoked_reason"`) in that release's `manifest.json` and re-uploading *only the manifest* -- never the binaries. Clients that see `revoked: true` skip the download and fall back to compile-from-source, surfacing the reason.

### 2. Publish a release with assets for every installable version and don't micro-release

The client resolves its download by version string, not by whether the Stan models actually changed. So any version a user can install (every CRAN) needs a matching release + `manifest.json`, or `bvarnet_setup_models(method = "download")` 404s and silently falls back to compile-from-source for that version's entire user base.

| Change | Publish a release + assets? | New model binaries? |
| --- | --- | --- |
| Update that edits a `.stan` file | **Yes** | Yes -- new `source_hash` |
| Patch (R code only) shipped to users | **Yes** | No -- same `source_hash`, rebuilt copies |
| Intermediate commit / `.9000` dev version, never shipped | No (404 -> compile-from-source) | n/a |

- **Don't micro-release.** Batch small fixes into fewer, meaningful versions. Every released version costs one full asset set (~150-220 MB across all platforms) whether or not the Stan models changed, so frequent tiny releases are the main bloat risk. This also follows CRAN's own guidance against submitting more than roughly every 1-2 months.
- **Dev/pre-release tags get no assets.** Only versions users actually install need them.
- **The release tag MUST equal `paste0("v", packageVersion())` exactly** -- e.g. `v1.0.2`, not `1.0.2` or `v1.0.2.9000`. The client builds the download URL from this string; any mismatch 404s the whole platform silently. Double-check the tag against DESCRIPTION`'s `Version` field before publishing.
