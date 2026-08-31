# Remove bvarnet's cached Stan model binaries

Remove bvarnet's cached Stan model binaries

## Usage

``` r
bvarnet_clear_model_cache(all_versions = FALSE, quiet = FALSE)
```

## Arguments

- all_versions:

  Logical. Remove the entire cache (all versions), not just the entry
  for the currently installed version.

- quiet:

  Logical. Suppress the confirmation message.

## Value

Invisibly, the path removed.

## See also

[`bvarnet_setup_models()`](https://flo1met.github.io/bvarnet/reference/bvarnet_setup_models.md)
to set models up again afterward,
[`bvarnet_model_cache_dir()`](https://flo1met.github.io/bvarnet/reference/bvarnet_model_cache_dir.md)
to inspect the cache path.
