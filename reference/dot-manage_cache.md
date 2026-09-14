# Manage Caching and Offline Fallback

Checks if cached data is available and valid by comparing local metadata
against a remote timestamp. Implements a three-tier priority logic:

1.  Remote: If online and cache is stale, signals a re-download.

2.  Local Cache: If offline or cache is fresh, returns user-downloaded
    data.

3.  Bundled Data: If no local cache exists and offline, returns internal
    package data.

## Usage

``` r
.manage_cache(
  mapping_file,
  log_file,
  remote_stamp,
  sub_dir = NULL,
  force_update = FALSE,
  exact_version = FALSE,
  offline = getOption("geneRosetta.offline", FALSE),
  bundled_root = .bundled_extdata_root()
)
```

## Arguments

- mapping_file:

  Character string path to the cached mapping TSV file.

- log_file:

  Character string path to the JSON log file containing metadata.

- remote_stamp:

  Date/Character representing the remote data's version/date.

- sub_dir:

  Character string specifying the subdirectory within `inst/extdata`
  where bundled fallback data resides.

- force_update:

  Logical; if `TRUE`, bypasses the local cache and forces a re-download.

- exact_version:

  Logical; if `TRUE`, requires exact equality between cached and remote
  version.

- offline:

  Logical; if `TRUE`, forces offline mode and uses cached/bundled data.

- bundled_root:

  Path to root extdata directory.

## Value

A list with two elements:

- data:

  A dataframe (from local cache or bundled internal data) or NULL if a
  re-download is required.

- download:

  Logical; `TRUE` if the server is reachable and the remote data is
  newer than the cache, `FALSE` otherwise.
