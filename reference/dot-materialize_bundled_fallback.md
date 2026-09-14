# Transactionally Materialize Bundled Fallback into Local Cache

Transactionally Materialize Bundled Fallback into Local Cache

## Usage

``` r
.materialize_bundled_fallback(
  mapping_file,
  log_file,
  sub_dir = NULL,
  bundled_root = .bundled_extdata_root()
)
```

## Arguments

- mapping_file:

  Destination mapping TSV path.

- log_file:

  Destination JSON log path.

- sub_dir:

  Subdirectory within extdata.

- bundled_root:

  Root extdata directory.

## Value

A list with `data` (tibble) and `success` (logical).
