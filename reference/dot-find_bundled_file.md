# Locate a Bundled File Within Extdata (preferring .gz over raw)

Locate a Bundled File Within Extdata (preferring .gz over raw)

## Usage

``` r
.find_bundled_file(sub_dir, file_name, bundled_root = .bundled_extdata_root())
```

## Arguments

- sub_dir:

  Subdirectory within extdata (e.g., "Ensembl_Genes").

- file_name:

  Base file name to find (e.g., "human_biomart_mapping.tsv").

- bundled_root:

  Root extdata directory. Defaults to
  [`.bundled_extdata_root()`](https://alexankc.github.io/geneRosetta/reference/dot-bundled_extdata_root.md).

## Value

A list with `path` and `is_compressed`, or NULL if not found.
