# Save Annotation Data and Metadata

Writes downloaded annotation data to a TSV file and saves associated
metadata (including download timestamp and source information) to a JSON
log file.

## Usage

``` r
.save_annotation(df, mapping_file, log_file, metadata)
```

## Arguments

- df:

  A dataframe containing the gene annotation data to save.

- mapping_file:

  Character string path where the TSV file will be written.

- log_file:

  Character string path where the JSON metadata log will be written.

- metadata:

  A named list of metadata to include in the JSON log (e.g., version
  info, URLs, package versions).

## Value

NULL (invisibly). Called for its side effect of writing files to disk.
