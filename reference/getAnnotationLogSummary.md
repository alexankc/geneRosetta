# Summarize and Export Annotation Build Logs

Provides a human-readable summary of the gene annotation build process
for a specific species. It can print a detailed report to the console,
export the data as a formatted text file, generate manuscript-ready
Methods text, or copy the raw JSON log.

## Usage

``` r
getAnnotationLogSummary(
  species = c("human", "mouse", "rat"),
  ensembl_annot = FALSE,
  export_format = c("none", "txt", "methods", "json"),
  export_path = NULL,
  annotation_dir = .get_data_path()
)
```

## Arguments

- species:

  Character. One of `"human"`, `"mouse"`, or `"rat"`. Specifies which
  log file to retrieve. Defaults to `"human"`.

- ensembl_annot:

  Logical. If `TRUE`, summarises the Ensembl BioMart source log for the
  species instead of the merged annotation log. Defaults to `FALSE`.

- export_format:

  Character. The format to export the log:

  - `"none"`: Just print the summary to the R console (default).

  - `"txt"`: Create a formatted, human-readable text file.

  - `"methods"`: Create a prose summary suitable for a scientific
    Methods section.

  - `"json"`: Copy the raw JSON log file to a local destination.

- export_path:

  Character (optional). The file path where the log should be saved. If
  the directory does not exist, it will be created recursively.

- annotation_dir:

  Character string specifying the directory holding the annotation logs.
  Defaults to the persistent user data directory, which is where
  [`buildAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildAnnotation.md)
  writes unless it was given a custom directory.

## Value

The log data as a nested list, returned invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
# Print human log to console
getAnnotationLogSummary(species = "human")

# Export mouse log to a specific folder
getAnnotationLogSummary(
  species = "mouse",
  export_format = "methods",
  export_path = "logs/mouse_annotation_methods.txt"
)
} # }
```
