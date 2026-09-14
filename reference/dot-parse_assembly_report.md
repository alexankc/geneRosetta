# Parse an NCBI Assembly Report

Assembly reports are tab-separated with a leading '#' comment block
whose last line carries the column names. Line endings are CRLF in the
comment block, so carriage returns are stripped before parsing.

## Usage

``` r
.parse_assembly_report(path)
```

## Arguments

- path:

  Character string path to the downloaded report.

## Value

A tibble with the ten `.ASSEMBLY_REPORT_COLS`, all character.
