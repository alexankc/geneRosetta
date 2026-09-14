# Coerce a Cached Assembly Report Back to Character

[`.manage_cache()`](https://alexankc.github.io/geneRosetta/reference/dot-manage_cache.md)
reads TSVs with
[`read.delim()`](https://rdrr.io/r/utils/read.table.html), which
type-guesses columns; sequence names such as "1" come back as integers.
Every field of an assembly report is an identifier, so force them all to
character.

## Usage

``` r
.as_assembly_report(df)
```

## Arguments

- df:

  A dataframe read from cache or bundled data.

## Value

A tibble with all columns as character.
