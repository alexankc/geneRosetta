# Restrict an Annotation Table to a Chromosome Scope

Restrict an Annotation Table to a Chromosome Scope

## Usage

``` r
.filter_chromosome_scope(df, scope, lookup)
```

## Arguments

- df:

  A dataframe with a `chromosome` column of formatted labels.

- scope:

  One of "all", "primary", or "standard".

- lookup:

  A list from
  [`.chromosome_lookup()`](https://alexankc.github.io/geneRosetta/reference/dot-chromosome_lookup.md),
  or NULL to skip filtering entirely.

## Value

`df`, filtered to the sequence roles the scope admits.

## Details

Without an assembly report there is no way to tell an unplaced scaffold
from an alt-haplotype, so a NULL `lookup` returns `df` untouched rather
than dropping rows it cannot classify. This mirrors
[`.format_chromosome()`](https://alexankc.github.io/geneRosetta/reference/dot-format_chromosome.md),
which leaves names as published when the report is unavailable.
