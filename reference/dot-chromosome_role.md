# Recover the Sequence Role of a Formatted Chromosome Label

Recover the Sequence Role of a Formatted Chromosome Label

## Usage

``` r
.chromosome_role(labels, lookup)
```

## Arguments

- labels:

  Character vector of chromosome labels as produced by
  [`.format_chromosome()`](https://alexankc.github.io/geneRosetta/reference/dot-format_chromosome.md).

- lookup:

  A list from
  [`.chromosome_lookup()`](https://alexankc.github.io/geneRosetta/reference/dot-chromosome_lookup.md),
  or NULL. With NULL only assembled molecules can be recognised;
  everything else is `"unknown"`.

## Value

A character vector of sequence roles; `"unknown"` where the label is
neither an assembled molecule nor a scaffold present in the report.
