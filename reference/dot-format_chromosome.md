# Label Non-Standard Chromosomes

Rewrites scaffold names to `chr<molecule>(<GenBank accession>)`, e.g.
`GL000194.1` becomes `chr14(GL000194.1)`. Scaffolds with no assigned
molecule become `chrUn(<accession>)`. Assembled molecules (1-22, X, Y,
MT) and names absent from the assembly report are returned unchanged.

## Usage

``` r
.format_chromosome(x, lookup)
```

## Arguments

- x:

  Character vector of chromosome names as reported by the source.

- lookup:

  A list from
  [`.chromosome_lookup()`](https://alexankc.github.io/geneRosetta/reference/dot-chromosome_lookup.md),
  or NULL to pass `x` through unchanged.

## Value

A character vector the same length as `x`.
