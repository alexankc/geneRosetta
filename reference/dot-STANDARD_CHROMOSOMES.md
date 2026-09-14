# Standard (Assembled-Molecule) Chromosome Names

The chromosome tokens that Ensembl reports for assembled molecules
across the species supported by this package. Used as a fast path by
[`.format_chromosome()`](https://alexankc.github.io/geneRosetta/reference/dot-format_chromosome.md):
anything in this set is passed through unchanged, anything else is
looked up in an NCBI assembly report.

## Usage

``` r
.STANDARD_CHROMOSOMES
```

## Format

An object of class `character` of length 25.
