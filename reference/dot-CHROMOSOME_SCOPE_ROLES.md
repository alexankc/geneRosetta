# Sequence Roles Retained by Each Chromosome Scope

`"standard"` keeps only assembled molecules (the historic behaviour),
`"primary"` additionally keeps the unlocalized and unplaced scaffolds
that make up the rest of the primary assembly, and `"all"` keeps
everything, including alternate haplotypes (alt-scaffolds) and patch
releases.

## Usage

``` r
.CHROMOSOME_SCOPE_ROLES
```

## Format

An object of class `list` of length 3.
