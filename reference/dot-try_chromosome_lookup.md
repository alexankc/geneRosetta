# Build a Chromosome Lookup Without Failing

Sources whose chromosome names are already plain (MGI, RGD) should not
lose an entire download because the assembly report is unreachable.
Returns NULL on any failure, which
[`.format_chromosome()`](https://alexankc.github.io/geneRosetta/reference/dot-format_chromosome.md)
treats as "leave names untouched".

## Usage

``` r
.try_chromosome_lookup(
  species,
  annotation_dir = .get_data_path(),
  genome_build = NULL,
  offline = getOption("geneRosetta.offline", FALSE)
)
```

## Arguments

- species:

  Character string. One of "human", "mouse", or "rat".

- annotation_dir:

  Character string of the cache directory.

- genome_build:

  Character string of the Ensembl genome build (e.g. "GRCh38.p14"). Used
  to pick the matching NCBI assembly and to invalidate a cache built
  against a different build.

- offline:

  Logical.

## Value

A list from
[`.chromosome_lookup()`](https://alexankc.github.io/geneRosetta/reference/dot-chromosome_lookup.md),
or NULL.
