# Resolve the URL of an NCBI Assembly Report

Ensembl reports its genome build as the same string NCBI uses to suffix
the assembly directory ("GRCh38.p14", "GRCm39", "GRCr8"), so the report
for the exact build Ensembl is serving can be located by listing the
species' accession directory. If the listing is unreachable or the build
is absent, falls back to the pinned assembly and warns.

## Usage

``` r
.assembly_report_url(species, genome_build = NULL)
```

## Arguments

- species:

  Character string. One of "human", "mouse", or "rat".

- genome_build:

  Character string of the Ensembl genome build, or NULL/NA to go
  straight to the pinned assembly.

## Value

Character string URL of the `*_assembly_report.txt` file.
