# Get the NCBI Assembly Report for a Species

Downloads the assembly report matching the Ensembl genome build, caching
it alongside the annotation data with the same offline fallback
behaviour as the other sources. The report is the reference that turns
an Ensembl scaffold name into a human-readable chromosome label such as
`chr14(GL000194.1)`.

## Usage

``` r
.get_assembly_report(
  species,
  annotation_dir = .get_data_path(),
  genome_build = NULL,
  force_update = FALSE,
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

- force_update:

  Logical.

- offline:

  Logical.

## Value

A tibble of the assembly report.
