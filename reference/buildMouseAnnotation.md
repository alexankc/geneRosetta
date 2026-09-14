# Build Merged Mouse Gene Annotation

Merges gene annotation data from multiple sources (BioMart/Ensembl, MGI,
and NCBI) into a unified mouse gene annotation dataset. Implements
multi-tier matching with preference for MGI nomenclature as the official
mouse gene naming authority.

## Usage

``` r
buildMouseAnnotation(
  annotation_dir = .get_data_path(),
  force_update = FALSE,
  offline = getOption("geneRosetta.offline", FALSE)
)
```

## Arguments

- annotation_dir:

  Character string specifying the directory. Defaults to the persistent
  user data directory.

- force_update:

  Logical. If TRUE, forces re-download and rebuild of annotations even
  if a valid cache exists. Defaults to FALSE.

## Value

A dataframe with the following columns:

- ensembl_id:

  Ensembl gene identifier (primary key)

- gene_symbol:

  Official gene symbol (prioritizes MGI nomenclature)

- gene_biotype:

  Gene biotype classification from BioMart

- chromosome:

  Chromosome location

- gene_start:

  Genomic start position from BioMart/Ensembl

- gene_end:

  Genomic end position from BioMart/Ensembl

- strand:

  Strand orientation (+ or -)

- metadata_source:

  Metadata source: Ensembl or NCBI

- gene_synonyms:

  Pipe-separated merged gene synonyms from all sources

- ncbi_id:

  NCBI Gene ID

- mgi_id:

  MGI accession ID

- gene_description:

  Gene description

## Details

The merging strategy uses multiple tiers:

1.  BioMart -\> MGI (perfect match via Ensembl ID AND Symbol)

2.  BioMart -\> MGI (Ensembl ID match with updated symbol from MGI)

3.  BioMart -\> NCBI (rescue for MGI orphans via Ensembl ID)

4.  Rescue missing MGI IDs via NCBI ID matching

5.  Rescue missing MGI IDs via symbol matching

MGI provides the official mouse nomenclature and is preferred for gene
symbols when available. Genomic coordinates are sourced from
BioMart/Ensembl for consistency. Results are cached locally as TSV with
accompanying JSON metadata log.

## Examples

``` r
if (FALSE) { # \dontrun{
mouse_annotation <- buildMouseAnnotation()
mouse_annotation <- buildMouseAnnotation(force_update = TRUE)
} # }
```
