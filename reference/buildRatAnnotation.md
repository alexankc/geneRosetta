# Build Merged Rat Gene Annotation

Merges gene annotation data from multiple sources (BioMart/Ensembl, RGD,
and NCBI) into a unified rat gene annotation dataset. Implements
multi-tier matching with RGD as the preferred source for rat
nomenclature and identity validation.

## Usage

``` r
buildRatAnnotation(
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

  Official gene symbol (prioritizes RGD nomenclature)

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

  NCBI Gene ID (validated against NCBI reference)

- rgd_id:

  RGD identifier

- gene_description:

  Gene description

## Details

The merging strategy uses multiple tiers:

1.  BioMart -\> RGD (perfect match via Ensembl ID AND Symbol)

2.  BioMart -\> RGD (Ensembl ID match with updated symbol from RGD)

3.  BioMart -\> NCBI (rescue for RGD orphans via Ensembl ID)

4.  Validate multiple NCBI IDs against the NCBI reference database

RGD provides the official rat nomenclature and is preferred for gene
symbols when available. Genomic coordinates are sourced from
BioMart/Ensembl for consistency. NCBI IDs are validated and cleaned when
multiple values are present (separated by semicolons in the source
data). Results are cached locally as TSV with accompanying JSON metadata
log.

## Examples

``` r
if (FALSE) { # \dontrun{
rat_annotation <- buildRatAnnotation()
rat_annotation <- buildRatAnnotation(force_update = TRUE)
} # }
```
