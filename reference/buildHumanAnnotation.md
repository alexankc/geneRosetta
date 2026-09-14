# Build Merged Human Gene Annotation

Merges gene annotation data from multiple sources (BioMart/Ensembl,
NCBI, and HGNC) into a unified human gene annotation dataset. Implements
multi-tier matching strategy and intelligent symbol prioritization.

## Usage

``` r
buildHumanAnnotation(
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

  Official gene symbol (priority: HGNC \> NCBI \> BioMart)

- gene_biotype:

  Gene biotype classification

- chromosome:

  Chromosome location

- gene_start:

  Genomic start position

- gene_end:

  Genomic end position

- strand:

  Strand orientation (+ or -)

- metadata_source:

  Metadata source: Ensembl or NCBI

- gene_synonyms:

  Pipe-separated merged gene synonyms from all sources

- ncbi_id:

  NCBI Gene ID

- hgnc_id:

  HGNC identifier

- gene_description:

  Gene description

## Details

The merging strategy uses the following tiers:

1.  BioMart -\> HGNC (via Ensembl ID)

2.  BioMart -\> NCBI (via Ensembl ID) for unmatched BioMart records

3.  Attempt to rescue NCBI IDs via symbol matching

4.  Attempt to rescue HGNC IDs via symbol matching

Gene symbols are prioritized as: HGNC (official) \> NCBI \> BioMart.
Gene synonyms are deduplicated, sorted, and merged with pipes (\|) as
separators. Results are cached locally as TSV with accompanying JSON
metadata log.

## Examples

``` r
if (FALSE) { # \dontrun{
human_annotation <- buildHumanAnnotation()
human_annotation <- buildHumanAnnotation(force_update = TRUE)
} # }
```
