# Get Human Gene Annotation from HGNC Database

Downloads and processes human gene annotation data from the HGNC (Human
Genome Nomenclature Committee) database. Includes caching and offline
fallback functionality.

## Usage

``` r
getAnnotationHGNChuman(
  annotation_dir = .get_data_path(),
  force_update = FALSE,
  offline = getOption("geneRosetta.offline", FALSE)
)
```

## Arguments

- annotation_dir:

  Character string specifying the directory where annotation data will
  be cached. Defaults to the persistent user data directory.

- force_update:

  Logical; if `TRUE`, forces re-download and re-parsing.

- offline:

  Logical; if `TRUE`, forces offline mode and uses cached/bundled data.

## Value

A dataframe with the following columns for human genes:

- hgnc_id:

  HGNC identifier

- gene_symbol:

  Official HGNC gene symbol

- gene_synonyms:

  Pipe-separated alternative gene names

- ncbi_id:

  NCBI Gene ID (Entrez ID)

- ensembl_id:

  Ensembl gene identifier

## Details

The function downloads the complete HGNC dataset in TSV format. Only
essential mapping columns are retained. Multiple synonyms and previous
symbols are merged and separated by pipes (\|).

## Examples

``` r
if (FALSE) { # \dontrun{
  hgnc_genes <- getAnnotationHGNChuman()
} # }
```
