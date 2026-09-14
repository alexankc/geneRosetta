# Get Gene Annotation from NCBI Gene Database

Downloads and processes gene annotation data from NCBI Gene database for
human, mouse, or rat species. Includes caching and offline fallback.

## Usage

``` r
getAnnotationNCBI(
  species = "human",
  annotation_dir = .get_data_path(),
  force_update = FALSE,
  offline = getOption("geneRosetta.offline", FALSE)
)
```

## Arguments

- species:

  Character string specifying the species. One of "human" (default),
  "mouse", or "rat".

- annotation_dir:

  Character string specifying the directory where annotation data will
  be cached. Defaults to the persistent user data directory.

- force_update:

  Logical; if `TRUE`, forces re-download and re-parsing.

- offline:

  Logical; if `TRUE`, forces offline mode and uses cached/bundled data.

## Value

A dataframe with the following columns:

- ncbi_id:

  NCBI Gene ID (GeneID)

- gene_symbol:

  Official gene symbol

- gene_synonyms:

  Pipe-separated alternative gene names

- ensembl_id:

  Corresponding Ensembl gene identifier (if available)

- gene_biotype:

  NCBI gene type from `type_of_gene`

- gene_description:

  NCBI gene description

## Details

Ensembl IDs are extracted from the dbXrefs field where available.
Multiple synonyms are combined and separated by pipes (\|). Data is
cached with modification timestamps for incremental updates.

## Examples

``` r
if (FALSE) { # \dontrun{
  human_genes <- getAnnotationNCBI("human")
  mouse_genes <- getAnnotationNCBI("mouse")
} # }
```
