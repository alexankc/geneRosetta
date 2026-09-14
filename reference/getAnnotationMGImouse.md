# Get Mouse Gene Annotation from MGI Database

Downloads and processes mouse gene annotation data from the Mouse Genome
Informatics (MGI) database. Combines genetic markers with coordinate
data.

## Usage

``` r
getAnnotationMGImouse(
  annotation_dir = .get_data_path(),
  offline = getOption("geneRosetta.offline", FALSE),
  force_update = FALSE
)
```

## Arguments

- annotation_dir:

  Character string specifying the directory where annotation data will
  be cached. Defaults to the persistent user data directory.

- offline:

  Logical. If TRUE, uses cached/bundled package data instead of
  downloading. Defaults to FALSE.

- force_update:

  Logical; if `TRUE`, forces re-download and re-parsing.

## Value

A dataframe with the following columns for mouse genes:

- mgi_id:

  MGI marker accession ID

- gene_symbol:

  Official gene symbol

- gene_synonyms:

  Pipe-separated alternative gene names

- ensembl_id:

  Ensembl gene identifier

- ncbi_id:

  NCBI Gene ID (Entrez ID)

- chromosome:

  Chromosome location

- gene_start:

  Genomic start position (cM)

- gene_end:

  Genomic end position (cM)

- strand:

  Strand (+ or -)

- gene_biotype:

  Marker type (e.g., "Gene", "Pseudogene")

- gene_description:

  Marker name/description

## Details

The function downloads two files from MGI:

- MRK_List2.rpt: Marker catalog with symbols, names, biotypes, and
  synonyms

- MGI_Gene_Model_Coord.rpt: Genomic coordinates and Ensembl/NCBI
  cross-references

Data is combined and deduplicated. Only markers with type "Gene" or
"Pseudogene" are retained.

## Examples

``` r
if (FALSE) { # \dontrun{
  mouse_genes <- getAnnotationMGImouse()
} # }
```
