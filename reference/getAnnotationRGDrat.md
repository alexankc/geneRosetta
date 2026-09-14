# Get Rat Gene Annotation from RGD Database

Downloads and processes rat gene annotation data from the Rat Genome
Database (RGD). Handles assembly-specific data selection and cleaning.

## Usage

``` r
getAnnotationRGDrat(
  annotation_dir = .get_data_path(),
  genome_build = "GRCr8",
  force_update = FALSE,
  offline = getOption("geneRosetta.offline", FALSE)
)
```

## Arguments

- annotation_dir:

  Character string specifying the directory where annotation data will
  be cached. Defaults to the persistent user data directory.

- genome_build:

  Character string specifying the rat genome build version. Defaults to
  "GRCr8". Must match column names in the RGD data file.

- force_update:

  Logical; if `TRUE`, forces re-download and re-parsing.

- offline:

  Logical; if `TRUE`, forces offline mode and uses cached/bundled data.

## Value

A dataframe with the following columns for rat genes:

- rgd_id:

  RGD gene ID

- gene_symbol:

  Official gene symbol

- gene_description:

  Gene name/description

- gene_biotype:

  Gene type classification

- chromosome:

  Chromosome location

- gene_start:

  Genomic start position

- gene_end:

  Genomic end position

- strand:

  Strand (+ or -)

- gene_synonyms:

  Pipe-separated alternative gene names

- ensembl_id:

  Ensembl gene identifier

- ncbi_id:

  NCBI Gene ID (Entrez ID)

## Details

The function downloads the primary GENES_RAT.txt file from RGD. Column
names are dynamically adjusted based on the specified `genome_build`.
Semicolon-separated Ensembl IDs are expanded into separate rows.

## Examples

``` r
if (FALSE) { # \dontrun{
  rat_genes <- getAnnotationRGDrat(genome_build = "GRCr8")
} # }
```
