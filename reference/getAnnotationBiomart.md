# Get Gene Annotation from Ensembl BioMart

Downloads and processes gene annotation data from Ensembl BioMart for
human, mouse, or rat species. Includes caching and offline fallback
functionality.

## Usage

``` r
getAnnotationBiomart(
  species = "human",
  annotation_dir = .get_data_path(),
  version = "current",
  chromosome_scope = c("all", "primary", "standard"),
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

- version:

  Character string specifying the Ensembl version to use. Defaults to
  "current" for the latest Ensembl version.

- chromosome_scope:

  Character string controlling which sequences are kept:

  - `"all"` (default): every gene Ensembl reports, including those on
    alternate haplotypes and patch scaffolds.

  - `"primary"`: the primary assembly only, i.e. assembled molecules
    plus unlocalized and unplaced scaffolds.

  - `"standard"`: assembled molecules only (1-22, X, Y, MT).

  Filtering is applied on the way out, so a single cache serves every
  scope.

- force_update:

  Logical; if `TRUE`, forces re-download and re-parsing.

- offline:

  Logical; if `TRUE`, forces offline mode and uses cached/bundled data.

## Value

A dataframe with the following columns:

- ensembl_id:

  Ensembl gene identifier

- gene_symbol:

  Official gene symbol

- gene_biotype:

  Biotype classification (protein_coding, lncRNA, etc.)

- chromosome:

  Chromosome location (see Details)

- gene_start:

  Genomic start position

- gene_end:

  Genomic end position

- strand:

  Strand (+ or -)

- gene_description:

  Gene description

- gene_synonyms:

  Pipe-separated alternative gene names

## Details

Assembled molecules are reported as bare names ("14", "X", "MT"). Every
other sequence is labelled `chr<molecule>(<GenBank accession>)` using
the NCBI assembly report for the build Ensembl is serving, so the
unlocalized scaffold carrying `MAFIP` is reported as
`chr14(GL000194.1)`. Scaffolds with no assigned molecule are labelled
`chrUn(<accession>)`.

Multiple synonyms are combined and separated by pipes (\|). Results are
cached locally with metadata about the Ensembl version and genome build
used.

## Examples

``` r
if (FALSE) { # \dontrun{
  human_genes <- getAnnotationBiomart("human")
  mouse_genes <- getAnnotationBiomart("mouse")

  # Drop alternate haplotypes and patches, keep scaffolds like chr14(GL000194.1)
  primary_only <- getAnnotationBiomart("human", chromosome_scope = "primary")
} # }
```
