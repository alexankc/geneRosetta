#' @keywords internal
"_PACKAGE"

#' Imports used across the package
#'
#' \code{\%||\%} is taken from rlang rather than base R: base only gained it in
#' R 4.4.0, and the package supports R (>= 4.1).
#'
#' @importFrom rlang %||% :=
#' @importFrom stats na.omit setNames
#' @importFrom tools R_user_dir
#' @importFrom utils download.file packageVersion read.delim write.table
#' @name geneRosetta-imports
#' @keywords internal
NULL

# Column names referenced through tidy evaluation. Declaring them keeps
# 'R CMD check' from reporting them as undefined global variables.
utils::globalVariables(c(
  # magrittr's pipe placeholder, used in getAnnotationMGImouse()
  ".",
  # BioMart / Ensembl
  "ensembl_gene_id", "external_gene_name", "external_synonym", "chromosome_name",
  "start_position", "end_position", "description",
  # NCBI gene_info
  "dbXrefs", "GeneID", "Symbol", "Synonyms",
  # HGNC
  "alias_symbol", "prev_symbol", "symbol", "entrez_id",
  # MGI
  "MGI_accession_id", "marker_symbol", "marker_name", "Feature_Type",
  "Ensembl_gene_chromosome", "Ensembl_gene_start", "Ensembl_gene_end",
  "Ensembl_gene_strand", "Marker_Synonyms__pipe_separated_", "Entrez_gene_id",
  "Ensembl_gene_id",
  # RGD
  "GENE_RGD_ID", "SYMBOL", "OLD_SYMBOL", "GENE_TYPE", "ENSEMBL_ID",
  "NCBI_GENE_ID", "NAME",
  # shared annotation schema
  "ensembl_id", "gene_symbol", "gene_synonyms", "gene_biotype", "chromosome",
  "gene_start", "gene_end", "strand", "gene_description",
  "ncbi_id", "hgnc_id", "mgi_id", "rgd_id",
  # merge intermediates
  "symbol_bm", "symbol_ncbi", "symbol_hgnc",
  "syn_bm", "syn_ncbi", "syn_hgnc", "syn_mgi", "syn_rgd",
  "syn_ncbi.x", "syn_ncbi.y", "syn_hgnc.x", "syn_hgnc.y",
  "gene_symbol.x", "gene_symbol.y", "temp_synonyms",
  # mapping helpers
  "query_gene_symbol", "query_gene_symbol_clean", "row_id_internal", "temp_id"
))
