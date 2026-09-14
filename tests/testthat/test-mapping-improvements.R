# Tests for mapping robustness, row preservation, case/version normalization,
# and caching version/force_update improvements.

test_that("exact Ensembl version caching prevents loading newer cached versions", {
  dir <- withr::local_tempdir()
  m_file_116 <- file.path(dir, "human_biomart_mapping.tsv")
  l_file_116 <- file.path(dir, "human_biomart_mapping_log.json")
  
  df_116 <- data.frame(
    ensembl_id = "ENSG001",
    gene_symbol = "GENE1",
    stringsAsFactors = FALSE
  )
  .save_annotation(df_116, m_file_116, l_file_116, list(ensembl_version = "Ensembl_Genes_116"))
  
  # When version is "current", newer cache is accepted
  cache_curr <- .manage_cache(m_file_116, l_file_116, "Ensembl_Genes_114", exact_version = FALSE)
  expect_false(cache_curr$download)
  
  # When an explicit version is requested, exact match is required
  cache_exact <- .manage_cache(m_file_116, l_file_116, "Ensembl_Genes_114", exact_version = TRUE)
  expect_true(cache_exact$download)
})

test_that("force_update = TRUE signals re-download regardless of cache status", {
  dir <- withr::local_tempdir()
  m_file <- file.path(dir, "test_mapping.tsv")
  l_file <- file.path(dir, "test_log.json")
  
  df <- data.frame(ensembl_id = "ENSG001", gene_symbol = "GENE1", stringsAsFactors = FALSE)
  .save_annotation(df, m_file, l_file, list(data_last_modified = "2026-01-01"))
  
  # Normal check -> hit
  cache_normal <- .manage_cache(m_file, l_file, "2026-01-01", force_update = FALSE)
  expect_false(cache_normal$download)
  
  # Forced check -> re-download
  cache_forced <- .manage_cache(m_file, l_file, "2026-01-01", force_update = TRUE)
  expect_true(cache_forced$download)
})

test_that("mapIDToSymbol preserves distinct input dataframe rows with identical IDs", {
  annot <- data.frame(
    ensembl_id = c("ENSG001", "ENSG002"),
    gene_symbol = c("TP53", "BRCA1"),
    stringsAsFactors = FALSE
  )

  # Multiple experimental measurements for the same gene ID
  exp_data <- data.frame(
    id = c("ENSG001", "ENSG001", "ENSG002"),
    sample = c("S1", "S2", "S1"),
    logFC = c(1.5, 2.3, -0.8),
    stringsAsFactors = FALSE
  )

  res_collapse <- suppressMessages(mapIDToSymbol(exp_data, col_id = "id", annotation = annot, multi_handling = "collapse"))
  expect_equal(nrow(res_collapse$mapped), 3L)
  expect_equal(res_collapse$mapped$logFC, c(1.5, 2.3, -0.8))
  expect_equal(res_collapse$mapped$gene_symbol, c("TP53", "TP53", "BRCA1"))

  res_unique <- suppressMessages(mapIDToSymbol(exp_data, col_id = "id", annotation = annot, multi_handling = "unique"))
  expect_equal(nrow(res_unique$mapped), 3L)
  expect_equal(res_unique$mapped$sample, c("S1", "S2", "S1"))
})

test_that("mapIDToSymbol multi_mapped report hides internal row IDs", {
  annot <- data.frame(
    ensembl_id = c("ENSG001", "ENSG001"),
    gene_symbol = c("GENE_A", "GENE_B"),
    stringsAsFactors = FALSE
  )

  result <- suppressMessages(mapIDToSymbol(
    "ENSG001",
    annotation = annot,
    multi_handling = "keepAll"
  ))

  expect_equal(nrow(result$multi_mapped), 1L)
  expect_false("row_id_internal" %in% names(result$multi_mapped))
  expect_named(result$multi_mapped, c("temp_id", "n_symbols", "symbols"))
})

test_that("mapSymbolToID preserves exact input row order even with synonym rescue", {
  annot <- data.frame(
    ensembl_id = c("ENSG001", "ENSG002", "ENSG003"),
    gene_symbol = c("TP53", "BRCA1", "MYC"),
    gene_synonyms = c("P53", "BRCAI", "C-MYC"),
    stringsAsFactors = FALSE
  )

  # Input where first element is a synonym and second is a primary symbol
  input_symbols <- c("P53", "BRCA1", "C-MYC")
  res <- suppressMessages(mapSymbolToID(input_symbols, annotation = annot, keep_unmapped = TRUE))

  expect_equal(res$mapped$query_gene_symbol, c("P53", "BRCA1", "C-MYC"))
  expect_equal(res$mapped$ensembl_id, c("ENSG001", "ENSG002", "ENSG003"))
  expect_equal(res$mapped$gene_symbol, c("TP53", "BRCA1", "MYC"))
})

test_that("mapSymbolToID returns primary and synonym-linked records", {
  annot <- data.frame(
    ensembl_id = c("ENSG_AS2", NA_character_),
    gene_symbol = c("FRMD6-AS2", "FRMD6-AS1"),
    gene_synonyms = c("FRMD6-AS1", NA_character_),
    ncbi_id = c("100874185", "145438"),
    hgnc_id = c("43637", "20129"),
    chromosome = c("14", "14"),
    gene_start = c("51379509", "51649516"),
    gene_end = c("51651744", "51651744"),
    strand = c("-1", "-1"),
    metadata_source = c("Ensembl", "NCBI"),
    stringsAsFactors = FALSE
  )

  kept <- suppressMessages(mapSymbolToID(
    "FRMD6-AS1", annotation = annot, multi_handling = "keepAll"
  ))$mapped
  expect_equal(nrow(kept), 2L)
  expect_setequal(kept$gene_symbol, c("FRMD6-AS1", "FRMD6-AS2"))
  expect_setequal(kept$metadata_source, c("NCBI", "Ensembl"))

  collapsed <- suppressMessages(mapSymbolToID(
    "FRMD6-AS1", annotation = annot, multi_handling = "collapse"
  ))$mapped
  expect_equal(nrow(collapsed), 1L)
  expect_equal(collapsed$ensembl_id, "ENSG_AS2")
  expect_setequal(strsplit(collapsed$ncbi_id, "|", fixed = TRUE)[[1]], c("100874185", "145438"))
  expect_setequal(strsplit(collapsed$metadata_source, "|", fixed = TRUE)[[1]], c("Ensembl", "NCBI"))
})

test_that("mapSymbolToID retains unmatched rows without a synonyms column", {
  annot <- data.frame(
    ensembl_id = "ENSG001",
    gene_symbol = "TP53",
    stringsAsFactors = FALSE
  )

  res <- suppressMessages(mapSymbolToID(
    c("TP53", "MISSING", NA_character_),
    annotation = annot,
    keep_unmapped = TRUE
  ))

  expect_equal(nrow(res$mapped), 3L)
  expect_equal(res$mapped$query_gene_symbol, c("TP53", "MISSING", NA_character_))
  expect_equal(res$mapped$ensembl_id, c("ENSG001", NA_character_, NA_character_))
  expect_equal(res$unmapped_final, "MISSING")
})

test_that("NCBI GFF parsing prefers gene-level features without losing loci", {
  gff <- data.frame(
    ncbi_id = c("1", "1", "2", "2", "1"),
    seqid = c("chr1", "chr1", "chr2", "chr2", "alt1"),
    type = c("gene", "exon", "transcript", "exon", "exon"),
    stringsAsFactors = FALSE
  )

  retained <- .prefer_gene_level_features(gff)

  expect_equal(retained$type, c("gene", "transcript", "exon", "exon"))
  expect_true(all(c("1", "2") %in% retained$ncbi_id))
  expect_true("alt1" %in% retained$seqid)
})

test_that("mapIDToSymbol handles lowercase, whitespace, and versioned Ensembl IDs", {
  annot <- data.frame(
    ensembl_id = c("ENSG00000141510", "ENSG00000012048"),
    gene_symbol = c("TP53", "BRCA1"),
    stringsAsFactors = FALSE
  )

  input_ids <- c("  ensg00000141510.5  ", "ENSG00000012048.12")
  res <- suppressMessages(mapIDToSymbol(input_ids, annotation = annot))

  expect_equal(res$mapped$gene_symbol, c("TP53", "BRCA1"))
})

test_that("mapping functions validate multi_handling and annotation columns", {
  annot <- data.frame(
    ensembl_id = "ENSG001",
    gene_symbol = "TP53",
    stringsAsFactors = FALSE
  )

  # Invalid multi_handling value
  expect_error(
    suppressMessages(mapSymbolToID("TP53", annotation = annot, multi_handling = "invalid_option")),
    "should be one of"
  )
  expect_error(
    suppressMessages(mapIDToSymbol("ENSG001", annotation = annot, multi_handling = "invalid_option")),
    "should be one of"
  )

  # Missing required column in annotation
  bad_annot <- data.frame(gene_name = "TP53", stringsAsFactors = FALSE)
  expect_error(
    suppressMessages(mapSymbolToID("TP53", annotation = bad_annot)),
    "missing required column"
  )
  expect_error(
    suppressMessages(mapIDToSymbol("ENSG001", annotation = bad_annot)),
    "must contain column"
  )
})

test_that("mapIDToSymbol defaults to keep_unmapped = FALSE and handles keep_unmapped = TRUE", {
  annot <- data.frame(
    ensembl_id = c("ENSG001"),
    gene_symbol = c("TP53"),
    stringsAsFactors = FALSE
  )

  input_ids <- c("ENSG001", "ENSG_UNMAPPED")

  # Default: keep_unmapped = FALSE
  res_default <- suppressMessages(mapIDToSymbol(input_ids, annotation = annot))
  expect_equal(nrow(res_default$mapped), 1L)
  expect_equal(res_default$mapped$ensembl_id, "ENSG001")
  expect_equal(res_default$unmapped, "ENSG_UNMAPPED")

  # Explicit: keep_unmapped = TRUE
  res_kept <- suppressMessages(mapIDToSymbol(input_ids, annotation = annot, keep_unmapped = TRUE))
  expect_equal(nrow(res_kept$mapped), 2L)
  expect_true("ENSG_UNMAPPED" %in% res_kept$mapped$ensembl_id)
})

