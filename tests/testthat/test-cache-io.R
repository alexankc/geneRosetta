# Regression tests for the cache round-trip and the NULL-lookup guards.
# All offline: they exercise the helpers directly, no network involved.

test_that("cached annotations survive a quote in a text field", {
  # .save_annotation() writes unquoted, so .manage_cache() must read unquoted.
  # Reading with read.delim()'s default quote treated a stray '"' as an opening
  # quote and silently swallowed every row after it.
  dir <- withr::local_tempdir()
  m_file <- file.path(dir, "quoted_mapping.tsv")
  l_file <- file.path(dir, "quoted_mapping_log.json")

  df <- data.frame(
    ensembl_id = paste0("ENSG", 1:5),
    gene_description = c(
      "plain one",
      'a 5" repeat region',
      "plain three",
      'flanked "both" sides',
      "plain five"
    ),
    stringsAsFactors = FALSE
  )

  .save_annotation(df, m_file, l_file, list(data_last_modified = "2024-01-01"))
  cache <- .manage_cache(m_file, l_file, "2024-01-01")

  expect_false(cache$download)
  expect_equal(nrow(cache$data), nrow(df))
  expect_equal(cache$data$ensembl_id, df$ensembl_id)
  expect_equal(cache$data$gene_description, df$gene_description)
})

test_that("a tab-free description with an apostrophe also round-trips", {
  dir <- withr::local_tempdir()
  m_file <- file.path(dir, "apos_mapping.tsv")
  l_file <- file.path(dir, "apos_mapping_log.json")

  df <- data.frame(
    ensembl_id = c("ENSG1", "ENSG2"),
    gene_description = c("5'-nucleotidase, cytosolic II",
                         "3' repair exonuclease"),
    stringsAsFactors = FALSE
  )

  .save_annotation(df, m_file, l_file, list(data_last_modified = "2024-01-01"))
  cache <- .manage_cache(m_file, l_file, "2024-01-01")

  expect_equal(cache$data$gene_description, df$gene_description)
})

test_that("cached annotations preserve embedded tabs and newlines", {
  dir <- withr::local_tempdir()
  m_file <- file.path(dir, "delimited_mapping.tsv")
  l_file <- file.path(dir, "delimited_mapping_log.json")

  df <- data.frame(
    ensembl_id = c("ENSG1", "ENSG2", "ENSG3"),
    gene_description = c("contains\ta tab", "contains\na newline", 'contains "quotes"'),
    stringsAsFactors = FALSE
  )

  .save_annotation(df, m_file, l_file, list(data_last_modified = "2024-01-01"))
  cache <- .manage_cache(m_file, l_file, "2024-01-01")

  expect_identical(cache$data, df)
})

test_that("a merged annotation reloads exactly as it was written", {
  # The builders emit an all-character table and use "" to mean "no symbol".
  # readr's default na = c("", "NA") turned those into NA on reload, so
  # buildAnnotation() returned different data on a cache hit than on the
  # initial build.
  combined <- data.frame(
    ensembl_id  = c("ENSG1", "ENSG2", "ENSG3"),
    gene_symbol = c("TP53", "", "MAFIP"),
    ncbi_id     = c("7157", NA, "727764"),
    chromosome  = c("17", "17", "chr14(GL000194.1)"),
    stringsAsFactors = FALSE
  )

  f <- withr::local_tempfile(fileext = ".tsv")
  readr::write_tsv(combined, f)

  reloaded <- readr::read_tsv(
    f,
    na = "NA",
    col_types = readr::cols(.default = readr::col_character())
  )

  # "" must stay "", NA must stay NA, and nothing may become numeric.
  expect_identical(reloaded$gene_symbol, combined$gene_symbol)
  expect_true(is.na(reloaded$ncbi_id[2]))
  expect_identical(as.data.frame(reloaded), combined)
})

test_that(".chromosome_role degrades to 'unknown' on a NULL lookup", {
  labels <- c("14", "chr14(GL000194.1)", "NOT_A_SCAFFOLD")

  expect_equal(
    .chromosome_role(labels, NULL),
    c("assembled-molecule", "unknown", "unknown")
  )
})

test_that(".filter_chromosome_scope keeps every row with no report", {
  # Without the report there is no way to separate an unplaced scaffold from an
  # alt-haplotype, so the safe answer is to drop nothing.
  df <- data.frame(
    gene = c("primary", "scaffold"),
    chromosome = c("14", "chr14(GL000194.1)"),
    stringsAsFactors = FALSE
  )

  # cli alerts are signalled as messages, not warning conditions.
  expect_message(out <- .filter_chromosome_scope(df, "standard", NULL))
  expect_equal(nrow(out), nrow(df))

  # "all" short-circuits before the lookup is ever consulted.
  expect_equal(nrow(.filter_chromosome_scope(df, "all", NULL)), nrow(df))
})

test_that("mapSymbolToID handles input with nothing mappable", {
  annot <- data.frame(
    ensembl_id = "ENSG001",
    gene_symbol = "TP53",
    gene_synonyms = "P53",
    stringsAsFactors = FALSE
  )

  # Zero-length input used to abort in the required-column backfill.
  expect_no_error(
    res <- suppressMessages(mapSymbolToID(character(0), annotation = annot))
  )
  expect_equal(nrow(res$mapped), 0L)
})

test_that("mapIDToSymbol works on an annotation without a synonyms column", {
  annot <- data.frame(
    ensembl_id = c("ENSG001", "ENSG002"),
    gene_symbol = c("TP53", "BRCA1"),
    stringsAsFactors = FALSE
  )

  expect_no_error(
    res <- suppressMessages(mapIDToSymbol("ENSG001", annotation = annot))
  )
  expect_equal(res$mapped$gene_symbol[1], "TP53")
})
