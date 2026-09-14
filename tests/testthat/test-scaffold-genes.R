# Live tests: they hit Ensembl BioMart, NCBI, HGNC and the NCBI assembly reports.
# Skipped on CRAN and when offline. All tests share one cache directory so the
# sources are downloaded once for the whole file.

cache <- new.env(parent = emptyenv())

live_annotation_dir <- function() {
  skip_on_cran()
  skip_if_offline()

  if (is.null(cache$dir)) {
    cache$dir <- tempfile("generosetta-live-")
    dir.create(cache$dir, recursive = TRUE)
    withr::defer(unlink(cache$dir, recursive = TRUE), teardown_env())
  }
  cache$dir
}

live_human_merged <- function(annotation_dir) {
  if (is.null(cache$merged)) {
    cache$merged <- buildHumanAnnotation(annotation_dir = annotation_dir)
  }
  cache$merged
}

test_that("getAnnotationBiomart places MAFIP on chr14(GL000194.1)", {
  human <- getAnnotationBiomart("human", annotation_dir = live_annotation_dir())
  mafip <- human[human$gene_symbol == "MAFIP", ]

  expect_equal(nrow(mafip), 1L)
  expect_equal(mafip$ensembl_id, "ENSG00000274847")
  expect_equal(mafip$chromosome, "chr14(GL000194.1)")
})

test_that("chromosome_scope narrows the table without a second download", {
  annotation_dir <- live_annotation_dir()

  all_genes <- getAnnotationBiomart("human", annotation_dir = annotation_dir)
  primary <- getAnnotationBiomart("human",
    annotation_dir = annotation_dir, chromosome_scope = "primary"
  )
  standard <- getAnnotationBiomart("human",
    annotation_dir = annotation_dir, chromosome_scope = "standard"
  )

  expect_gt(nrow(all_genes), nrow(primary))
  expect_gt(nrow(primary), nrow(standard))

  # MAFIP sits on a primary-assembly scaffold, so only "standard" may drop it.
  expect_true("MAFIP" %in% all_genes$gene_symbol)
  expect_true("MAFIP" %in% primary$gene_symbol)
  expect_false("MAFIP" %in% standard$gene_symbol)

  # "standard" must reproduce the pre-scaffold behaviour exactly.
  expect_setequal(
    unique(standard$chromosome),
    c(as.character(1:22), "X", "Y", "MT")
  )

  # "primary" adds unlocalized and unplaced scaffolds, never alt loci or patches.
  added <- setdiff(primary$chromosome, standard$chromosome)
  expect_true(all(grepl("^chr([0-9]+|X|Y|Un)\\(.+\\)$", added)))
})

test_that("getAnnotationHGNChuman includes alternate reference locus genes", {
  hgnc <- getAnnotationHGNChuman(annotation_dir = live_annotation_dir())

  expect_true("MAFIP" %in% hgnc$gene_symbol)
  # Only present in hgnc_complete_set, not in non_alt_loci_set.
  expect_true(all(c("HLA-DRB3", "HLA-DRB4") %in% hgnc$gene_symbol))
})

test_that("buildHumanAnnotation carries MAFIP through the merge with its scaffold", {
  human <- live_human_merged(live_annotation_dir())
  mafip <- human[human$gene_symbol == "MAFIP", ]

  expect_equal(nrow(mafip), 1L)
  expect_equal(mafip$ensembl_id, "ENSG00000274847")
  expect_equal(mafip$chromosome, "chr14(GL000194.1)")
  # The scaffold gene must still pick up its cross-references.
  expect_equal(mafip$hgnc_id, "31102")
  expect_equal(mafip$ncbi_id, "727764")
})

test_that("buildHumanAnnotation keeps every BioMart gene and the merged schema", {
  annotation_dir <- live_annotation_dir()
  human <- live_human_merged(annotation_dir)
  biomart <- getAnnotationBiomart("human", annotation_dir = annotation_dir)

  expect_true(all(unique(biomart$ensembl_id) %in% unique(human$ensembl_id)))
  expect_named(
    human,
    c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome", "gene_start",
      "gene_end", "strand", "metadata_source", "gene_synonyms", "ncbi_id", "hgnc_id",
      "gene_description")
  )
})

test_that("buildHumanAnnotation records each source log it merged", {
  annotation_dir <- live_annotation_dir()
  live_human_merged(annotation_dir)

  log_file <- file.path(annotation_dir, "human_merged_annotation_log.json")
  expect_true(file.exists(log_file))

  sources <- jsonlite::read_json(log_file)$source_details
  # A source that resolved records its metadata; a missing one records a string.
  expect_true(all(vapply(sources, is.list, logical(1))))
  expect_named(sources, c("Ensembl", "NCBI", "HGNC"))
})

test_that("mapSymbolToID resolves MAFIP against the merged annotation", {
  human <- live_human_merged(live_annotation_dir())
  result <- mapSymbolToID("MAFIP", annotation = human)

  expect_equal(result$mapped$ensembl_id, "ENSG00000274847")
  expect_equal(result$mapped$chromosome, "chr14(GL000194.1)")
})
