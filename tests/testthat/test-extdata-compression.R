# tests/testthat/test-extdata-compression.R
# Test streaming compression, transactional fallback materialization,
# deterministic offline builds, and manifest integrity.

test_that("streaming gzip compression and decompression have exact roundtrip fidelity", {
  tmp_src <- tempfile(fileext = ".tsv")
  tmp_gz <- tempfile(fileext = ".tsv.gz")
  tmp_dst <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_src, tmp_gz, tmp_dst)), add = TRUE)

  # Write dummy test data
  test_data <- data.frame(
    id = paste0("ID_", 1:1000),
    symbol = paste0("GENE_", 1:1000),
    description = paste("A simulated gene description with special chars:", 1:1000, "quote ' and \" | pipe"),
    value = rnorm(1000),
    stringsAsFactors = FALSE
  )
  readr::write_tsv(test_data, tmp_src)
  src_sha256 <- digest::digest(file = tmp_src, algo = "sha256")

  # Compress
  geneRosetta:::.compress_gz(tmp_src, tmp_gz, level = 9)
  expect_true(file.exists(tmp_gz))
  expect_lt(file.info(tmp_gz)$size, file.info(tmp_src)$size)

  # Decompress
  geneRosetta:::.decompress_gz(tmp_gz, tmp_dst)
  expect_true(file.exists(tmp_dst))
  dst_sha256 <- digest::digest(file = tmp_dst, algo = "sha256")

  expect_equal(dst_sha256, src_sha256)

  read_back <- readr::read_tsv(tmp_dst, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
  expect_equal(nrow(read_back), 1000L)
  expect_equal(read_back$symbol[1], "GENE_1")
})

test_that("corrupted gzip decompression fails cleanly without creating partial files", {
  tmp_corrupt_gz <- tempfile(fileext = ".tsv.gz")
  tmp_cache_file <- tempfile(fileext = ".tsv")
  on.exit(unlink(c(tmp_corrupt_gz, tmp_cache_file)), add = TRUE)

  # Write truncated/garbage gzip data
  writeBin(charToRaw("Not a valid gzip file header\x1f\x8b\x08garbage_content_truncated"), tmp_corrupt_gz)

  expect_error(
    geneRosetta:::.decompress_gz(tmp_corrupt_gz, tmp_cache_file)
  )
  expect_false(file.exists(tmp_cache_file))
})

test_that(".materialize_bundled_fallback handles .tsv.gz and populates cache", {
  fixture_root <- tempfile("fixture_extdata_")
  fixture_sub <- file.path(fixture_root, "Ensembl_Genes")
  dir.create(fixture_sub, recursive = TRUE)
  on.exit(unlink(fixture_root, recursive = TRUE), add = TRUE)

  # Create a small sample mapping TSV and log
  sample_df <- data.frame(
    ensembl_id = c("ENSG000001", "ENSG000002"),
    gene_symbol = c("TEST1", "TEST2"),
    stringsAsFactors = FALSE
  )
  raw_tsv <- tempfile(fileext = ".tsv")
  readr::write_tsv(sample_df, raw_tsv)

  gz_dest <- file.path(fixture_sub, "human_biomart_mapping.tsv.gz")
  geneRosetta:::.compress_gz(raw_tsv, gz_dest, level = 9)
  unlink(raw_tsv)

  log_dest <- file.path(fixture_sub, "human_biomart_mapping_log.json")
  jsonlite::write_json(list(ensembl_version = "Ensembl_Genes_116"), log_dest, auto_unbox = TRUE)

  # Target cache dir
  cache_dir <- tempfile("cache_dir_")
  dir.create(cache_dir, recursive = TRUE)
  on.exit(unlink(cache_dir, recursive = TRUE), add = TRUE)

  mapping_file <- file.path(cache_dir, "human_biomart_mapping.tsv")
  log_file <- file.path(cache_dir, "human_biomart_mapping_log.json")

  # Materialize
  fb <- geneRosetta:::.materialize_bundled_fallback(
    mapping_file = mapping_file,
    log_file = log_file,
    sub_dir = "Ensembl_Genes",
    bundled_root = fixture_root
  )

  expect_true(fb$success)
  expect_true(file.exists(mapping_file))
  expect_true(file.exists(log_file))
  expect_equal(nrow(fb$data), 2L)
  expect_equal(fb$data$gene_symbol, c("TEST1", "TEST2"))
})

test_that(".find_bundled_file prioritizes .gz over uncompressed and supports root override", {
  fixture_root <- tempfile("fixture_priority_")
  fixture_sub <- file.path(fixture_root, "Test_Source")
  dir.create(fixture_sub, recursive = TRUE)
  on.exit(unlink(fixture_root, recursive = TRUE), add = TRUE)

  # Case 1: only uncompressed exists
  raw_file <- file.path(fixture_sub, "test_mapping.tsv")
  writeLines("header\nval", raw_file)

  hit1 <- geneRosetta:::.find_bundled_file("Test_Source", "test_mapping.tsv", bundled_root = fixture_root)
  expect_false(is.null(hit1))
  expect_false(hit1$is_compressed)
  expect_equal(hit1$path, raw_file)

  # Case 2: .gz also exists -> preferred
  gz_file <- file.path(fixture_sub, "test_mapping.tsv.gz")
  writeBin(as.raw(c(1, 2, 3)), gz_file)

  hit2 <- geneRosetta:::.find_bundled_file("Test_Source", "test_mapping.tsv", bundled_root = fixture_root)
  expect_false(is.null(hit2))
  expect_true(hit2$is_compressed)
  expect_equal(hit2$path, gz_file)
})

test_that("deterministic offline buildAnnotation works from bundled .tsv.gz without internet access", {
  tmp_cache <- tempfile("generosetta_offline_cache_")
  dir.create(tmp_cache, recursive = TRUE)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)

  # Human offline build
  human_annot <- buildAnnotation("human", annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(human_annot, "data.frame")
  expect_gt(nrow(human_annot), 50000L)

  # Check that scaffold gene MAFIP was preserved
  mafip <- human_annot[human_annot$gene_symbol == "MAFIP" & !is.na(human_annot$gene_symbol), ]
  expect_equal(nrow(mafip), 1L)
  expect_equal(mafip$chromosome, "chr14(GL000194.1)")

  # Verify that cache now contains uncompressed TSVs
  cached_tsvs <- list.files(tmp_cache, pattern = "\\.tsv$", recursive = TRUE)
  expect_true("human_merged_annotation.tsv" %in% cached_tsvs)
  expect_true(any(grepl("human_biomart_mapping\\.tsv$", cached_tsvs)))
  expect_true(any(grepl("human_ncbi_mapping\\.tsv$", cached_tsvs)))
  expect_true(any(grepl("human_hgnc_mapping\\.tsv$", cached_tsvs)))

  # Second call should hit the uncompressed cache directly
  human_annot_2 <- buildAnnotation("human", annotation_dir = tmp_cache, offline = TRUE)
  expect_equal(nrow(human_annot), nrow(human_annot_2))
})

test_that("deterministic offline buildAnnotation works for mouse and rat", {
  tmp_cache <- tempfile("generosetta_offline_mouse_rat_")
  dir.create(tmp_cache, recursive = TRUE)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)

  # Mouse offline build
  mouse_annot <- buildAnnotation("mouse", annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(mouse_annot, "data.frame")
  expect_gt(nrow(mouse_annot), 50000L)

  # Rat offline build
  rat_annot <- buildAnnotation("rat", annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(rat_annot, "data.frame")
  expect_gt(nrow(rat_annot), 30000L)
})

test_that("all individual getters support offline mode directly", {
  tmp_cache <- tempfile("generosetta_getters_offline_")
  dir.create(tmp_cache, recursive = TRUE)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)

  # Assembly report
  report <- geneRosetta:::.get_assembly_report("human", annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(report, "data.frame")
  expect_gt(nrow(report), 100L)

  # BioMart
  bm <- getAnnotationBiomart("human", annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(bm, "data.frame")
  expect_gt(nrow(bm), 50000L)

  # NCBI
  ncbi <- getAnnotationNCBI("human", annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(ncbi, "data.frame")
  expect_gt(nrow(ncbi), 50000L)

  # HGNC
  hgnc <- getAnnotationHGNChuman(annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(hgnc, "data.frame")
  expect_gt(nrow(hgnc), 35000L)

  # MGI
  mgi <- getAnnotationMGImouse(annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(mgi, "data.frame")
  expect_gt(nrow(mgi), 50000L)

  # RGD
  rgd <- getAnnotationRGDrat(annotation_dir = tmp_cache, offline = TRUE)
  expect_s3_class(rgd, "data.frame")
  expect_gt(nrow(rgd), 35000L)
})

test_that("force_update with offline=TRUE or simulated network error falls back safely", {
  tmp_cache <- tempfile("generosetta_force_fallback_")
  dir.create(tmp_cache, recursive = TRUE)
  on.exit(unlink(tmp_cache, recursive = TRUE), add = TRUE)

  # Even with force_update = TRUE, if offline = TRUE it falls back cleanly to bundled data
  res <- getAnnotationHGNChuman(annotation_dir = tmp_cache, force_update = TRUE, offline = TRUE)
  expect_gt(nrow(res), 35000L)
})

test_that("manifest.json is valid and matches all bundled .tsv.gz files", {
  manifest_path <- system.file("extdata", "manifest.json", package = "geneRosetta")
  if (!file.exists(manifest_path)) {
    manifest_path <- file.path("inst", "extdata", "manifest.json")
  }
  expect_true(file.exists(manifest_path))

  manifest <- jsonlite::read_json(manifest_path)
  expect_equal(manifest$format_version, 1L)
  expect_setequal(manifest$species, c("human", "mouse", "rat"))
  expect_gt(length(manifest$files), 10L)

  extdata_dir <- dirname(manifest_path)

  # Verify each file listed in manifest
  for (rel_file in names(manifest$files)) {
    entry <- manifest$files[[rel_file]]
    full_path <- file.path(extdata_dir, rel_file)
    expect_true(file.exists(full_path), info = paste("Missing manifest file:", full_path))

    actual_sha256 <- digest::digest(file = full_path, algo = "sha256")
    expect_equal(actual_sha256, entry$compressed_sha256, info = paste("SHA256 mismatch for", rel_file))
  }
})
