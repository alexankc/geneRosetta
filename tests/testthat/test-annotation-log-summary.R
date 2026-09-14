test_that("getAnnotationLogSummary exports manuscript-ready methods text", {
  annotation_dir <- tempfile("annotation-logs-")
  dir.create(annotation_dir)
  log_file <- file.path(annotation_dir, "human_merged_annotation_log.json")

  jsonlite::write_json(list(
    project_name = "geneRosetta human merged annotation",
    gene_count = 12345,
    source_details = list(
      Ensembl = list(
        data_access_timestamp = "2026-07-10T16:56:20Z",
        ensembl_version = "Ensembl_Genes_116",
        genome_build = "GRCh38.p14"
      ),
      NCBI = list(
        data_access_timestamp = "2026-08-26T14:36:52Z",
        data_last_modified = "2026-08-25"
      )
    )
  ), log_file, auto_unbox = TRUE)

  export_file <- file.path(annotation_dir, "methods", "human.txt")
  result <- getAnnotationLogSummary(
    species = "human",
    export_format = "methods",
    export_path = export_file,
    annotation_dir = annotation_dir
  )

  text <- readLines(export_file, warn = FALSE)
  expect_match(text, "Human gene annotations were generated")
  expect_match(text, "geneRosetta \\(v1.0.0;")
  expect_match(text, "Ensembl release 116")
  expect_match(text, "GRCh38.p14")
  expect_match(text, "NCBI Gene database \\(accessed August 2026\\)")
  expect_match(text, "HGNC-approved gene symbols were preferred as primary symbols")
  expect_match(text, "non-coding RNAs, and pseudogenes were retained")
  expect_match(text, "records lacking both an Ensembl and an NCBI identifier were excluded")
  expect_match(text, "12,345 unique entries", fixed = TRUE)
  expect_equal(result$gene_count, 12345)
})
