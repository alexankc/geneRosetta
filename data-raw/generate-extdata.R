# data-raw/generate-extdata.R
# Maintainer script to generate, validate, compress, and manifest bundled reference extdata.
#
# Usage:
#   source("data-raw/generate-extdata.R")
#   generate_bundled_extdata(species = "all", output_dir = "inst/extdata")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(digest)
})

#' Generate and Compress Bundled Extdata for geneRosetta
#'
#' @param species Character vector of species ("human", "mouse", "rat") or "all"
#'   to update all species simultaneously. Defaults to "all".
#' @param output_dir Character string path to destination extdata directory.
#'   Must be explicitly provided.
#' @param ensembl_version Character string specifying Ensembl release to pin across
#'   all species, or "current". Defaults to "current".
#' @param force_update Logical; if TRUE, re-downloads all live source data.
#' @param compress Logical; if TRUE, compresses TSV files to .tsv.gz (gzip level 9).
#' @param remove_uncompressed Logical; if TRUE, removes uncompressed TSVs in output_dir.
#'
#' @return A list containing the generated manifest and status summary.
generate_bundled_extdata <- function(
    species = c("all", "human", "mouse", "rat"),
    output_dir = NULL,
    ensembl_version = "current",
    force_update = TRUE,
    compress = TRUE,
    remove_uncompressed = TRUE
) {
  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("An explicit 'output_dir' must be specified (e.g., 'inst/extdata').")
  }

  if ("all" %in% species) {
    species <- c("human", "mouse", "rat")
  } else {
    species <- match.arg(species, c("human", "mouse", "rat"), several.ok = TRUE)
  }

  cli::cli_h1("Generating Bundled Extdata for geneRosetta")
  cli::cli_alert_info("Target species: {paste(species, collapse = ', ')}")
  cli::cli_alert_info("Target output directory: {.file {output_dir}}")

  # Create a clean temporary staging directory
  staging_dir <- tempfile("generosetta_extdata_staging_")
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(staging_dir, recursive = TRUE), add = TRUE)

  # Define expected structure and minimum row count thresholds
  min_rows_thresholds <- list(
    human_assembly_report.tsv   = 100L,
    mouse_assembly_report.tsv   = 50L,
    rat_assembly_report.tsv     = 30L,
    human_biomart_mapping.tsv   = 50000L,
    mouse_biomart_mapping.tsv   = 40000L,
    rat_biomart_mapping.tsv     = 25000L,
    human_hgnc_mapping.tsv      = 35000L,
    mouse_mgi_mapping.tsv       = 45000L,
    human_ncbi_mapping.tsv      = 40000L,
    mouse_ncbi_mapping.tsv      = 35000L,
    rat_ncbi_mapping.tsv        = 25000L,
    rat_rgd_mapping.tsv         = 35000L
  )

  required_columns <- list(
    Assembly_Reports = c("sequence_name", "sequence_role", "assigned_molecule", "genbank_accn", "refseq_accn"),
    Ensembl_Genes    = c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome", "gene_start", "gene_end", "strand"),
    HGNC_human_Genes = c("hgnc_id", "gene_symbol", "ncbi_id", "ensembl_id"),
    MGI_mouse_genes  = c("mgi_id", "gene_symbol", "chromosome", "ensembl_id", "ncbi_id"),
    NCBI_Genes       = c("ncbi_id", "gene_symbol", "ensembl_id", "gene_biotype", "chromosome", "gene_start", "gene_end", "strand"),
    RGD_rat_genes    = c("rgd_id", "gene_symbol", "ensembl_id", "ncbi_id", "chromosome", "gene_start", "gene_end", "strand")
  )

  # Explicit allowlist of target subdirectories and files
  expected_files <- list()

  for (sp in species) {
    cli::cli_h2("Fetching and processing data for: {.val {sp}}")

    # 1. Assembly Report
    cli::cli_alert_info("Fetching assembly report for {sp}...")
    geneRosetta:::.get_assembly_report(sp, annotation_dir = staging_dir)
    expected_files[[length(expected_files) + 1]] <- list(
      sub_dir = "Assembly_Reports",
      tsv = paste0(sp, "_assembly_report.tsv"),
      log = paste0(sp, "_assembly_report_log.json"),
      type = "Assembly_Reports"
    )

    # 2. Ensembl BioMart
    cli::cli_alert_info("Fetching Ensembl BioMart for {sp} (version = {ensembl_version})...")
    geneRosetta::getAnnotationBiomart(
      species = sp,
      annotation_dir = staging_dir,
      version = ensembl_version,
      chromosome_scope = "all",
      force_update = force_update
    )
    expected_files[[length(expected_files) + 1]] <- list(
      sub_dir = "Ensembl_Genes",
      tsv = paste0(sp, "_biomart_mapping.tsv"),
      log = paste0(sp, "_biomart_mapping_log.json"),
      type = "Ensembl_Genes"
    )

    # 3. NCBI Genes
    cli::cli_alert_info("Fetching NCBI Genes for {sp}...")
    geneRosetta::getAnnotationNCBI(
      species = sp,
      annotation_dir = staging_dir,
      force_update = force_update
    )
    expected_files[[length(expected_files) + 1]] <- list(
      sub_dir = "NCBI_Genes",
      tsv = paste0(sp, "_ncbi_mapping.tsv"),
      log = paste0(sp, "_ncbi_mapping_log.json"),
      type = "NCBI_Genes"
    )

    # 4. Species-specific Authority
    if (sp == "human") {
      cli::cli_alert_info("Fetching HGNC human genes...")
      geneRosetta::getAnnotationHGNChuman(annotation_dir = staging_dir, force_update = force_update)
      expected_files[[length(expected_files) + 1]] <- list(
        sub_dir = "HGNC_human_Genes",
        tsv = "human_hgnc_mapping.tsv",
        log = "human_hgnc_mapping_log.json",
        type = "HGNC_human_Genes"
      )
    } else if (sp == "mouse") {
      cli::cli_alert_info("Fetching MGI mouse genes...")
      geneRosetta::getAnnotationMGImouse(annotation_dir = staging_dir, force_update = force_update)
      expected_files[[length(expected_files) + 1]] <- list(
        sub_dir = "MGI_mouse_genes",
        tsv = "mouse_mgi_mapping.tsv",
        log = "mouse_mgi_mapping_log.json",
        type = "MGI_mouse_genes"
      )
    } else if (sp == "rat") {
      cli::cli_alert_info("Fetching RGD rat genes...")
      geneRosetta::getAnnotationRGDrat(annotation_dir = staging_dir, force_update = force_update)
      expected_files[[length(expected_files) + 1]] <- list(
        sub_dir = "RGD_rat_genes",
        tsv = "rat_rgd_mapping.tsv",
        log = "rat_rgd_mapping_log.json",
        type = "RGD_rat_genes"
      )
    }
  }

  cli::cli_h2("Validating Generated Datasets")

  manifest_entries <- list()

  for (item in expected_files) {
    tsv_path <- file.path(staging_dir, item$sub_dir, item$tsv)
    log_path <- file.path(staging_dir, item$sub_dir, item$log)

    if (!file.exists(tsv_path)) {
      stop(sprintf("Expected TSV file missing in staging: %s", tsv_path))
    }
    if (!file.exists(log_path)) {
      stop(sprintf("Expected log file missing in staging: %s", log_path))
    }

    # Read and validate schema and row count
    df <- readr::read_tsv(tsv_path, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE)
    min_r <- min_rows_thresholds[[item$tsv]] %||% 10L
    if (nrow(df) < min_r) {
      stop(sprintf("Validation failed for %s: got %d rows, expected at least %d.", item$tsv, nrow(df), min_r))
    }

    req_cols <- required_columns[[item$type]]
    missing_c <- setdiff(req_cols, colnames(df))
    if (length(missing_c) > 0) {
      stop(sprintf("Validation failed for %s: missing required columns: %s", item$tsv, paste(missing_c, collapse = ", ")))
    }

    raw_bytes <- file.info(tsv_path)$size
    raw_sha256 <- digest::digest(file = tsv_path, algo = "sha256")

    # Compress if requested
    final_file_name <- item$tsv
    dest_sub_dir <- file.path(output_dir, item$sub_dir)
    dir.create(dest_sub_dir, recursive = TRUE, showWarnings = FALSE)

    if (isTRUE(compress)) {
      gz_file_name <- paste0(item$tsv, ".gz")
      gz_staging_path <- file.path(staging_dir, item$sub_dir, gz_file_name)
      geneRosetta:::.compress_gz(tsv_path, gz_staging_path, level = 9)

      gz_bytes <- file.info(gz_staging_path)$size
      gz_sha256 <- digest::digest(file = gz_staging_path, algo = "sha256")

      # Verify round-trip decompression
      verify_tmp <- tempfile(fileext = ".tsv")
      geneRosetta:::.decompress_gz(gz_staging_path, verify_tmp)
      verify_sha256 <- digest::digest(file = verify_tmp, algo = "sha256")
      unlink(verify_tmp)

      if (verify_sha256 != raw_sha256) {
        stop(sprintf("Decompression round-trip checksum mismatch for %s", item$tsv))
      }

      # Copy .gz to destination
      file.copy(gz_staging_path, file.path(dest_sub_dir, gz_file_name), overwrite = TRUE)
      if (!isTRUE(remove_uncompressed)) {
        file.copy(tsv_path, file.path(dest_sub_dir, item$tsv), overwrite = TRUE)
      } else {
        # Ensure any old uncompressed TSV in destination is removed
        old_tsv <- file.path(dest_sub_dir, item$tsv)
        if (file.exists(old_tsv)) unlink(old_tsv)
      }

      manifest_entries[[file.path(item$sub_dir, gz_file_name)]] <- list(
        rows = nrow(df),
        columns = colnames(df),
        uncompressed_size_bytes = raw_bytes,
        compressed_size_bytes = gz_bytes,
        uncompressed_sha256 = raw_sha256,
        compressed_sha256 = gz_sha256
      )
      cli::cli_alert_success("{item$sub_dir}/{gz_file_name}: {format(raw_bytes, big.mark=',')} B -> {format(gz_bytes, big.mark=',')} B ({nrow(df)} rows)")
    } else {
      file.copy(tsv_path, file.path(dest_sub_dir, item$tsv), overwrite = TRUE)
      manifest_entries[[file.path(item$sub_dir, item$tsv)]] <- list(
        rows = nrow(df),
        columns = colnames(df),
        uncompressed_size_bytes = raw_bytes,
        uncompressed_sha256 = raw_sha256
      )
      cli::cli_alert_success("{item$sub_dir}/{item$tsv}: {format(raw_bytes, big.mark=',')} B ({nrow(df)} rows)")
    }

    # Copy log file uncompressed
    file.copy(log_path, file.path(dest_sub_dir, item$log), overwrite = TRUE)
  }

  pkg_ver <- tryCatch(as.character(utils::packageVersion("geneRosetta")), error = function(e) {
    if (file.exists("DESCRIPTION")) unname(read.dcf("DESCRIPTION")[, "Version"]) else "1.0.0"
  })

  # Build manifest
  manifest <- list(
    snapshot_version = format(Sys.time(), "%Y-%m-%d"),
    created_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ"),
    format_version = 1L,
    species = species,
    r_version = R.version.string,
    package_version = pkg_ver,
    files = manifest_entries
  )

  manifest_path <- file.path(output_dir, "manifest.json")
  jsonlite::write_json(manifest, path = manifest_path, auto_unbox = TRUE, pretty = TRUE)
  cli::cli_alert_success("Manifest written to {.file {manifest_path}}")

  cli::cli_h1("Extdata Generation Completed Successfully")
  invisible(manifest)
}

# If invoked directly via Rscript:
if (!interactive() && identical(environment(), globalenv())) {
  args <- commandArgs(trailingOnly = TRUE)
  out_dir <- if (length(args) >= 1) args[1] else file.path("inst", "extdata")
  sp_arg <- if (length(args) >= 2) strsplit(args[2], ",")[[1]] else "all"
  generate_bundled_extdata(species = sp_arg, output_dir = out_dir)
}
