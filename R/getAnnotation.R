#' Get Remote File Modification Date
#'
#' Retrieves the last modification date of a remote file from its HTTP headers.
#' Useful for checking if cached data needs updating.
#'
#' @param url Character string specifying the URL to check. If NULL, returns NULL.
#'
#' @return A Date object representing the last modification time, or NA if the request fails.
#'   Returns NA if the URL is NULL or if the remote server is unreachable.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#'   .get_remote_date("https://example.com/data.txt")
#' }
.get_remote_date <- function(url) {
  if (is.null(url)) return(NULL)
  tryCatch({
    h <- httr::HEAD(url, httr::timeout(10))
    lm <- httr::headers(h)[["last-modified"]]
    if (is.null(lm)) return(as.Date(Sys.time()))
    
    old_locale <- Sys.getlocale("LC_TIME")
    on.exit(Sys.setlocale("LC_TIME", old_locale))
    Sys.setlocale("LC_TIME", "C") 
    as.Date(as.POSIXct(lm, format = "%a, %d %b %Y %H:%M:%S", tz = "GMT"))
  }, error = function(e) NA)
}

#' Get Persistent Data Path
#' @return Character string of the user data directory
.get_data_path <- function() {
  path <- tools::R_user_dir("geneRosetta", which = "data")
  if (!dir.exists(path)) dir.create(path, recursive = TRUE)
  return(path)
}

# Read cached TSVs with the same NA and quoting conventions used by
# .save_annotation(). All fields stay character until a public getter applies
# its schema-specific conversions.
.read_annotation_tsv <- function(path) {
  as.data.frame(readr::read_tsv(
    path,
    na = "NA",
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ), stringsAsFactors = FALSE)
}


#' Default Root Directory for Bundled Extdata
#'
#' @param override Optional directory path for tests or fixtures.
#' @return Character string path to the bundled extdata root directory.
#' @keywords internal
.bundled_extdata_root <- function(override = NULL) {
  if (!is.null(override) && dir.exists(override)) return(override)
  opt <- getOption("geneRosetta.bundled_extdata_root")
  if (!is.null(opt) && dir.exists(opt)) return(opt)
  system.file("extdata", package = "geneRosetta")
}

#' Compress a File Using Streaming Gzip
#'
#' @param src_file Path to source uncompressed file.
#' @param gz_file Path to destination .gz file.
#' @param level Gzip compression level (1-9). Defaults to 9.
#' @return The destination path, invisibly.
#' @keywords internal
.compress_gz <- function(src_file, gz_file, level = 9) {
  if (!file.exists(src_file)) {
    cli::cli_abort("Source file {.file {src_file}} does not exist.")
  }
  dir.create(dirname(gz_file), recursive = TRUE, showWarnings = FALSE)

  con_in <- file(src_file, "rb")
  on.exit(close(con_in), add = TRUE)

  con_out <- gzfile(gz_file, "wb", compression = level)
  on.exit(close(con_out), add = TRUE)

  while (length(buf <- readBin(con_in, "raw", n = 1048576L)) > 0L) {
    writeBin(buf, con_out)
  }
  invisible(gz_file)
}

#' Streaming Gzip Decompressor
#'
#' @param gz_file Path to input .gz file.
#' @param dest_file Path to output decompressed file.
#' @return The destination file path invisibly.
#' @keywords internal
.decompress_gz <- function(gz_file, dest_file) {
  if (!file.exists(gz_file)) {
    cli::cli_abort("Compressed file {.file {gz_file}} does not exist.")
  }
  dir.create(dirname(dest_file), recursive = TRUE, showWarnings = FALSE)

  # Check magic bytes for gzip (\x1f\x8b)
  raw_header <- readBin(gz_file, "raw", n = 2)
  if (length(raw_header) < 2 || raw_header[1] != as.raw(0x1f) || raw_header[2] != as.raw(0x8b)) {
    cli::cli_abort("File {.file {gz_file}} is not a valid gzip file (invalid magic number).")
  }

  tmp_dest <- tempfile(pattern = "decompress_", tmpdir = dirname(dest_file), fileext = ".tmp")
  on.exit(if (file.exists(tmp_dest)) unlink(tmp_dest), add = TRUE)

  con_in <- gzfile(gz_file, "rb")
  con_out <- file(tmp_dest, "wb")

  bytes_written <- 0L
  tryCatch({
    while (length(buf <- readBin(con_in, "raw", n = 1048576L)) > 0L) {
      writeBin(buf, con_out)
      bytes_written <- bytes_written + length(buf)
    }
  }, finally = {
    close(con_in)
    close(con_out)
  })

  if (bytes_written == 0L) {
    cli::cli_abort("Decompressed file {.file {dest_file}} is empty.")
  }

  if (file.exists(dest_file)) unlink(dest_file)
  if (!file.rename(tmp_dest, dest_file)) {
    file.copy(tmp_dest, dest_file, overwrite = TRUE)
    unlink(tmp_dest)
  }

  invisible(dest_file)
}

#' Locate a Bundled File Within Extdata (preferring .gz over raw)
#'
#' @param sub_dir Subdirectory within extdata (e.g., "Ensembl_Genes").
#' @param file_name Base file name to find (e.g., "human_biomart_mapping.tsv").
#' @param bundled_root Root extdata directory. Defaults to `.bundled_extdata_root()`.
#' @return A list with `path` and `is_compressed`, or NULL if not found.
#' @keywords internal
.find_bundled_file <- function(sub_dir, file_name, bundled_root = .bundled_extdata_root()) {
  if (is.null(bundled_root) || !nzchar(bundled_root) || !dir.exists(bundled_root)) {
    return(NULL)
  }

  base_clean <- sub("\\.gz$", "", file_name)
  target_dir <- if (!is.null(sub_dir) && nzchar(sub_dir)) file.path(bundled_root, sub_dir) else bundled_root

  candidate_gz <- file.path(target_dir, paste0(base_clean, ".gz"))
  candidate_raw <- file.path(target_dir, base_clean)

  if (file.exists(candidate_gz)) {
    return(list(path = candidate_gz, is_compressed = TRUE))
  }
  if (file.exists(candidate_raw)) {
    return(list(path = candidate_raw, is_compressed = FALSE))
  }
  return(NULL)
}

#' Transactionally Materialize Bundled Fallback into Local Cache
#'
#' @param mapping_file Destination mapping TSV path.
#' @param log_file Destination JSON log path.
#' @param sub_dir Subdirectory within extdata.
#' @param bundled_root Root extdata directory.
#' @return A list with `data` (tibble) and `success` (logical).
#' @keywords internal
.materialize_bundled_fallback <- function(mapping_file, log_file, sub_dir = NULL, bundled_root = .bundled_extdata_root()) {
  bundled_map <- .find_bundled_file(sub_dir, basename(mapping_file), bundled_root = bundled_root)
  if (is.null(bundled_map)) {
    return(list(data = NULL, success = FALSE))
  }

  cli::cli_alert_info("Materializing bundled reference data into cache (may be outdated).")
  dest_dir <- dirname(mapping_file)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  # Decompress/copy transactionally using a temporary file in the same directory (enables atomic rename)
  tmp_dest <- tempfile(pattern = "tmp_map_", tmpdir = dest_dir, fileext = ".tsv")

  decompress_ok <- tryCatch({
    if (bundled_map$is_compressed) {
      .decompress_gz(bundled_map$path, tmp_dest)
    } else {
      file.copy(bundled_map$path, tmp_dest, overwrite = TRUE)
    }
    TRUE
  }, error = function(e) {
    unlink(tmp_dest)
    FALSE
  })

  if (!decompress_ok || !file.exists(tmp_dest) || file.info(tmp_dest)$size == 0) {
    if (file.exists(tmp_dest)) unlink(tmp_dest)
    cli::cli_abort("Failed to decompress bundled fallback file {.file {bundled_map$path}}.")
  }

  # Atomic rename
  renamed <- file.rename(tmp_dest, mapping_file)
  if (!renamed) {
    file.copy(tmp_dest, mapping_file, overwrite = TRUE)
    unlink(tmp_dest)
  }

  # Copy/decompress log file if available
  bundled_log <- .find_bundled_file(sub_dir, basename(log_file), bundled_root = bundled_root)
  if (!is.null(bundled_log)) {
    if (bundled_log$is_compressed) {
      .decompress_gz(bundled_log$path, log_file)
    } else {
      file.copy(bundled_log$path, log_file, overwrite = TRUE)
    }
  } else {
    cli::cli_warn("Bundled metadata log not found \u2014 cache may be incomplete.")
  }

  data <- .read_annotation_tsv(mapping_file)
  return(list(data = data, success = TRUE))
}

#' Manage Caching and Offline Fallback
#'
#' Checks if cached data is available and valid by comparing local metadata against 
#' a remote timestamp. Implements a three-tier priority logic:
#' 1. Remote: If online and cache is stale, signals a re-download.
#' 2. Local Cache: If offline or cache is fresh, returns user-downloaded data.
#' 3. Bundled Data: If no local cache exists and offline, returns internal package data.
#'
#' @param mapping_file Character string path to the cached mapping TSV file.
#' @param log_file Character string path to the JSON log file containing metadata.
#' @param remote_stamp Date/Character representing the remote data's version/date.
#' @param sub_dir Character string specifying the subdirectory within \code{inst/extdata} 
#'   where bundled fallback data resides.
#' @param force_update Logical; if \code{TRUE}, bypasses the local cache and forces a re-download.
#' @param exact_version Logical; if \code{TRUE}, requires exact equality between cached and remote version.
#' @param offline Logical; if \code{TRUE}, forces offline mode and uses cached/bundled data.
#' @param bundled_root Path to root extdata directory.
#'
#' @return A list with two elements:
#'   \item{data}{A dataframe (from local cache or bundled internal data) 
#'     or NULL if a re-download is required.}
#'   \item{download}{Logical; \code{TRUE} if the server is reachable and 
#'     the remote data is newer than the cache, \code{FALSE} otherwise.}
#'
#' @keywords internal
.manage_cache <- function(mapping_file, log_file, remote_stamp, sub_dir = NULL, force_update = FALSE, exact_version = FALSE, offline = getOption("geneRosetta.offline", FALSE), bundled_root = .bundled_extdata_root()) {
  if (isTRUE(force_update) && !isTRUE(offline)) {
    return(list(data = NULL, download = TRUE))
  }

  cache_exists <- file.exists(mapping_file) && file.exists(log_file)
  
  # Offline or server unreachable: check local cache or fallback
  if (isTRUE(offline) || is.null(remote_stamp) || is.na(remote_stamp)) {
    if (cache_exists) {
      cli::cli_alert_info("Offline: Using cached version from {basename(mapping_file)}")
      return(list(
        data = .read_annotation_tsv(mapping_file),
        download = FALSE
      ))
    }
    
    # Materialize bundled fallback
    fb <- .materialize_bundled_fallback(mapping_file, log_file, sub_dir = sub_dir, bundled_root = bundled_root)
    if (fb$success) {
      return(list(data = fb$data, download = FALSE))
    } else {
      cli::cli_abort(c(
        "x" = "Internet connection unavailable and no cached/bundled data found.",
        "i" = "Please connect to the internet for the initial data download."
      ))
    }
  }

  # Online: check if local cache matches remote version/date
  if (cache_exists) {
    log_info <- tryCatch(jsonlite::read_json(log_file), error = function(e) NULL)
    if (!is.null(log_info)) {
      cached_val <- log_info$ensembl_version %||% log_info$data_last_modified
      remote_val <- as.character(remote_stamp)
      
      if (!is.null(cached_val)) {
        is_up_to_date <- .compare_versions(cached_val, remote_val, exact = exact_version)
        
        if (is_up_to_date) {
          cli::cli_alert_success("Cache is up-to-date ({cached_val}). Loading...")
          return(list(
            data = .read_annotation_tsv(mapping_file),
            download = FALSE
          ))
        }
      }
    }
  }
  
  return(list(data = NULL, download = TRUE))
}


.compare_versions <- function(cached, remote, exact = FALSE) {
  if (exact) {
    return(identical(as.character(cached), as.character(remote)))
  }

  # If both are Ensembl version strings
  if (grepl("Ensembl_Genes_", cached) && grepl("Ensembl_Genes_", remote)) {
    cached_num <- as.numeric(gsub(".*_(\\d+)$", "\\1", cached))
    remote_num <- as.numeric(gsub(".*_(\\d+)$", "\\1", remote))
    return(!is.na(cached_num) && !is.na(remote_num) && cached_num >= remote_num)
  }
  
  # Otherwise treat as dates/strings
  return(cached >= remote)
}

#' Save Annotation Data and Metadata
#'
#' Writes downloaded annotation data to a TSV file and saves associated metadata
#' (including download timestamp and source information) to a JSON log file.
#'
#' @param df A dataframe containing the gene annotation data to save.
#' @param mapping_file Character string path where the TSV file will be written.
#' @param log_file Character string path where the JSON metadata log will be written.
#' @param metadata A named list of metadata to include in the JSON log
#'   (e.g., version info, URLs, package versions).
#'
#' @return NULL (invisibly). Called for its side effect of writing files to disk.
#'
#' @keywords internal
.save_annotation <- function(df, mapping_file, log_file, metadata) {
  dir.create(dirname(mapping_file), recursive = TRUE, showWarnings = FALSE)
  readr::write_tsv(df, mapping_file, na = "NA", quote = "needed")
  
  log_info <- c(
    list(data_access_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")),
    metadata
  )
  jsonlite::write_json(log_info, path = log_file, auto_unbox = TRUE, pretty = TRUE)
  cli::cli_alert_success("Data and metadata log saved to disk.")
}

#' Get Gene Annotation from Ensembl BioMart
#'
#' Downloads and processes gene annotation data from Ensembl BioMart for human,
#' mouse, or rat species. Includes caching and offline fallback functionality.
#'
#' @param species Character string specifying the species. One of "human" (default),
#'   "mouse", or "rat".
#' @param annotation_dir Character string specifying the directory where annotation
#'   data will be cached. Defaults to the persistent user data directory.
#' @param version Character string specifying the Ensembl version to use.
#'   Defaults to "current" for the latest Ensembl version.
#' @param chromosome_scope Character string controlling which sequences are kept:
#'   \itemize{
#'     \item \code{"all"} (default): every gene Ensembl reports, including those on
#'       alternate haplotypes and patch scaffolds.
#'     \item \code{"primary"}: the primary assembly only, i.e. assembled molecules
#'       plus unlocalized and unplaced scaffolds.
#'     \item \code{"standard"}: assembled molecules only (1-22, X, Y, MT).
#'   }
#'   Filtering is applied on the way out, so a single cache serves every scope.
#' @param force_update Logical; if \code{TRUE}, forces re-download and re-parsing.
#' @param offline Logical; if \code{TRUE}, forces offline mode and uses cached/bundled data.
#'
#' @return A dataframe with the following columns:
#'   \item{ensembl_id}{Ensembl gene identifier}
#'   \item{gene_symbol}{Official gene symbol}
#'   \item{gene_biotype}{Biotype classification (protein_coding, lncRNA, etc.)}
#'   \item{chromosome}{Chromosome location (see Details)}
#'   \item{gene_start}{Genomic start position}
#'   \item{gene_end}{Genomic end position}
#'   \item{strand}{Strand (+ or -)}
#'   \item{gene_description}{Gene description}
#'   \item{gene_synonyms}{Pipe-separated alternative gene names}
#'
#' @details
#' Assembled molecules are reported as bare names ("14", "X", "MT"). Every other
#' sequence is labelled \code{chr<molecule>(<GenBank accession>)} using the NCBI
#' assembly report for the build Ensembl is serving, so the unlocalized scaffold
#' carrying \code{MAFIP} is reported as \code{chr14(GL000194.1)}. Scaffolds with no
#' assigned molecule are labelled \code{chrUn(<accession>)}.
#'
#' Multiple synonyms are combined and separated by pipes (|).
#' Results are cached locally with metadata about the Ensembl version
#' and genome build used.
#'
#' @examples
#' \dontrun{
#'   human_genes <- getAnnotationBiomart("human")
#'   mouse_genes <- getAnnotationBiomart("mouse")
#'
#'   # Drop alternate haplotypes and patches, keep scaffolds like chr14(GL000194.1)
#'   primary_only <- getAnnotationBiomart("human", chromosome_scope = "primary")
#' }
#'
#' @import dplyr
#' @export
#'
getAnnotationBiomart <- function(
    species = "human",
    annotation_dir = .get_data_path(),
    version = "current",
    chromosome_scope = c("all", "primary", "standard"),
    force_update = FALSE,
    offline = getOption("geneRosetta.offline", FALSE)
) {
  chromosome_scope <- match.arg(chromosome_scope)
  species <- tolower(species)

  dataset <- switch(
    species,
    human = "hsapiens_gene_ensembl",
    mouse = "mmusculus_gene_ensembl",
    rat   = "rnorvegicus_gene_ensembl",
    cli::cli_abort("Unsupported species: {.val {species}}")
  )

  folder <- file.path(annotation_dir, "Ensembl_Genes")
  version_suffix <- if (version != "current") paste0("_v", version) else ""
  m_file <- file.path(folder, paste0(species, "_biomart_mapping", version_suffix, ".tsv"))
  l_file <- file.path(folder, paste0(species, "_biomart_mapping", version_suffix, "_log.json"))

  mart <- NULL
  if (!isTRUE(offline)) {
    cli::cli_progress_step("Connecting to Ensembl...")
    mart <- tryCatch({
      httr::with_config(httr::timeout(30), {
        v_arg <- if (version == "current") NULL else version
        biomaRt::useEnsembl(biomart = "genes", dataset = dataset, version = v_arg)
      })
    }, error = function(e) NULL)
  }

  ensembl_version <- if (!is.null(mart)) {
    tryCatch({
      v <- biomaRt::listEnsembl()$version[1]
      if (version != "current") paste0("Ensembl_Genes_", version) else gsub(" ", "_", v)
    }, error = function(e) NA)
  } else {
    if (version != "current") paste0("Ensembl_Genes_", version) else NA
  }

  genome_build <- if (!is.null(mart)) {
    tryCatch({
      dataset_list <- biomaRt::listDatasets(mart)
      as.character(dataset_list$version[dataset_list$dataset == dataset])
    }, error = function(e) NA_character_)
  } else {
    NA_character_
  }

  cache <- .manage_cache(
    m_file, l_file, ensembl_version,
    sub_dir = "Ensembl_Genes",
    force_update = force_update,
    exact_version = (version != "current"),
    offline = offline
  )

  if (!cache$download) {
    cached_df <- as_tibble(cache$data) %>%
      mutate(
        across(any_of(c("gene_start", "gene_end", "strand")), as.integer),
        across(any_of(c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome", "gene_synonyms", "gene_description")), as.character)
      )
    if (chromosome_scope != "all") {
      lookup <- .chromosome_lookup(.get_assembly_report(species, annotation_dir, genome_build, offline = offline))
      cached_df <- .filter_chromosome_scope(cached_df, chromosome_scope, lookup)
    }
    return(cached_df)
  }

  lookup <- .chromosome_lookup(.get_assembly_report(species, annotation_dir, genome_build, offline = offline))

  cli::cli_progress_step("Connecting to Ensembl and querying...")

  if (is.null(mart)) {
    v_arg <- if (version == "current") NULL else version
    mart <- tryCatch({
      biomaRt::useEnsembl(biomart = "genes", dataset = dataset, version = v_arg)
    }, error = function(e) NULL)
  }

  attrs <- c(
    "ensembl_gene_id", "external_gene_name", "gene_biotype",
    "external_synonym", "chromosome_name", "start_position",
    "end_position", "strand", "description"
  )

  bm_df <- if (!is.null(mart)) {
    tryCatch({
      biomaRt::getBM(attributes = attrs, mart = mart)
    }, error = function(e) NULL)
  } else {
    NULL
  }

  if (is.null(bm_df)) {
    if (file.exists(m_file)) {
      cli::cli_alert_warning("Ensembl is currently unresponsive; using cached data from {basename(m_file)}")
      cached_df <- as_tibble(.read_annotation_tsv(m_file)) %>%
        mutate(
          across(any_of(c("gene_start", "gene_end", "strand")), as.integer),
          across(any_of(c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome", "gene_synonyms", "gene_description")), as.character)
        )
      if (chromosome_scope != "all") {
        lookup <- .chromosome_lookup(.get_assembly_report(species, annotation_dir, genome_build, offline = offline))
        cached_df <- .filter_chromosome_scope(cached_df, chromosome_scope, lookup)
      }
      return(cached_df)
    }
    
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "Ensembl_Genes")
    if (fb$success) {
      cli::cli_alert_warning("Ensembl is currently unresponsive; using bundled package data.")
      cached_df <- as_tibble(fb$data) %>%
        mutate(
          across(any_of(c("gene_start", "gene_end", "strand")), as.integer),
          across(any_of(c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome", "gene_synonyms", "gene_description")), as.character)
        )
      if (chromosome_scope != "all") {
        lookup <- .chromosome_lookup(.get_assembly_report(species, annotation_dir, genome_build, offline = offline))
        cached_df <- .filter_chromosome_scope(cached_df, chromosome_scope, lookup)
      }
      return(cached_df)
    }

    cli::cli_abort("Unable to connect to Ensembl and no cached or bundled fallback data found.")
  }

  mapping_df <- bm_df %>%
    mutate(external_synonym = if_else(external_synonym == "", NA_character_, external_synonym)) %>%
    group_by(
      ensembl_id = ensembl_gene_id,
      gene_symbol = external_gene_name,
      gene_biotype,
      chromosome = .format_chromosome(chromosome_name, lookup),
      gene_start = start_position,
      gene_end = end_position,
      strand,
      gene_description = description
    ) %>%
    summarise(
      gene_synonyms = {
        vals <- unique(na.omit(external_synonym))
        if (length(vals) == 0L) NA_character_ else paste(vals, collapse = "|")
      },
      .groups = "drop"
    ) %>%
    as_tibble()

  .save_annotation(mapping_df, m_file, l_file, list(
    ensembl_version = ensembl_version,
    genome_build = genome_build,
    mart_dataset = dataset,
    biomart_pkg_version = as.character(utils::packageVersion("biomaRt"))
  ))

  return(.filter_chromosome_scope(mapping_df, chromosome_scope, lookup))
}

#' Get Gene Annotation from NCBI Gene Database
#'
#' Downloads and processes gene annotation data from NCBI Gene database for
#' human, mouse, or rat species. Includes caching and offline fallback.
#'
#' @param species Character string specifying the species. One of "human" (default),
#'   "mouse", or "rat".
#' @param annotation_dir Character string specifying the directory where annotation
#'   data will be cached. Defaults to the persistent user data directory.
#' @param force_update Logical; if \code{TRUE}, forces re-download and re-parsing.
#' @param offline Logical; if \code{TRUE}, forces offline mode and uses cached/bundled data.
#'
#' @return A dataframe with the following columns:
#'   \item{ncbi_id}{NCBI Gene ID (GeneID)}
#'   \item{gene_symbol}{Official gene symbol}
#'   \item{gene_synonyms}{Pipe-separated alternative gene names}
#'   \item{ensembl_id}{Corresponding Ensembl gene identifier (if available)}
#'   \item{gene_biotype}{NCBI gene type from \code{type_of_gene}}
#'   \item{gene_description}{NCBI gene description}
#'
#' @details
#' Ensembl IDs are extracted from the dbXrefs field where available.
#' Multiple synonyms are combined and separated by pipes (|).
#' Data is cached with modification timestamps for incremental updates.
#'
#' @examples
#' \dontrun{
#'   human_genes <- getAnnotationNCBI("human")
#'   mouse_genes <- getAnnotationNCBI("mouse")
#' }
#'
#' @import dplyr
#' @export
#' 
getAnnotationNCBI <- function(
    species = "human",
    annotation_dir = .get_data_path(),
    force_update = FALSE,
    offline = getOption("geneRosetta.offline", FALSE)
) {
  species <- tolower(species)

  urls <- list(
    human = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz",
    mouse = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Mus_musculus.gene_info.gz",
    rat   = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Rattus_norvegicus.gene_info.gz"
  )

  query_url <- urls[[species]] %||% cli::cli_abort("Unsupported species: {.val {species}}")
  
  folder <- file.path(annotation_dir, "NCBI_Genes")
  m_file <- file.path(folder, paste0(species, "_ncbi_mapping.tsv"))
  l_file <- file.path(folder, paste0(species, "_ncbi_mapping_log.json"))

  remote_stamp <- if (isTRUE(offline)) NA else .get_remote_date(query_url)
  cache <- .manage_cache(m_file, l_file, remote_stamp, sub_dir = "NCBI_Genes", force_update = force_update, offline = offline)
  if (!cache$download) {
    cached_df <- as_tibble(cache$data) %>%
      mutate(
        across(any_of(c("gene_start", "gene_end", "strand")), as.integer),
        across(any_of(c("ncbi_id", "gene_symbol", "gene_synonyms", "ensembl_id", "chromosome", "gene_biotype", "gene_description")), as.character)
      )
    required_gene_info_cols <- c("gene_biotype", "gene_description")
    if (all(required_gene_info_cols %in% names(cached_df))) {
      return(cached_df)
    }
    if (!isTRUE(offline)) {
      cli::cli_alert_warning(
        "Cached NCBI annotation predates retained gene type/description fields; refreshing it."
      )
    } else {
      return(cached_df)
    }
  }

  cli::cli_progress_step("Downloading and parsing NCBI data...")
  
  mapping_df <- tryCatch({
    response <- httr::GET(query_url)
    if (httr::status_code(response) != 200) {
      cli::cli_abort("Failed to download from NCBI: HTTP {httr::status_code(response)}")
    }

    temp_file <- tempfile(fileext = ".gz")
    on.exit(unlink(temp_file), add = TRUE)
    writeBin(httr::content(response, "raw"), temp_file)
    
    parsed_df <- read.delim(
      temp_file, 
      header = TRUE, 
      sep = "\t", 
      stringsAsFactors = FALSE, 
      quote = ""
    )

    clean_df <- parsed_df %>%
      mutate(
        ensembl_id = sapply(strsplit(as.character(dbXrefs), "\\|"), function(x) {
          ids <- grep("^Ensembl:", x, value = TRUE)
          if (length(ids) > 0) sub("^Ensembl:", "", ids[1]) else NA_character_
        })
      ) %>%
      dplyr::select(
        ncbi_id = GeneID, 
        gene_symbol = Symbol, 
        gene_synonyms = Synonyms, 
        ensembl_id,
        gene_biotype = type_of_gene,
        gene_description = description
      ) %>%
      dplyr::mutate(ncbi_id = as.character(ncbi_id)) %>%
      as_tibble()

    coords_df <- .getAnnotationNCBICoordsGFF(species, annotation_dir, force_update = force_update, offline = offline)
    
    clean_df %>% dplyr::left_join(coords_df, by = "ncbi_id")
  }, error = function(e) {
    if (file.exists(m_file)) {
      cli::cli_alert_warning("NCBI download failed; using cached data from {basename(m_file)}")
      return(as_tibble(.read_annotation_tsv(m_file)))
    }
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "NCBI_Genes")
    if (fb$success) {
      cli::cli_alert_warning("NCBI download failed; using bundled package data.")
      return(as_tibble(fb$data))
    }
    cli::cli_abort("Failed to download from NCBI and no cached/bundled data found: {conditionMessage(e)}")
  })

  if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
    .save_annotation(mapping_df, m_file, l_file, list(
      data_last_modified = as.character(remote_stamp),
      data_url = query_url,
      httr_pkg_version = as.character(utils::packageVersion("httr"))
    ))
  }

  return(mapping_df)
}
  
#' Get Human Gene Annotation from HGNC Database
#'
#' Downloads and processes human gene annotation data from the HGNC
#' (Human Genome Nomenclature Committee) database. Includes caching and
#' offline fallback functionality.
#'
#' @param annotation_dir Character string specifying the directory where annotation
#'   data will be cached. Defaults to the persistent user data directory.
#' @param force_update Logical; if \code{TRUE}, forces re-download and re-parsing.
#' @param offline Logical; if \code{TRUE}, forces offline mode and uses cached/bundled data.
#'
#' @return A dataframe with the following columns for human genes:
#'   \item{hgnc_id}{HGNC identifier}
#'   \item{gene_symbol}{Official HGNC gene symbol}
#'   \item{gene_synonyms}{Pipe-separated alternative gene names}
#'   \item{ncbi_id}{NCBI Gene ID (Entrez ID)}
#'   \item{ensembl_id}{Ensembl gene identifier}
#'
#' @details
#' The function downloads the complete HGNC dataset in TSV format.
#' Only essential mapping columns are retained.
#' Multiple synonyms and previous symbols are merged and separated by pipes (|).
#'
#' @examples
#' \dontrun{
#'   hgnc_genes <- getAnnotationHGNChuman()
#' }
#'
#' @import dplyr
#' @export
#' 
getAnnotationHGNChuman <- function(
    annotation_dir = .get_data_path(),
    force_update = FALSE,
    offline = getOption("geneRosetta.offline", FALSE)
) {
  hgnc_url <- "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/hgnc_complete_set.txt"
  folder   <- file.path(annotation_dir, "HGNC_human_Genes")
  m_file   <- file.path(folder, "human_hgnc_mapping.tsv")
  l_file   <- file.path(folder, "human_hgnc_mapping_log.json")

  remote_stamp <- if (isTRUE(offline)) NA else .get_remote_date(hgnc_url)
  
  cache <- .manage_cache(m_file, l_file, remote_stamp, sub_dir = "HGNC_human_Genes", force_update = force_update, offline = offline)
  if (!cache$download) {
    return(as_tibble(cache$data) %>%
      mutate(across(everything(), as.character)))
  }

  cli::cli_progress_step("Downloading and processing HGNC data...")
  
  mapping_df <- tryCatch({
    tmp <- tempfile()
    on.exit(unlink(tmp), add = TRUE)
    download.file(hgnc_url, destfile = tmp, mode = "wb", quiet = TRUE)
    
    df_raw <- read.delim(tmp, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")

    df_raw %>%
      dplyr::select(
        hgnc_id,
        gene_symbol = symbol,
        alias_symbol,
        prev_symbol,
        ncbi_id = entrez_id,
        ensembl_id = ensembl_gene_id
      ) %>%
      mutate(hgnc_id = sub("^HGNC:", "", as.character(hgnc_id))) %>%
      tidyr::unite("gene_synonyms", alias_symbol, prev_symbol, sep = "|", na.rm = TRUE) %>%
      mutate(gene_synonyms = purrr::map_chr(stringi::stri_split_fixed(gene_synonyms, "|"), function(x) {
        x <- stringi::stri_trim_both(x)
        x <- x[x != "" & !is.na(x) & x != "NA" & x != "-"]
        unique_syns <- stringi::stri_sort(unique(x))
        if (length(unique_syns) == 0) return(NA_character_)
        stringi::stri_join(unique_syns, collapse = "|")
      })) %>%
      mutate(across(everything(), as.character)) %>%
      as_tibble()
  }, error = function(e) {
    if (file.exists(m_file)) {
      cli::cli_alert_warning("HGNC download failed; using cached data from {basename(m_file)}")
      return(as_tibble(.read_annotation_tsv(m_file)))
    }
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "HGNC_human_Genes")
    if (fb$success) {
      cli::cli_alert_warning("HGNC download failed; using bundled package data.")
      return(as_tibble(fb$data))
    }
    cli::cli_abort("Failed to download from HGNC and no cached/bundled data found: {conditionMessage(e)}")
  })

  if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
    .save_annotation(mapping_df, m_file, l_file, list(
      data_last_modified = as.character(remote_stamp),
      data_url = hgnc_url
    ))
  }

  return(mapping_df)
}

#' Get Mouse Gene Annotation from MGI Database
#'
#' Downloads and processes mouse gene annotation data from the Mouse Genome
#' Informatics (MGI) database. Combines genetic markers with coordinate data.
#'
#' @param annotation_dir Character string specifying the directory where annotation
#'   data will be cached. Defaults to the persistent user data directory.
#' @param offline Logical. If TRUE, uses cached/bundled package data instead of downloading.
#'   Defaults to FALSE.
#' @param force_update Logical; if \code{TRUE}, forces re-download and re-parsing.
#'
#' @return A dataframe with the following columns for mouse genes:
#'   \item{mgi_id}{MGI marker accession ID}
#'   \item{gene_symbol}{Official gene symbol}
#'   \item{gene_synonyms}{Pipe-separated alternative gene names}
#'   \item{ensembl_id}{Ensembl gene identifier}
#'   \item{ncbi_id}{NCBI Gene ID (Entrez ID)}
#'   \item{chromosome}{Chromosome location}
#'   \item{gene_start}{Genomic start position (cM)}
#'   \item{gene_end}{Genomic end position (cM)}
#'   \item{strand}{Strand (+ or -)}
#'   \item{gene_biotype}{Marker type (e.g., "Gene", "Pseudogene")}
#'   \item{gene_description}{Marker name/description}
#'
#' @details
#' The function downloads two files from MGI:
#' \itemize{
#'   \item MRK_List2.rpt: Marker catalog with symbols, names, biotypes, and synonyms
#'   \item MGI_Gene_Model_Coord.rpt: Genomic coordinates and Ensembl/NCBI cross-references
#' }
#' Data is combined and deduplicated. Only markers with type "Gene" or "Pseudogene"
#' are retained.
#'
#' @examples
#' \dontrun{
#'   mouse_genes <- getAnnotationMGImouse()
#' }
#'
#' @import dplyr
#' @export
#' 
getAnnotationMGImouse <- function(
    annotation_dir = .get_data_path(),
    offline = getOption("geneRosetta.offline", FALSE),
    force_update = FALSE
) {
  base_url <- "https://www.informatics.jax.org/downloads/reports/"
  url_catalog <- paste0(base_url, "MRK_List2.rpt")
  folder <- file.path(annotation_dir, "MGI_mouse_genes")
  m_file <- file.path(folder, "mouse_mgi_mapping.tsv")
  l_file <- file.path(folder, "mouse_mgi_mapping_log.json")

  remote_stamp <- if (isTRUE(offline)) NA else .get_remote_date(url_catalog)
  cache <- .manage_cache(m_file, l_file, remote_stamp, sub_dir = "MGI_mouse_genes", force_update = force_update, offline = offline)
  
  if (!cache$download) {
    return(as_tibble(cache$data) %>%
      mutate(across(everything(), as.character)))
  }

  cli::cli_progress_step("Downloading and merging MGI tables...")
  
  mapping_df <- tryCatch({
    tmp_cat <- tempfile()
    tmp_crd <- tempfile()
    on.exit(unlink(c(tmp_cat, tmp_crd)), add = TRUE)
    download.file(url_catalog, destfile = tmp_cat, mode = "wb", quiet = TRUE)
    download.file(paste0(base_url, "MGI_Gene_Model_Coord.rpt"), destfile = tmp_crd, mode = "wb", quiet = TRUE)

    lookup <- .try_chromosome_lookup("mouse", annotation_dir, offline = offline)

    # Clean coordinates
    mouse_coords <- read.delim(tmp_crd, header = TRUE, sep = "\t", check.names = TRUE, row.names = NULL, stringsAsFactors = FALSE, quote = "")
    colnames(mouse_coords) <- c(colnames(mouse_coords)[-1], "extra_column") # Fix MGI header shift
    clean_names <- gsub("^X\\d+\\.\\.", "", colnames(mouse_coords)) %>%
                   gsub("\\.", "_", .)
    colnames(mouse_coords) <- clean_names
    mouse_coords$extra_column <- NULL

    # Clean markers
    mouse_markers <- read.delim(tmp_cat, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "")
    names(mouse_markers) <- gsub("\\.", "_", names(mouse_markers))

    # Merge and rename
    mouse_coords %>%
      left_join(mouse_markers, by = c("MGI_accession_id" = "MGI_Accession_ID")) %>%
      dplyr::select(
        mgi_id = MGI_accession_id,
        gene_symbol = marker_symbol,
        gene_biotype = Feature_Type,
        chromosome = Ensembl_gene_chromosome,
        gene_start = Ensembl_gene_start,
        gene_end = Ensembl_gene_end,
        strand = Ensembl_gene_strand,
        gene_synonyms = Marker_Synonyms__pipe_separated_,
        ncbi_id = Entrez_gene_id,
        ensembl_id = Ensembl_gene_id,
        gene_description = marker_name
      ) %>%
      mutate(
        mgi_id = stringr::str_remove(mgi_id, "^MGI:"),
        gene_synonyms = dplyr::na_if(gene_synonyms, ""),
        strand = case_when(strand == "+" ~ 1L, strand == "-" ~ -1L, TRUE ~ NA_integer_),
        chromosome = .format_chromosome(chromosome, lookup),
        ensembl_id = dplyr::na_if(ensembl_id, "null"),
        ncbi_id = dplyr::na_if(ncbi_id, "null"),
        ensembl_id = dplyr::na_if(ensembl_id, ""),
        ncbi_id = dplyr::na_if(ncbi_id, "")
      ) %>%
      mutate(across(everything(), as.character)) %>%
      as_tibble()
  }, error = function(e) {
    if (file.exists(m_file)) {
      cli::cli_alert_warning("MGI download failed; using cached data from {basename(m_file)}")
      return(as_tibble(.read_annotation_tsv(m_file)))
    }
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "MGI_mouse_genes")
    if (fb$success) {
      cli::cli_alert_warning("MGI download failed; using bundled package data.")
      return(as_tibble(fb$data))
    }
    cli::cli_abort("Failed to download from MGI and no cached/bundled data found: {conditionMessage(e)}")
  })

  if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
    .save_annotation(mapping_df, m_file, l_file, list(
      data_last_modified = as.character(remote_stamp),
      catalog_url = url_catalog
    ))
  }

  return(mapping_df)
}

#' Get Rat Gene Annotation from RGD Database
#'
#' Downloads and processes rat gene annotation data from the Rat Genome
#' Database (RGD). Handles assembly-specific data selection and cleaning.
#'
#' @param annotation_dir Character string specifying the directory where annotation
#'   data will be cached. Defaults to the persistent user data directory.
#' @param genome_build Character string specifying the rat genome build version.
#'   Defaults to "GRCr8". Must match column names in the RGD data file.
#' @param force_update Logical; if \code{TRUE}, forces re-download and re-parsing.
#' @param offline Logical; if \code{TRUE}, forces offline mode and uses cached/bundled data.
#'
#' @return A dataframe with the following columns for rat genes:
#'   \item{rgd_id}{RGD gene ID}
#'   \item{gene_symbol}{Official gene symbol}
#'   \item{gene_description}{Gene name/description}
#'   \item{gene_biotype}{Gene type classification}
#'   \item{chromosome}{Chromosome location}
#'   \item{gene_start}{Genomic start position}
#'   \item{gene_end}{Genomic end position}
#'   \item{strand}{Strand (+ or -)}
#'   \item{gene_synonyms}{Pipe-separated alternative gene names}
#'   \item{ensembl_id}{Ensembl gene identifier}
#'   \item{ncbi_id}{NCBI Gene ID (Entrez ID)}
#'
#' @details
#' The function downloads the primary GENES_RAT.txt file from RGD.
#' Column names are dynamically adjusted based on the specified \code{genome_build}.
#' Semicolon-separated Ensembl IDs are expanded into separate rows.
#'
#' @examples
#' \dontrun{
#'   rat_genes <- getAnnotationRGDrat(genome_build = "GRCr8")
#' }
#'
#' @import dplyr
#' @export
#' 
getAnnotationRGDrat <- function(
    annotation_dir = .get_data_path(),
    genome_build = "GRCr8",
    force_update = FALSE,
    offline = getOption("geneRosetta.offline", FALSE)
) {
  rgd_url <- "https://download.rgd.mcw.edu/data_release/GENES_RAT.txt"
  folder  <- file.path(annotation_dir, "RGD_rat_genes")
  m_file  <- file.path(folder, "rat_rgd_mapping.tsv")
  l_file  <- file.path(folder, "rat_rgd_mapping_log.json")

  remote_stamp <- if (isTRUE(offline)) NA else .get_remote_date(rgd_url)
  cache <- .manage_cache(m_file, l_file, remote_stamp, sub_dir = "RGD_rat_genes", force_update = force_update, offline = offline)
  
  if (!cache$download) {
    return(as_tibble(cache$data) %>%
      mutate(across(everything(), as.character)))
  }

  cli::cli_progress_step("Processing RGD table (Build: {genome_build})...")

  old_to <- getOption("timeout")
  on.exit(options(timeout = old_to), add = TRUE)
  options(timeout = 600)

  mapping_df <- tryCatch({
    tmp_rgd <- tempfile()
    on.exit(unlink(tmp_rgd), add = TRUE)
    download.file(rgd_url, destfile = tmp_rgd, mode = "wb", quiet = TRUE)

    lookup <- .try_chromosome_lookup("rat", annotation_dir, genome_build, offline = offline)

    read.delim(tmp_rgd, header = TRUE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE) %>%
      dplyr::select(
        rgd_id = GENE_RGD_ID, gene_symbol = SYMBOL, gene_synonyms = OLD_SYMBOL,
        gene_biotype = GENE_TYPE, ensembl_id = ENSEMBL_ID, ncbi_id = NCBI_GENE_ID, 
        gene_description = NAME,
        contains(genome_build)
      ) %>%
      rename_with(~ "chromosome", contains(paste0("CHROMOSOME_", genome_build))) %>%
      rename_with(~ "gene_start", contains(paste0("START_POS_", genome_build))) %>%
      rename_with(~ "gene_end",   contains(paste0("STOP_POS_", genome_build))) %>%
      rename_with(~ "strand",     contains(paste0("STRAND_", genome_build))) %>%
      tidyr::separate_rows(ensembl_id, sep = ";") %>%
      mutate(gene_synonyms = stringr::str_replace_all(gene_synonyms, ";\\s*", "|")) %>%
      mutate(across(everything(), as.character)) %>%
      mutate(chromosome = .format_chromosome(chromosome, lookup)) %>%
      as_tibble()
  }, error = function(e) {
    if (file.exists(m_file)) {
      cli::cli_alert_warning("RGD download failed; using cached data from {basename(m_file)}")
      return(as_tibble(.read_annotation_tsv(m_file)))
    }
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "RGD_rat_genes")
    if (fb$success) {
      cli::cli_alert_warning("RGD download failed; using bundled package data.")
      return(as_tibble(fb$data))
    }
    cli::cli_abort("Failed to download from RGD and no cached/bundled data found: {conditionMessage(e)}")
  })

  if (!is.null(mapping_df) && nrow(mapping_df) > 0) {
    .save_annotation(mapping_df, m_file, l_file, list(
      data_last_modified = as.character(remote_stamp),
      data_url = rgd_url,
      genome_build = genome_build
    ))
  }

  return(mapping_df)
}

#' Prefer Gene-Level GFF Features Without Dropping Unrepresented Loci
#'
#' @keywords internal
.prefer_gene_level_features <- function(df) {
  df %>%
    dplyr::group_by(ncbi_id, seqid) %>%
    dplyr::filter(
      !any(type %in% c("gene", "pseudogene")) |
        type %in% c("gene", "pseudogene")
    ) %>%
    dplyr::ungroup()
}

#' Get NCBI Genomic Coordinates from GFF3
#'
#' @param species Species name.
#' @param annotation_dir Directory path.
#' @param force_update Logical.
#' @param offline Logical.
#' @keywords internal
.getAnnotationNCBICoordsGFF <- function(species, annotation_dir, force_update = FALSE, offline = getOption("geneRosetta.offline", FALSE)) {
  folder <- file.path(annotation_dir, "NCBI_Genes")
  coords_file <- file.path(folder, paste0(species, "_ncbi_coords.tsv"))
  coords_log  <- file.path(folder, paste0(species, "_ncbi_coords_log.json"))

  if (file.exists(coords_file) && (!isTRUE(force_update) || isTRUE(offline))) {
    return(readr::read_tsv(
      coords_file,
      col_types = readr::cols(
        ncbi_id = readr::col_character(),
        chromosome = readr::col_character(),
        gene_start = readr::col_integer(),
        gene_end = readr::col_integer(),
        strand = readr::col_integer()
      ),
      show_col_types = FALSE
    ))
  }

  if (isTRUE(offline)) {
    return(tibble::tibble(
      ncbi_id = character(),
      chromosome = character(),
      gene_start = integer(),
      gene_end = integer(),
      strand = integer()
    ))
  }

  org_map <- list(
    human = "Homo_sapiens",
    mouse = "Mus_musculus",
    rat   = "Rattus_norvegicus"
  )
  org <- org_map[[species]]
  
  base_url <- paste0("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/", org, "/latest_assembly_versions/")
  listing <- tryCatch({
    resp <- httr::GET(base_url, httr::timeout(30))
    httr::content(resp, as = "text", encoding = "UTF-8")
  }, error = function(e) NULL)
  
  pattern <- "GCF_[0-9]+\\.[0-9]+_[A-Za-z0-9._-]+"
  hits <- unique(stringr::str_extract_all(listing, pattern)[[1]])

  preferred_build <- switch(
    species,
    human = "GRCh38",
    mouse = "GRCm39",
    rat   = "GRCr8"
  )

  preferred_hits <- hits[
    stringr::str_detect(hits, stringr::fixed(preferred_build))
  ]

  if (length(preferred_hits) == 0) {
    if (file.exists(coords_file)) {
      return(readr::read_tsv(coords_file, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE))
    }
    return(tibble::tibble(
      ncbi_id = character(),
      chromosome = character(),
      gene_start = integer(),
      gene_end = integer(),
      strand = integer()
    ))
  }

  dir_name <- preferred_hits[1]
  gff_url <- paste0(base_url, dir_name, "/", dir_name, "_genomic.gff.gz")

  cli::cli_progress_step(paste0("Downloading and parsing NCBI GFF3 (", dir_name, ")..."))
  temp_gff <- tempfile(fileext = ".gff.gz")
  on.exit(unlink(temp_gff), add = TRUE)
  
  old_to <- getOption("timeout")
  options(timeout = max(old_to, 600))
  on.exit(options(timeout = old_to), add = TRUE)
  
  download_res <- tryCatch({
    utils::download.file(gff_url, destfile = temp_gff, mode = "wb", quiet = TRUE)
    TRUE
  }, error = function(e) FALSE)

  if (!download_res || !file.exists(temp_gff) || file.info(temp_gff)$size == 0) {
    if (file.exists(coords_file)) {
      return(readr::read_tsv(coords_file, col_types = readr::cols(.default = readr::col_character()), show_col_types = FALSE))
    }
    return(tibble::tibble(
      ncbi_id = character(),
      chromosome = character(),
      gene_start = integer(),
      gene_end = integer(),
      strand = integer()
    ))
  }
  
  gff_cols <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
  gff_df <- utils::read.delim(
    temp_gff,
    sep = "\t",
    comment.char = "#",
    header = FALSE,
    col.names = gff_cols,
    stringsAsFactors = FALSE,
    quote = ""
  ) %>% as_tibble()
  
  gene_features <- gff_df %>%
    dplyr::filter(
      !is.na(attributes),
      stringr::str_detect(attributes, "Dbxref=GeneID:")
    ) %>%
    dplyr::mutate(
      ncbi_id = as.character(stringr::str_extract(attributes, "(?<=Dbxref=GeneID:)\\d+")),
      strand_num = dplyr::case_when(
        strand == "+" ~ 1L,
        strand == "-" ~ -1L,
        TRUE ~ NA_integer_
      )
    ) %>%
    dplyr::filter(!is.na(ncbi_id)) %>%
    .prefer_gene_level_features()

  locus_spans <- gene_features %>%
    dplyr::group_by(ncbi_id, seqid) %>%
    dplyr::summarise(
      gene_start = min(start, na.rm = TRUE),
      gene_end = max(end, na.rm = TRUE),
      strand = if (all(is.na(strand_num))) NA_integer_ else na.omit(strand_num)[1],
      .groups = "drop"
    )

  report <- .get_assembly_report(species, annotation_dir, offline = offline)
  refseq_to_base <- stats::setNames(
    ifelse(report$sequence_role == "assembled-molecule", report$sequence_name, report$genbank_accn),
    report$refseq_accn
  )
  lookup <- .try_chromosome_lookup(species, annotation_dir, offline = offline)
  
  coords_clean <- locus_spans %>%
    dplyr::mutate(
      base_acc = unname(refseq_to_base[seqid]),
      chromosome = .format_chromosome(base_acc, lookup),
      is_primary_chrom = !is.na(chromosome) & !grepl("\\(", chromosome),
      span_length = gene_end - gene_start
    ) %>%
    dplyr::arrange(ncbi_id, desc(is_primary_chrom), desc(span_length)) %>%
    dplyr::distinct(ncbi_id, .keep_all = TRUE) %>%
    dplyr::select(ncbi_id, chromosome, gene_start, gene_end, strand)
    
  if (!dir.exists(folder)) {
    dir.create(folder, recursive = TRUE, showWarnings = FALSE)
  }
  readr::write_tsv(coords_clean, coords_file)
  jsonlite::write_json(
    list(assembly = dir_name, timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ")),
    path = coords_log,
    auto_unbox = TRUE,
    pretty = TRUE
  )
  return(coords_clean)
}

#' Internal Maintainer Tooling: Generate and Compress Extdata
#'
#' @param species Character vector of species or "all" to update all species.
#' @param output_dir Path to destination extdata directory.
#' @param force_update Logical.
#' @param compress Logical.
#' @param remove_uncompressed Logical.
#' @keywords internal
#' @noRd
.generate_extdata <- function(
    species = c("all", "human", "mouse", "rat"),
    output_dir = NULL,
    force_update = TRUE,
    compress = TRUE,
    remove_uncompressed = TRUE
) {
  if (is.null(output_dir) || !nzchar(output_dir)) {
    stop("An explicit 'output_dir' must be specified.")
  }
  # Source data-raw generator if present, or invoke directly
  gen_script <- file.path("data-raw", "generate-extdata.R")
  if (file.exists(gen_script)) {
    source(gen_script, local = TRUE)
    generate_bundled_extdata(
      species = species,
      output_dir = output_dir,
      force_update = force_update,
      compress = compress,
      remove_uncompressed = remove_uncompressed
    )
  } else {
    stop("data-raw/generate-extdata.R not found.")
  }
}

