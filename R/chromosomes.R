#' Standard (Assembled-Molecule) Chromosome Names
#'
#' The chromosome tokens that Ensembl reports for assembled molecules across the
#' species supported by this package. Used as a fast path by
#' \code{.format_chromosome()}: anything in this set is passed through unchanged,
#' anything else is looked up in an NCBI assembly report.
#'
#' @keywords internal
.STANDARD_CHROMOSOMES <- c(as.character(1:22), "X", "Y", "MT")

#' Sequence Roles Retained by Each Chromosome Scope
#'
#' \code{"standard"} keeps only assembled molecules (the historic behaviour),
#' \code{"primary"} additionally keeps the unlocalized and unplaced scaffolds that
#' make up the rest of the primary assembly, and \code{"all"} keeps everything,
#' including alternate haplotypes (alt-scaffolds) and patch releases.
#'
#' @keywords internal
.CHROMOSOME_SCOPE_ROLES <- list(
  standard = "assembled-molecule",
  primary  = c("assembled-molecule", "unlocalized-scaffold", "unplaced-scaffold"),
  all      = NULL
)

#' NCBI GenBank Assembly Locations
#'
#' The GenBank accession prefix directory for each species on the NCBI genomes
#' FTP site, plus a pinned assembly directory used when the live listing cannot
#' be read (offline, or Ensembl has moved to a build NCBI has not published yet).
#'
#' @keywords internal
.ASSEMBLY_SPECIES <- list(
  human = list(root = "GCA/000/001/405", fallback = "GCA_000001405.29_GRCh38.p14"),
  mouse = list(root = "GCA/000/001/635", fallback = "GCA_000001635.9_GRCm39"),
  rat   = list(root = "GCA/036/323/735", fallback = "GCA_036323735.1_GRCr8")
)

#' Base URL of the NCBI genomes FTP tree
#' @keywords internal
.NCBI_GENOMES_URL <- "https://ftp.ncbi.nlm.nih.gov/genomes/all/"

#' Column names of an NCBI assembly report
#' @keywords internal
.ASSEMBLY_REPORT_COLS <- c(
  "sequence_name", "sequence_role", "assigned_molecule",
  "assigned_molecule_location_type", "genbank_accn", "relationship",
  "refseq_accn", "assembly_unit", "sequence_length", "ucsc_style_name"
)

#' Resolve the URL of an NCBI Assembly Report
#'
#' Ensembl reports its genome build as the same string NCBI uses to suffix the
#' assembly directory ("GRCh38.p14", "GRCm39", "GRCr8"), so the report for the
#' exact build Ensembl is serving can be located by listing the species'
#' accession directory. If the listing is unreachable or the build is absent,
#' falls back to the pinned assembly and warns.
#'
#' @param species Character string. One of "human", "mouse", or "rat".
#' @param genome_build Character string of the Ensembl genome build, or NULL/NA
#'   to go straight to the pinned assembly.
#'
#' @return Character string URL of the \code{*_assembly_report.txt} file.
#'
#' @keywords internal
.assembly_report_url <- function(species, genome_build = NULL) {
  cfg <- .ASSEMBLY_SPECIES[[tolower(species)]]
  if (is.null(cfg)) cli::cli_abort("Unsupported species: {.val {species}}")

  dir_name <- NULL
  have_build <- !is.null(genome_build) && !is.na(genome_build) && nzchar(genome_build)

  if (have_build) {
    dir_name <- tryCatch({
      resp <- httr::GET(paste0(.NCBI_GENOMES_URL, cfg$root, "/"), httr::timeout(30))
      listing <- httr::content(resp, as = "text", encoding = "UTF-8")
      pattern <- paste0("GCA_[0-9]+\\.[0-9]+_", stringr::str_escape(genome_build), "(?=/)")
      hits <- unique(stringr::str_extract_all(listing, pattern)[[1]])
      if (length(hits) > 0) hits[[1]] else NULL
    }, error = function(e) NULL)
  }

  if (is.null(dir_name)) {
    dir_name <- cfg$fallback
    if (have_build && !stringr::str_detect(dir_name, stringr::fixed(genome_build))) {
      cli::cli_alert_warning(
        "No NCBI assembly report for build {.val {genome_build}}; using {.val {dir_name}}."
      )
    }
  }

  paste0(.NCBI_GENOMES_URL, cfg$root, "/", dir_name, "/", dir_name, "_assembly_report.txt")
}

#' Parse an NCBI Assembly Report
#'
#' Assembly reports are tab-separated with a leading '#' comment block whose last
#' line carries the column names. Line endings are CRLF in the comment block, so
#' carriage returns are stripped before parsing.
#'
#' @param path Character string path to the downloaded report.
#'
#' @return A tibble with the ten \code{.ASSEMBLY_REPORT_COLS}, all character.
#'
#' @keywords internal
.parse_assembly_report <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- sub("\r$", "", lines)
  lines <- lines[!startsWith(lines, "#") & nzchar(lines)]

  if (length(lines) == 0) cli::cli_abort("Assembly report at {.file {path}} contains no records.")

  df <- utils::read.delim(
    text = paste(lines, collapse = "\n"),
    header = FALSE, sep = "\t", quote = "",
    colClasses = "character", stringsAsFactors = FALSE
  )

  if (ncol(df) < length(.ASSEMBLY_REPORT_COLS)) {
    cli::cli_abort(
      "Assembly report has {ncol(df)} column{?s}; expected {length(.ASSEMBLY_REPORT_COLS)}."
    )
  }

  df <- df[, seq_along(.ASSEMBLY_REPORT_COLS), drop = FALSE]
  colnames(df) <- .ASSEMBLY_REPORT_COLS
  as_tibble(df)
}

#' Coerce a Cached Assembly Report Back to Character
#'
#' \code{.manage_cache()} reads TSVs with \code{read.delim()}, which type-guesses
#' columns; sequence names such as "1" come back as integers. Every field of an
#' assembly report is an identifier, so force them all to character.
#'
#' @param df A dataframe read from cache or bundled data.
#'
#' @return A tibble with all columns as character.
#'
#' @keywords internal
.as_assembly_report <- function(df) {
  df <- as_tibble(df)
  missing_cols <- setdiff(.ASSEMBLY_REPORT_COLS, names(df))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Cached assembly report is missing column{?s}: {.field {missing_cols}}.")
  }
  df %>% mutate(across(everything(), as.character))
}

#' Get the NCBI Assembly Report for a Species
#'
#' Downloads the assembly report matching the Ensembl genome build, caching it
#' alongside the annotation data with the same offline fallback behaviour as the
#' other sources. The report is the reference that turns an Ensembl scaffold name
#' into a human-readable chromosome label such as \code{chr14(GL000194.1)}.
#'
#' @param species Character string. One of "human", "mouse", or "rat".
#' @param annotation_dir Character string of the cache directory.
#' @param genome_build Character string of the Ensembl genome build (e.g.
#'   "GRCh38.p14"). Used to pick the matching NCBI assembly and to invalidate a
#'   cache built against a different build.
#' @param force_update Logical.
#' @param offline Logical.
#'
#' @return A tibble of the assembly report.
#'
#' @keywords internal
.get_assembly_report <- function(species,
                                 annotation_dir = .get_data_path(),
                                 genome_build = NULL,
                                 force_update = FALSE,
                                 offline = getOption("geneRosetta.offline", FALSE)) {
  species <- tolower(species)
  folder <- file.path(annotation_dir, "Assembly_Reports")
  m_file <- file.path(folder, paste0(species, "_assembly_report.tsv"))
  l_file <- file.path(folder, paste0(species, "_assembly_report_log.json"))

  query_url <- .assembly_report_url(species, genome_build)
  remote_stamp <- if (isTRUE(offline)) NA else .get_remote_date(query_url)

  # A cache built against a different genome build describes different scaffolds.
  build_matches <- TRUE
  if (file.exists(l_file) && !is.null(genome_build) && !is.na(genome_build)) {
    cached_build <- tryCatch(jsonlite::read_json(l_file)$genome_build, error = function(e) NULL)
    build_matches <- is.null(cached_build) ||
      identical(as.character(cached_build), as.character(genome_build))
  }

  online <- !is.null(remote_stamp) && !is.na(remote_stamp)

  cache <- if (!build_matches && online && !isTRUE(offline)) {
    list(data = NULL, download = TRUE)
  } else {
    .manage_cache(m_file, l_file, remote_stamp, sub_dir = "Assembly_Reports", force_update = force_update, offline = offline)
  }

  if (!cache$download) return(.as_assembly_report(cache$data))

  if (isTRUE(offline)) {
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "Assembly_Reports")
    if (fb$success) return(.as_assembly_report(fb$data))
    cli::cli_abort("Offline mode enabled and no cached/bundled assembly report found for {species}.")
  }

  cli::cli_progress_step("Downloading NCBI assembly report for {species}...")

  tmp <- tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  
  report <- tryCatch({
    utils::download.file(query_url, destfile = tmp, mode = "wb", quiet = TRUE)
    parsed <- .parse_assembly_report(tmp)
    .save_annotation(parsed, m_file, l_file, list(
      data_last_modified = as.character(remote_stamp),
      genome_build = if (is.null(genome_build) || is.na(genome_build)) NA else as.character(genome_build),
      data_url = query_url
    ))
    parsed
  }, error = function(e) {
    if (file.exists(m_file)) {
      cli::cli_alert_warning("Assembly report download failed; using cached data from {basename(m_file)}")
      return(.read_annotation_tsv(m_file))
    }
    fb <- .materialize_bundled_fallback(m_file, l_file, sub_dir = "Assembly_Reports")
    if (fb$success) {
      cli::cli_alert_warning("Assembly report download failed; using bundled package data.")
      return(fb$data)
    }
    cli::cli_abort("Failed to download or load assembly report for {species}: {conditionMessage(e)}")
  })

  .as_assembly_report(report)
}

#' Build Scaffold Name Lookups From an Assembly Report
#'
#' Ensembl names scaffolds inconsistently: unlocalized and unplaced scaffolds are
#' named by their GenBank accession ("GL000194.1"), while alt-scaffolds and
#' patches are named by their sequence name ("HG1_PATCH", "HSCHR6_MHC_APD_CTG1").
#' Indexing every alias a report offers therefore resolves all of them.
#'
#' @param report A tibble from \code{.get_assembly_report()}.
#'
#' @return A list with two named character vectors:
#'   \item{label}{alias -> \code{chr<molecule>(<accession>)} label}
#'   \item{acc_role}{GenBank accession -> sequence role}
#'
#' @keywords internal
.chromosome_lookup <- function(report) {
  scaffolds <- report[report$sequence_role != "assembled-molecule", , drop = FALSE]

  molecule <- scaffolds$assigned_molecule
  molecule[is.na(molecule) | molecule %in% c("na", "")] <- "Un"
  label <- paste0("chr", molecule, "(", scaffolds$genbank_accn, ")")

  alias_cols <- c("sequence_name", "genbank_accn", "refseq_accn", "ucsc_style_name")
  aliases <- unlist(scaffolds[alias_cols], use.names = FALSE)
  labels <- rep(label, times = length(alias_cols))

  keep <- !is.na(aliases) & nzchar(aliases) & aliases != "na"
  aliases <- aliases[keep]
  labels <- labels[keep]

  first <- !duplicated(aliases)

  list(
    label = stats::setNames(labels[first], aliases[first]),
    acc_role = stats::setNames(scaffolds$sequence_role, scaffolds$genbank_accn)
  )
}

#' Build a Chromosome Lookup Without Failing
#'
#' Sources whose chromosome names are already plain (MGI, RGD) should not lose an
#' entire download because the assembly report is unreachable. Returns NULL on any
#' failure, which \code{.format_chromosome()} treats as "leave names untouched".
#'
#' @inheritParams .get_assembly_report
#'
#' @return A list from \code{.chromosome_lookup()}, or NULL.
#'
#' @keywords internal
.try_chromosome_lookup <- function(species,
                                   annotation_dir = .get_data_path(),
                                   genome_build = NULL,
                                   offline = getOption("geneRosetta.offline", FALSE)) {
  tryCatch(
    .chromosome_lookup(.get_assembly_report(species, annotation_dir, genome_build, offline = offline)),
    error = function(e) {
      cli::cli_alert_warning(
        "Could not load the {species} assembly report; leaving chromosome names as published."
      )
      NULL
    }
  )
}

#' Label Non-Standard Chromosomes
#'
#' Rewrites scaffold names to \code{chr<molecule>(<GenBank accession>)}, e.g.
#' \code{GL000194.1} becomes \code{chr14(GL000194.1)}. Scaffolds with no assigned
#' molecule become \code{chrUn(<accession>)}. Assembled molecules (1-22, X, Y, MT)
#' and names absent from the assembly report are returned unchanged.
#'
#' @param x Character vector of chromosome names as reported by the source.
#' @param lookup A list from \code{.chromosome_lookup()}, or NULL to pass \code{x}
#'   through unchanged.
#'
#' @return A character vector the same length as \code{x}.
#'
#' @keywords internal
.format_chromosome <- function(x, lookup) {
  x <- as.character(x)
  if (is.null(lookup)) return(x)
  out <- x
  idx <- which(!is.na(x) & nzchar(x) & !(x %in% .STANDARD_CHROMOSOMES))
  if (length(idx) > 0) {
    hit <- unname(lookup$label[x[idx]])
    out[idx] <- ifelse(is.na(hit), x[idx], hit)
  }
  out
}

#' Recover the Sequence Role of a Formatted Chromosome Label
#'
#' @param labels Character vector of chromosome labels as produced by
#'   \code{.format_chromosome()}.
#' @param lookup A list from \code{.chromosome_lookup()}, or NULL. With NULL only
#'   assembled molecules can be recognised; everything else is \code{"unknown"}.
#'
#' @return A character vector of sequence roles; \code{"unknown"} where the label
#'   is neither an assembled molecule nor a scaffold present in the report.
#'
#' @keywords internal
.chromosome_role <- function(labels, lookup) {
  labels <- as.character(labels)
  role <- rep(NA_character_, length(labels))
  role[!is.na(labels) & labels %in% .STANDARD_CHROMOSOMES] <- "assembled-molecule"

  accession <- stringr::str_match(labels, "^chr[^(]*\\(([^)]+)\\)$")[, 2]
  idx <- which(is.na(role) & !is.na(accession))
  if (length(idx) > 0 && !is.null(lookup)) {
    role[idx] <- unname(lookup$acc_role[accession[idx]])
  }

  role[is.na(role)] <- "unknown"
  role
}

#' Restrict an Annotation Table to a Chromosome Scope
#'
#' @param df A dataframe with a \code{chromosome} column of formatted labels.
#' @param scope One of "all", "primary", or "standard".
#' @param lookup A list from \code{.chromosome_lookup()}, or NULL to skip
#'   filtering entirely.
#'
#' @return \code{df}, filtered to the sequence roles the scope admits.
#'
#' @details Without an assembly report there is no way to tell an unplaced
#'   scaffold from an alt-haplotype, so a NULL \code{lookup} returns \code{df}
#'   untouched rather than dropping rows it cannot classify. This mirrors
#'   \code{.format_chromosome()}, which leaves names as published when the report
#'   is unavailable.
#'
#' @keywords internal
.filter_chromosome_scope <- function(df, scope, lookup) {
  keep_roles <- .CHROMOSOME_SCOPE_ROLES[[scope]]
  if (is.null(keep_roles)) return(df)
  if (is.null(lookup)) {
    cli::cli_alert_warning(
      "No assembly report available; keeping every sequence instead of applying
       {.val {scope}} scope."
    )
    return(df)
  }
  df[.chromosome_role(df$chromosome, lookup) %in% keep_roles, , drop = FALSE]
}
