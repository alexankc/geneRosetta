#' Build Merged Gene Annotation for Multiple Species
#'
#' A unified entry point to build or load gene annotations for Human, Mouse, or Rat.
#' This function automatically routes the request to the species-specific
#' builder while maintaining a consistent interface.
#'
#' @param species Character string. One of "human", "mouse", or "rat".
#' @param annotation_dir Character string specifying the directory where annotation
#'   data files are located. Defaults to the persistent user data directory.
#' @param force_update Logical. If TRUE, forces re-download and rebuild of annotations
#'   even if a valid cache exists. Defaults to FALSE.
#'
#' @return A dataframe containing the merged annotation for the selected species.
#' @export
#'
#' @examples
#' \dontrun{
#' # The simple way
#' human_annot <- buildAnnotation("human")
#'
#' # Forcing an update for mouse
#' mouse_annot <- buildAnnotation("mouse", force_update = TRUE)
#' }
buildAnnotation <- function(species = c("human", "mouse", "rat"),
                            annotation_dir = .get_data_path(),
                            force_update = FALSE,
                            offline = getOption("geneRosetta.offline", FALSE)) {
  species <- match.arg(species)

  result <- switch(species,
    "human" = buildHumanAnnotation(annotation_dir = annotation_dir, force_update = force_update, offline = offline),
    "mouse" = buildMouseAnnotation(annotation_dir = annotation_dir, force_update = force_update, offline = offline),
    "rat"   = buildRatAnnotation(annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  )

  return(result)
}


# Collect authority symbols and synonyms from Ensembl-less rows before those
# rows are folded into Ensembl-linked records by NCBI ID. rows_patch() only
# fills missing cells, so without this sidecar a non-primary authority symbol
# is lost whenever the retained row already has authority annotation.
.collect_authority_orphan_synonyms <- function(df, symbol_col, synonym_col) {
  required <- c("ncbi_id", symbol_col, synonym_col)
  if (!all(required %in% names(df))) {
    cli::cli_abort("Cannot collect orphan synonyms; missing columns: {toString(setdiff(required, names(df)))}")
  }

  payload <- purrr::map2_chr(df[[symbol_col]], df[[synonym_col]], function(symbol, synonyms) {
    values <- c(symbol, synonyms)
    values <- values[!is.na(values) & values != ""]
    if (length(values) == 0) "" else paste(values, collapse = "|")
  })

  data.frame(
    ncbi_id = df$ncbi_id,
    syn_authority_orphan = payload,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(ncbi_id), ncbi_id != "", syn_authority_orphan != "") %>%
    dplyr::group_by(ncbi_id) %>%
    dplyr::summarise(
      syn_authority_orphan = paste(unique(syn_authority_orphan), collapse = "|"),
      .groups = "drop"
    )
}


# Supply NCBI-owned metadata through the stable NCBI ID. This avoids treating
# an Ensembl cross-reference as proof that two distinct NCBI genes are one gene.
.backfill_ncbi_metadata_by_id <- function(df, df_ncbi) {
  metadata <- df_ncbi %>%
    dplyr::filter(!is.na(ncbi_id)) %>%
    dplyr::group_by(ncbi_id) %>%
    dplyr::summarise(
      gene_biotype_ncbi_exact = dplyr::first(
        stats::na.omit(gene_biotype_ncbi), default = NA_character_
      ),
      gene_description_ncbi_exact = dplyr::first(
        stats::na.omit(gene_description_ncbi), default = NA_character_
      ),
      chromosome_ncbi_exact = dplyr::first(
        stats::na.omit(chromosome_ncbi), default = NA_character_
      ),
      gene_start_ncbi_exact = dplyr::first(
        stats::na.omit(gene_start_ncbi), default = NA_integer_
      ),
      gene_end_ncbi_exact = dplyr::first(
        stats::na.omit(gene_end_ncbi), default = NA_integer_
      ),
      strand_ncbi_exact = dplyr::first(
        stats::na.omit(strand_ncbi), default = NA_integer_
      ),
      .groups = "drop"
    )

  df %>%
    dplyr::left_join(metadata, by = "ncbi_id", na_matches = "never") %>%
    dplyr::mutate(
      gene_biotype = dplyr::coalesce(
        gene_biotype, gene_biotype_ncbi_exact
      ),
      gene_description = dplyr::coalesce(
        gene_description, gene_description_ncbi_exact
      ),
      chromosome = dplyr::coalesce(chromosome, chromosome_ncbi_exact),
      gene_start = dplyr::coalesce(gene_start, gene_start_ncbi_exact),
      gene_end = dplyr::coalesce(gene_end, gene_end_ncbi_exact),
      strand = dplyr::coalesce(strand, strand_ncbi_exact),
      metadata_source = dplyr::coalesce(
        metadata_source,
        dplyr::if_else(
          !is.na(chromosome_ncbi_exact) | !is.na(gene_start_ncbi_exact) |
            !is.na(gene_end_ncbi_exact) | !is.na(strand_ncbi_exact),
          "NCBI",
          NA_character_
        )
      )
    ) %>%
    dplyr::select(
      -gene_biotype_ncbi_exact,
      -gene_description_ncbi_exact,
      -chromosome_ncbi_exact,
      -gene_start_ncbi_exact,
      -gene_end_ncbi_exact,
      -strand_ncbi_exact
    )
}


# Translate the broad NCBI type_of_gene vocabulary to the closest Ensembl gene
# biotype spelling used elsewhere in the merged annotation. More specific RNA
# classes already shared by both sources pass through unchanged.
.harmonize_gene_biotype <- function(x) {
  dplyr::case_when(
    is.na(x) ~ NA_character_,
    x == "unknown" ~ NA_character_,
    x == "protein-coding" ~ "protein_coding",
    x == "ncRNA" ~ "lncRNA",
    x == "pseudo" ~ "pseudogene",
    x == "biological-region" ~ "biological_region",
    x == "other" ~ "misc_RNA",
    TRUE ~ x
  )
}


# A deprecated or missing Ensembl link is not a reason to discard a distinct
# authority record. Drop a non-current row only when a complete current row has
# the same exact NCBI-authority ID pair; otherwise retain the record and clear
# any obsolete Ensembl ID.
.resolve_deprecated_records <- function(df, current_ensembl_ids,
                                        authority_id = "hgnc_id") {
  if (!authority_id %in% names(df)) {
    cli::cli_abort("Authority ID column {.field {authority_id}} is missing.")
  }

  is_current <- !is.na(df$ensembl_id) &
    df$ensembl_id %in% current_ensembl_ids
  current_complete <- df %>%
    dplyr::filter(
      is_current,
      !is.na(gene_biotype),
      !is.na(gene_description)
    )

  authority_values <- df[[authority_id]]
  current_authority_values <- current_complete[[authority_id]]
  stable_identity <- paste(df$ncbi_id, authority_values, sep = "\r")
  current_stable_identity <- paste(
    current_complete$ncbi_id, current_authority_values, sep = "\r"
  )
  has_stable_identity <- !is.na(df$ncbi_id) & !is.na(authority_values)
  current_has_stable_identity <-
    !is.na(current_complete$ncbi_id) & !is.na(current_authority_values)
  duplicate_identity <- has_stable_identity &
    stable_identity %in% current_stable_identity[current_has_stable_identity]
  is_non_current <- !is_current

  df %>%
    dplyr::filter(!(is_non_current & duplicate_identity)) %>%
    dplyr::mutate(
      ensembl_id = dplyr::if_else(
        !is.na(ensembl_id) & !ensembl_id %in% current_ensembl_ids,
        NA_character_,
        ensembl_id
      )
    )
}


#' Build Merged Human Gene Annotation
#'
#' Merges gene annotation data from multiple sources (BioMart/Ensembl, NCBI, and HGNC)
#' into a unified human gene annotation dataset. Implements multi-tier matching strategy
#' and intelligent symbol prioritization.
#'
#' @param annotation_dir Character string specifying the directory.
#'   Defaults to the persistent user data directory.
#' @param force_update Logical. If TRUE, forces re-download and rebuild of annotations
#'   even if a valid cache exists. Defaults to FALSE.
#'
#' @return A dataframe with the following columns:
#'   \item{ensembl_id}{Ensembl gene identifier (primary key)}
#'   \item{gene_symbol}{Official gene symbol (priority: HGNC > NCBI > BioMart)}
#'   \item{gene_biotype}{Gene biotype classification}
#'   \item{chromosome}{Chromosome location}
#'   \item{gene_start}{Genomic start position}
#'   \item{gene_end}{Genomic end position}
#'   \item{strand}{Strand orientation (+ or -)}
#'   \item{metadata_source}{Metadata source: Ensembl or NCBI}
#'   \item{gene_synonyms}{Pipe-separated merged gene synonyms from all sources}
#'   \item{ncbi_id}{NCBI Gene ID}
#'   \item{hgnc_id}{HGNC identifier}
#'   \item{gene_description}{Gene description}
#'
#' @details
#' The merging strategy uses the following tiers:
#' 1. BioMart -> HGNC (via Ensembl ID)
#' 2. BioMart -> NCBI (via Ensembl ID) for unmatched BioMart records
#' 3. Attempt to rescue NCBI IDs via symbol matching
#' 4. Attempt to rescue HGNC IDs via symbol matching
#'
#' Gene symbols are prioritized as: HGNC (official) > NCBI > BioMart.
#' Gene synonyms are deduplicated, sorted, and merged with pipes (|) as separators.
#' Results are cached locally as TSV with accompanying JSON metadata log.
#'
#' @examples
#' \dontrun{
#' human_annotation <- buildHumanAnnotation()
#' human_annotation <- buildHumanAnnotation(force_update = TRUE)
#' }
#'
#' @import dplyr
#' @export
#'
buildHumanAnnotation <- function(annotation_dir = .get_data_path(), force_update = FALSE, offline = getOption("geneRosetta.offline", FALSE)) {
  annotation_file <- file.path(annotation_dir, "human_merged_annotation.tsv")
  summary_log_file <- file.path(annotation_dir, "human_merged_annotation_log.json")

  # Check cache
  cache_complete <- file.exists(annotation_file) && file.exists(summary_log_file)

  if (cache_complete && !force_update) {
    cli::cli_alert_info("Loading cached annotation from: {.file {annotation_file}}")
    combined_final <- readr::read_tsv(
      annotation_file,
      # The builders emit an all-character table and use "" for "no symbol", so
      # the cache must be read back the same way: na = "NA" keeps "" intact, and
      # the explicit col_types stops readr guessing ncbi_id/gene_start numeric.
      # Without this, buildAnnotation() returned different types and values on a
      # cache hit than it did on the initial build.
      na = "NA",
      col_types = readr::cols(.default = readr::col_character())
    )
    return(combined_final)
  }

  if (!cache_complete && !force_update) {
    cli::cli_alert_warning("No cache found (missing data or log). Creating merged annotation...")
  }

  message("Fetching annotations...")

  human_biomart_annot <- getAnnotationBiomart("human", annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  human_ncbi_annot <- getAnnotationNCBI("human", annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  human_hgnc_annot <- getAnnotationHGNChuman(annotation_dir = annotation_dir, force_update = force_update, offline = offline)

  # Standardize column types and rename to prevent collisions
  clean_na_str <- function(df) {
    df %>% mutate(across(where(is.character), ~ if_else(is.na(.x) | .x == "" | tolower(.x) %in% c("na", "null", "nan"), NA_character_, .x)))
  }

  df_bm <- human_biomart_annot %>%
    clean_na_str() %>%
    mutate(across(c(ensembl_id), as.character)) %>%
    mutate(across(any_of(c("gene_start", "gene_end", "strand")), as.integer)) %>%
    rename(symbol_bm = gene_symbol, syn_bm = gene_synonyms) %>%
    rename(any_of(c(
      chromosome_bm = "chromosome",
      gene_start_bm = "gene_start",
      gene_end_bm = "gene_end",
      strand_bm = "strand"
    )))

  df_ncbi <- human_ncbi_annot %>%
    clean_na_str() %>%
    mutate(across(c(ncbi_id, ensembl_id), as.character)) %>%
    mutate(across(any_of(c("gene_start", "gene_end", "strand")), as.integer)) %>%
    rename(any_of(c(
      symbol_ncbi = "gene_symbol",
      syn_ncbi = "gene_synonyms",
      gene_biotype_ncbi = "gene_biotype",
      gene_description_ncbi = "gene_description",
      chromosome_ncbi = "chromosome",
      gene_start_ncbi = "gene_start",
      gene_end_ncbi = "gene_end",
      strand_ncbi = "strand"
    )))

  if (!"symbol_ncbi" %in% colnames(df_ncbi)) df_ncbi$symbol_ncbi <- NA_character_
  if (!"syn_ncbi" %in% colnames(df_ncbi)) df_ncbi$syn_ncbi <- NA_character_
  if (!"gene_biotype_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_biotype_ncbi <- NA_character_
  if (!"gene_description_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_description_ncbi <- NA_character_

  if (!"chromosome_bm" %in% colnames(df_bm)) df_bm$chromosome_bm <- NA_character_
  if (!"gene_start_bm" %in% colnames(df_bm)) df_bm$gene_start_bm <- NA_integer_
  if (!"gene_end_bm" %in% colnames(df_bm)) df_bm$gene_end_bm <- NA_integer_
  if (!"strand_bm" %in% colnames(df_bm)) df_bm$strand_bm <- NA_integer_

  if (!"chromosome_ncbi" %in% colnames(df_ncbi)) df_ncbi$chromosome_ncbi <- NA_character_
  if (!"gene_start_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_start_ncbi <- NA_integer_
  if (!"gene_end_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_end_ncbi <- NA_integer_
  if (!"strand_ncbi" %in% colnames(df_ncbi)) df_ncbi$strand_ncbi <- NA_integer_

  df_hgnc <- human_hgnc_annot %>%
    clean_na_str() %>%
    mutate(across(c(hgnc_id, ensembl_id, ncbi_id), as.character)) %>%
    rename(
      symbol_hgnc = gene_symbol,
      syn_hgnc = gene_synonyms
    )

  message("Merging annotations...")

  # Merge database records using Ensembl IDs
  bm_hgnc_match <- inner_join(df_bm %>% filter(!is.na(ensembl_id)), df_hgnc %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_unmatched <- anti_join(df_bm, df_hgnc %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_ncbi_match <- inner_join(bm_unmatched %>% filter(!is.na(ensembl_id)), df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_orphans <- anti_join(bm_unmatched, df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")

  hgnc_orphans <- df_hgnc %>% filter(!ensembl_id %in% df_bm$ensembl_id | is.na(ensembl_id))
  hgnc_ncbi_match <- inner_join(hgnc_orphans %>% filter(!is.na(ensembl_id)), df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  hgnc_orphans_unmatched <- anti_join(hgnc_orphans, df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")

  ncbi_orphans <- df_ncbi %>%
    filter(!ensembl_id %in% df_bm$ensembl_id | is.na(ensembl_id)) %>%
    filter(!ensembl_id %in% df_hgnc$ensembl_id | is.na(ensembl_id))

  combined_raw <- bind_rows(
    bm_hgnc_match,
    bm_ncbi_match,
    bm_orphans,
    hgnc_ncbi_match,
    hgnc_orphans_unmatched,
    ncbi_orphans
  )

  # Prioritize BioMart coordinates, fallback to NCBI coordinates if BioMart is missing
  combined_raw <- combined_raw %>%
    mutate(
      chromosome = coalesce(chromosome_bm, chromosome_ncbi),
      gene_start = coalesce(gene_start_bm, gene_start_ncbi),
      gene_end = coalesce(gene_end_bm, gene_end_ncbi),
      strand = coalesce(strand_bm, strand_ncbi),
      metadata_source = case_when(
        !is.na(chromosome_bm) | !is.na(gene_start_bm) |
          !is.na(gene_end_bm) | !is.na(strand_bm) ~ "Ensembl",
        !is.na(chromosome_ncbi) | !is.na(gene_start_ncbi) |
          !is.na(gene_end_ncbi) | !is.na(strand_ncbi) ~ "NCBI",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-any_of(c("chromosome_bm", "chromosome_ncbi",
                            "gene_start_bm", "gene_start_ncbi",
                            "gene_end_bm", "gene_end_ncbi",
                            "strand_bm", "strand_ncbi")))

  # Merge Ensembl-unlinked records into existing matches via NCBI ID
  ens_present <- combined_raw %>% filter(!is.na(ensembl_id))
  ens_missing_ncbi_present <- combined_raw %>% filter(is.na(ensembl_id) & !is.na(ncbi_id))
  ens_missing_ncbi_missing <- combined_raw %>% filter(is.na(ensembl_id) & is.na(ncbi_id))

  authority_orphan_synonyms <- .collect_authority_orphan_synonyms(
    ens_missing_ncbi_present, "symbol_hgnc", "syn_hgnc"
  )

  ens_missing_ncbi_present_collapsed <- ens_missing_ncbi_present %>%
    group_by(ncbi_id) %>%
    summarise(across(everything(), ~ {
      if (all(is.na(.x))) .x[1] else na.omit(.x)[1]
    }), .groups = "drop")

  overlap_ncbi <- intersect(ens_present$ncbi_id, ens_missing_ncbi_present_collapsed$ncbi_id)
  overlap_ncbi <- overlap_ncbi[!is.na(overlap_ncbi)]

  if (length(overlap_ncbi) > 0) {
    to_merge_into <- ens_present %>% filter(ncbi_id %in% overlap_ncbi)
    to_merge_from <- ens_missing_ncbi_present_collapsed %>% filter(ncbi_id %in% overlap_ncbi)

    merged <- to_merge_into %>%
      rows_patch(to_merge_from, by = "ncbi_id", unmatched = "ignore") %>%
      left_join(authority_orphan_synonyms, by = "ncbi_id")

    ens_present <- bind_rows(
      ens_present %>% filter(!ncbi_id %in% overlap_ncbi | is.na(ncbi_id)),
      merged
    )

    ens_missing_ncbi_present_collapsed <- ens_missing_ncbi_present_collapsed %>%
      filter(!ncbi_id %in% overlap_ncbi)
  }

  combined_raw <- bind_rows(
    ens_present,
    ens_missing_ncbi_present_collapsed,
    ens_missing_ncbi_missing
  )

  # Retain conflicting NCBI symbols and synonyms for Ensembl-linked records
  ncbi_all_syns <- df_ncbi %>%
    filter(!is.na(ensembl_id)) %>%
    tidyr::unite("all_syns", symbol_ncbi, syn_ncbi, sep = "|", remove = FALSE, na.rm = TRUE) %>%
    group_by(ensembl_id) %>%
    summarise(syn_ncbi_extra = paste(unique(all_syns[all_syns != ""]), collapse = "|"), .groups = "drop")

  combined_raw <- combined_raw %>%
    left_join(ncbi_all_syns, by = "ensembl_id")

  # BioMart metadata applies to BioMart rows. For NCBI-derived rows, fill the
  # same fields directly from NCBI gene_info (type_of_gene and description).
  # Define primary gene symbol (Priority: HGNC -> NCBI -> BioMart)
  combined_symbol <- combined_raw %>%
    mutate(gene_symbol = coalesce(symbol_hgnc, symbol_ncbi, symbol_bm))

  # Match additional NCBI IDs by symbol for records lacking NCBI ID
  combined_symbol_ncbi_is_NA <- combined_symbol %>%
    filter(is.na(ncbi_id)) %>%
    dplyr::select(-ncbi_id)
  combined_symbol_1 <- combined_symbol %>%
    filter(!is.na(ncbi_id))
  df_ncbi_noEnsembl <- df_ncbi %>%
    dplyr::select(symbol_ncbi, syn_ncbi, ncbi_id)
  combined_symbol_extraNCBI <- combined_symbol_ncbi_is_NA %>%
    left_join(df_ncbi_noEnsembl, join_by("gene_symbol" == "symbol_ncbi")) %>%
    tidyr::unite("syn_ncbi", syn_ncbi.x, syn_ncbi.y, sep = "|", na.rm = TRUE)
  combined_symbol_2 <- bind_rows(combined_symbol_1, combined_symbol_extraNCBI)

  # Match additional HGNC IDs by symbol for records lacking HGNC ID
  combined_symbol_2_hgnc_is_NA <- combined_symbol_2 %>%
    filter(is.na(hgnc_id)) %>%
    dplyr::select(-hgnc_id)
  combined_symbol_3 <- combined_symbol_2 %>%
    filter(!is.na(hgnc_id))
  df_hgnc_noEnsembl <- df_hgnc %>%
    dplyr::select(symbol_hgnc, syn_hgnc, hgnc_id)
  combined_symbol_extraHGNC <- combined_symbol_2_hgnc_is_NA %>%
    left_join(df_hgnc_noEnsembl, join_by("gene_symbol" == "symbol_hgnc")) %>%
    tidyr::unite("syn_hgnc", syn_hgnc.x, syn_hgnc.y, sep = "|", na.rm = TRUE)
  combined_symbol_3 <- bind_rows(combined_symbol_3, combined_symbol_extraHGNC)

  # Normalize missing values
  combined_symbol_4 <- combined_symbol_3 %>%
    mutate(across(everything(), ~ ifelse(.x %in% c("null", "", "NA"), NA_character_, .x))) %>%
    mutate(gene_symbol = if_else(!is.na(ensembl_id) & gene_symbol == ensembl_id | is.na(gene_symbol), "", gene_symbol))

  # Concatenate, clean, and deduplicate gene synonyms across sources
  combined_final <- combined_symbol_4 %>%
    tidyr::unite("temp_synonyms", any_of(c("symbol_bm", "symbol_ncbi", "symbol_hgnc", "syn_bm", "syn_ncbi", "syn_hgnc", "syn_ncbi_extra", "syn_authority_orphan")), sep = "|", na.rm = TRUE) %>%
    mutate(gene_synonyms = purrr::map2_chr(stringi::stri_split_fixed(temp_synonyms, "|"), gene_symbol, function(x, sym) {
      x <- stringi::stri_replace_all_fixed(x, '"', "")
      x <- stringi::stri_trim_both(x)
      x <- x[x != "" & !is.na(x) & x != "NA" & x != "-" & x != sym]

      # Standardize open reading frame (orf) naming to uppercase
      orf_indices <- stringi::stri_detect_fixed(x, "orf", case_insensitive = TRUE)
      x[orf_indices] <- stringi::stri_trans_toupper(x[orf_indices])

      unique_syns <- stringi::stri_sort(unique(x))

      if (length(unique_syns) == 0) {
        return(NA_character_)
      }
      stringi::stri_join(unique_syns, collapse = "|")
    })) %>%
    dplyr::select(-temp_synonyms) %>%
    dplyr::select(
      ensembl_id, gene_symbol, gene_biotype, chromosome,
      gene_start, gene_end, strand, metadata_source, gene_synonyms,
      ncbi_id, hgnc_id, gene_description
    ) %>%
    .backfill_ncbi_metadata_by_id(df_ncbi) %>%
    .resolve_deprecated_records(df_bm$ensembl_id) %>%
    dplyr::filter(!(is.na(gene_biotype) & is.na(gene_description))) %>%
    dplyr::filter(!(is.na(ensembl_id) & is.na(ncbi_id))) %>%
    dplyr::mutate(gene_biotype = .harmonize_gene_biotype(gene_biotype)) %>%
    dplyr::distinct()

  if (!dir.exists(annotation_dir)) dir.create(annotation_dir, recursive = TRUE)
  readr::write_tsv(combined_final, annotation_file)

  # Save metadata log
  source_logs <- list(
    Ensembl = file.path(annotation_dir, "Ensembl_Genes", "human_biomart_mapping_log.json"),
    NCBI    = file.path(annotation_dir, "NCBI_Genes", "human_ncbi_mapping_log.json"),
    HGNC    = file.path(annotation_dir, "HGNC_human_Genes", "human_hgnc_mapping_log.json")
  )

  log_contents <- purrr::map(source_logs, function(path) {
    if (file.exists(path)) {
      return(jsonlite::read_json(path))
    } else {
      return("Log file not found")
    }
  })

  summary_log <- list(
    project_name = "geneRosetta human merged annotation",
    gene_count = nrow(combined_final),
    source_details = log_contents
  )

  jsonlite::write_json(summary_log, summary_log_file, auto_unbox = TRUE, pretty = TRUE)

  cli::cli_rule()
  cli::cli_alert_success("Build Complete")
  cli::cli_rule()

  return(combined_final)
}


#' Build Merged Mouse Gene Annotation
#'
#' Merges gene annotation data from multiple sources (BioMart/Ensembl, MGI, and NCBI)
#' into a unified mouse gene annotation dataset. Implements multi-tier matching with
#' preference for MGI nomenclature as the official mouse gene naming authority.
#'
#' @param annotation_dir Character string specifying the directory.
#'   Defaults to the persistent user data directory.
#' @param force_update Logical. If TRUE, forces re-download and rebuild of annotations
#'   even if a valid cache exists. Defaults to FALSE.
#'
#' @return A dataframe with the following columns:
#'   \item{ensembl_id}{Ensembl gene identifier (primary key)}
#'   \item{gene_symbol}{Official gene symbol (prioritizes MGI nomenclature)}
#'   \item{gene_biotype}{Gene biotype classification from BioMart}
#'   \item{chromosome}{Chromosome location}
#'   \item{gene_start}{Genomic start position from BioMart/Ensembl}
#'   \item{gene_end}{Genomic end position from BioMart/Ensembl}
#'   \item{strand}{Strand orientation (+ or -)}
#'   \item{metadata_source}{Metadata source: Ensembl or NCBI}
#'   \item{gene_synonyms}{Pipe-separated merged gene synonyms from all sources}
#'   \item{ncbi_id}{NCBI Gene ID}
#'   \item{mgi_id}{MGI accession ID}
#'   \item{gene_description}{Gene description}
#'
#' @details
#' The merging strategy uses multiple tiers:
#' 1. BioMart -> MGI (perfect match via Ensembl ID AND Symbol)
#' 2. BioMart -> MGI (Ensembl ID match with updated symbol from MGI)
#' 3. BioMart -> NCBI (rescue for MGI orphans via Ensembl ID)
#' 4. Rescue missing MGI IDs via NCBI ID matching
#' 5. Rescue missing MGI IDs via symbol matching
#'
#' MGI provides the official mouse nomenclature and is preferred for gene symbols
#' when available. Genomic coordinates are sourced from BioMart/Ensembl for consistency.
#' Results are cached locally as TSV with accompanying JSON metadata log.
#'
#' @examples
#' \dontrun{
#' mouse_annotation <- buildMouseAnnotation()
#' mouse_annotation <- buildMouseAnnotation(force_update = TRUE)
#' }
#'
#' @import dplyr
#' @export
#'
buildMouseAnnotation <- function(annotation_dir = .get_data_path(), force_update = FALSE, offline = getOption("geneRosetta.offline", FALSE)) {
  annotation_file <- file.path(annotation_dir, "mouse_merged_annotation.tsv")
  summary_log_file <- file.path(annotation_dir, "mouse_merged_annotation_log.json")

  # Check cache
  cache_complete <- file.exists(annotation_file) && file.exists(summary_log_file)

  if (cache_complete && !force_update) {
    cli::cli_alert_info("Loading cached mouse annotation from: {.file {annotation_file}}")
    combined_final <- readr::read_tsv(
      annotation_file,
      # The builders emit an all-character table and use "" for "no symbol", so
      # the cache must be read back the same way: na = "NA" keeps "" intact, and
      # the explicit col_types stops readr guessing ncbi_id/gene_start numeric.
      # Without this, buildAnnotation() returned different types and values on a
      # cache hit than it did on the initial build.
      na = "NA",
      col_types = readr::cols(.default = readr::col_character())
    )
    return(combined_final)
  }

  if (!cache_complete && !force_update) {
    cli::cli_alert_warning("No cache found (missing data or log). Creating merged annotation...")
  }

  cli::cli_alert_info("Fetching annotations...")

  mouse_mgi_annot <- getAnnotationMGImouse(annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  mouse_biomart_annot <- getAnnotationBiomart("mouse", annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  mouse_ncbi_annot <- getAnnotationNCBI("mouse", annotation_dir = annotation_dir, force_update = force_update, offline = offline)

  # Standardize column types and rename to prevent collisions
  clean_na_str <- function(df) {
    df %>% mutate(across(where(is.character), ~ if_else(is.na(.x) | .x == "" | tolower(.x) %in% c("na", "null", "nan"), NA_character_, .x)))
  }

  # Standardize IDs and rename synonyms to prevent collision during merge
  df_bm <- mouse_biomart_annot %>%
    clean_na_str() %>%
    mutate(across(c(ensembl_id), as.character)) %>%
    mutate(across(any_of(c("gene_start", "gene_end", "strand")), as.integer)) %>%
    rename(symbol_bm = gene_symbol, syn_bm = gene_synonyms) %>%
    rename(any_of(c(
      chromosome_bm = "chromosome",
      gene_start_bm = "gene_start",
      gene_end_bm = "gene_end",
      strand_bm = "strand"
    )))

  # MGI: Exclude duplicate genomic columns to keep BioMart/Ensembl as coordinate authority
  ensembl_dup_cols <- c(
    "gene_biotype", "chromosome", "gene_start", "gene_end",
    "strand", "gene_description"
  )

  df_mgi <- mouse_mgi_annot %>%
    clean_na_str() %>%
    dplyr::select(-any_of(ensembl_dup_cols)) %>%
    mutate(across(c(mgi_id, ensembl_id, ncbi_id), as.character)) %>%
    rename(symbol_mgi = gene_symbol, syn_mgi = gene_synonyms)

  df_ncbi <- mouse_ncbi_annot %>%
    clean_na_str() %>%
    mutate(across(c(ncbi_id, ensembl_id), as.character)) %>%
    mutate(across(any_of(c("gene_start", "gene_end", "strand")), as.integer)) %>%
    rename(any_of(c(
      symbol_ncbi = "gene_symbol",
      syn_ncbi = "gene_synonyms",
      gene_biotype_ncbi = "gene_biotype",
      gene_description_ncbi = "gene_description",
      chromosome_ncbi = "chromosome",
      gene_start_ncbi = "gene_start",
      gene_end_ncbi = "gene_end",
      strand_ncbi = "strand"
    )))

  if (!"symbol_ncbi" %in% colnames(df_ncbi)) df_ncbi$symbol_ncbi <- NA_character_
  if (!"syn_ncbi" %in% colnames(df_ncbi)) df_ncbi$syn_ncbi <- NA_character_
  if (!"gene_biotype_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_biotype_ncbi <- NA_character_
  if (!"gene_description_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_description_ncbi <- NA_character_

  if (!"chromosome_bm" %in% colnames(df_bm)) df_bm$chromosome_bm <- NA_character_
  if (!"gene_start_bm" %in% colnames(df_bm)) df_bm$gene_start_bm <- NA_integer_
  if (!"gene_end_bm" %in% colnames(df_bm)) df_bm$gene_end_bm <- NA_integer_
  if (!"strand_bm" %in% colnames(df_bm)) df_bm$strand_bm <- NA_integer_

  if (!"chromosome_ncbi" %in% colnames(df_ncbi)) df_ncbi$chromosome_ncbi <- NA_character_
  if (!"gene_start_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_start_ncbi <- NA_integer_
  if (!"gene_end_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_end_ncbi <- NA_integer_
  if (!"strand_ncbi" %in% colnames(df_ncbi)) df_ncbi$strand_ncbi <- NA_integer_

  cli::cli_alert_info("Merging mouse annotations...")

  # Merge database records using Ensembl IDs
  bm_mgi_match <- inner_join(df_bm %>% filter(!is.na(ensembl_id)), df_mgi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_unmatched <- anti_join(df_bm, df_mgi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_ncbi_match <- inner_join(bm_unmatched %>% filter(!is.na(ensembl_id)), df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_orphans <- anti_join(bm_unmatched, df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")

  mgi_orphans <- df_mgi %>% filter(!ensembl_id %in% df_bm$ensembl_id | is.na(ensembl_id))
  mgi_ncbi_match <- inner_join(mgi_orphans %>% filter(!is.na(ensembl_id)), df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  mgi_orphans_unmatched <- anti_join(mgi_orphans, df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")

  ncbi_orphans <- df_ncbi %>%
    filter(!ensembl_id %in% df_bm$ensembl_id | is.na(ensembl_id)) %>%
    filter(!ensembl_id %in% df_mgi$ensembl_id | is.na(ensembl_id))

  combined_raw <- bind_rows(
    bm_mgi_match,
    bm_ncbi_match,
    bm_orphans,
    mgi_ncbi_match,
    mgi_orphans_unmatched,
    ncbi_orphans
  )

  # Prioritize BioMart coordinates, fallback to NCBI coordinates if BioMart is missing
  combined_raw <- combined_raw %>%
    mutate(
      chromosome = coalesce(chromosome_bm, chromosome_ncbi),
      gene_start = coalesce(gene_start_bm, gene_start_ncbi),
      gene_end = coalesce(gene_end_bm, gene_end_ncbi),
      strand = coalesce(strand_bm, strand_ncbi),
      metadata_source = case_when(
        !is.na(chromosome_bm) | !is.na(gene_start_bm) |
          !is.na(gene_end_bm) | !is.na(strand_bm) ~ "Ensembl",
        !is.na(chromosome_ncbi) | !is.na(gene_start_ncbi) |
          !is.na(gene_end_ncbi) | !is.na(strand_ncbi) ~ "NCBI",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-any_of(c("chromosome_bm", "chromosome_ncbi",
                            "gene_start_bm", "gene_start_ncbi",
                            "gene_end_bm", "gene_end_ncbi",
                            "strand_bm", "strand_ncbi")))

  # Merge Ensembl-unlinked records into existing matches via NCBI ID
  ens_present <- combined_raw %>% filter(!is.na(ensembl_id))
  ens_missing_ncbi_present <- combined_raw %>% filter(is.na(ensembl_id) & !is.na(ncbi_id))
  ens_missing_ncbi_missing <- combined_raw %>% filter(is.na(ensembl_id) & is.na(ncbi_id))

  authority_orphan_synonyms <- .collect_authority_orphan_synonyms(
    ens_missing_ncbi_present, "symbol_mgi", "syn_mgi"
  )

  ens_missing_ncbi_present_collapsed <- ens_missing_ncbi_present %>%
    group_by(ncbi_id) %>%
    summarise(across(everything(), ~ {
      if (all(is.na(.x))) .x[1] else na.omit(.x)[1]
    }), .groups = "drop")

  overlap_ncbi <- intersect(ens_present$ncbi_id, ens_missing_ncbi_present_collapsed$ncbi_id)
  overlap_ncbi <- overlap_ncbi[!is.na(overlap_ncbi)]

  if (length(overlap_ncbi) > 0) {
    to_merge_into <- ens_present %>% filter(ncbi_id %in% overlap_ncbi)
    to_merge_from <- ens_missing_ncbi_present_collapsed %>% filter(ncbi_id %in% overlap_ncbi)

    merged <- to_merge_into %>%
      rows_patch(to_merge_from, by = "ncbi_id", unmatched = "ignore") %>%
      left_join(authority_orphan_synonyms, by = "ncbi_id")

    ens_present <- bind_rows(
      ens_present %>% filter(!ncbi_id %in% overlap_ncbi | is.na(ncbi_id)),
      merged
    )

    ens_missing_ncbi_present_collapsed <- ens_missing_ncbi_present_collapsed %>%
      filter(!ncbi_id %in% overlap_ncbi)
  }

  combined_raw <- bind_rows(
    ens_present,
    ens_missing_ncbi_present_collapsed,
    ens_missing_ncbi_missing
  )

  # Retain conflicting NCBI symbols and synonyms for Ensembl-linked records
  ncbi_all_syns <- df_ncbi %>%
    filter(!is.na(ensembl_id)) %>%
    tidyr::unite("all_syns", symbol_ncbi, syn_ncbi, sep = "|", remove = FALSE, na.rm = TRUE) %>%
    group_by(ensembl_id) %>%
    summarise(syn_ncbi_extra = paste(unique(all_syns[all_syns != ""]), collapse = "|"), .groups = "drop")

  combined_raw <- combined_raw %>%
    left_join(ncbi_all_syns, by = "ensembl_id") %>%
    mutate(
      gene_biotype = coalesce(gene_biotype, gene_biotype_ncbi),
      gene_description = coalesce(gene_description, gene_description_ncbi)
    )

  # Define primary gene symbol (Priority: MGI -> NCBI -> BioMart)
  combined_symbol <- combined_raw %>%
    mutate(gene_symbol = coalesce(symbol_mgi, symbol_ncbi, symbol_bm))

  # Match additional NCBI IDs by symbol for records lacking NCBI ID
  combined_symbol_ncbi_is_NA <- combined_symbol %>%
    filter(is.na(ncbi_id)) %>%
    dplyr::select(-ncbi_id)
  combined_symbol_1 <- combined_symbol %>%
    filter(!is.na(ncbi_id))
  df_ncbi_noEnsembl <- df_ncbi %>%
    dplyr::select(symbol_ncbi, syn_ncbi, ncbi_id)
  combined_symbol_extraNCBI <- combined_symbol_ncbi_is_NA %>%
    left_join(df_ncbi_noEnsembl, join_by("gene_symbol" == "symbol_ncbi")) %>%
    tidyr::unite("syn_ncbi", syn_ncbi.x, syn_ncbi.y, sep = "|", na.rm = TRUE)
  combined_symbol_2 <- bind_rows(combined_symbol_1, combined_symbol_extraNCBI)

  # Match additional MGI IDs by symbol for records lacking MGI ID
  combined_symbol_2_mgi_is_NA <- combined_symbol_2 %>%
    filter(is.na(mgi_id)) %>%
    dplyr::select(-mgi_id)
  combined_symbol_3 <- combined_symbol_2 %>%
    filter(!is.na(mgi_id))
  df_mgi_noEnsembl <- df_mgi %>%
    dplyr::select(symbol_mgi, syn_mgi, mgi_id)
  combined_symbol_extraMGI <- combined_symbol_2_mgi_is_NA %>%
    left_join(df_mgi_noEnsembl, join_by("gene_symbol" == "symbol_mgi")) %>%
    tidyr::unite("syn_mgi", syn_mgi.x, syn_mgi.y, sep = "|", na.rm = TRUE)
  combined_raw <- bind_rows(combined_symbol_3, combined_symbol_extraMGI)

  # Normalize missing values
  combined_raw <- combined_raw %>%
    mutate(across(everything(), ~ ifelse(.x %in% c("null", "", "NA"), NA_character_, .x))) %>%
    mutate(gene_symbol = if_else(!is.na(ensembl_id) & gene_symbol == ensembl_id | is.na(gene_symbol), "", gene_symbol))

  # Concatenate, clean, and deduplicate gene synonyms across sources
  cli::cli_alert_info("Cleaning synonyms and finalizing...")

  combined_final <- combined_raw %>%
    tidyr::unite("temp_synonyms", any_of(c("symbol_bm", "symbol_ncbi", "symbol_mgi", "syn_bm", "syn_ncbi", "syn_mgi", "syn_ncbi_extra", "syn_authority_orphan")), sep = "|", na.rm = TRUE) %>%
    mutate(gene_synonyms = purrr::map2_chr(stringi::stri_split_fixed(temp_synonyms, "|"), gene_symbol, function(x, sym) {
      x <- stringi::stri_trim_both(x)
      x <- x[x != "" & !is.na(x) & x != "NA" & x != "-" & x != sym]

      unique_syns <- stringi::stri_sort(unique(x))

      if (length(unique_syns) == 0) {
        return(NA_character_)
      }
      stringi::stri_join(unique_syns, collapse = "|")
    })) %>%
    dplyr::select(
      any_of(c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome",
      "gene_start", "gene_end", "strand", "metadata_source", "gene_synonyms",
      "ncbi_id", "mgi_id", "gene_description"))
    ) %>%
    .backfill_ncbi_metadata_by_id(df_ncbi) %>%
    .resolve_deprecated_records(
      df_bm$ensembl_id,
      authority_id = "mgi_id"
    ) %>%
    dplyr::filter(!(is.na(ensembl_id) & is.na(ncbi_id))) %>%
    dplyr::mutate(gene_biotype = .harmonize_gene_biotype(gene_biotype)) %>%
    dplyr::distinct()

  if (!dir.exists(annotation_dir)) dir.create(annotation_dir, recursive = TRUE)
  readr::write_tsv(combined_final, annotation_file)

  # Save metadata log
  source_logs <- list(
    Ensembl = file.path(annotation_dir, "Ensembl_Genes", "mouse_biomart_mapping_log.json"),
    NCBI = file.path(annotation_dir, "NCBI_Genes", "mouse_ncbi_mapping_log.json"),
    MGI = file.path(annotation_dir, "MGI_mouse_genes", "mouse_mgi_mapping_log.json")
  )

  log_contents <- purrr::map(source_logs, function(path) {
    if (file.exists(path)) {
      return(jsonlite::read_json(path))
    } else {
      return("Log file not found")
    }
  })

  summary_log <- list(
    project_name = "geneRosetta mouse merged annotation",
    gene_count = nrow(combined_final),
    source_details = log_contents
  )

  jsonlite::write_json(summary_log, summary_log_file, auto_unbox = TRUE, pretty = TRUE)

  cli::cli_rule()
  cli::cli_alert_success("Build Complete")
  cli::cli_rule()

  return(combined_final)
}

# =======================================

#' Build Merged Rat Gene Annotation
#'
#' Merges gene annotation data from multiple sources (BioMart/Ensembl, RGD, and NCBI)
#' into a unified rat gene annotation dataset. Implements multi-tier matching with
#' RGD as the preferred source for rat nomenclature and identity validation.
#'
#' @param annotation_dir Character string specifying the directory.
#'   Defaults to the persistent user data directory.
#' @param force_update Logical. If TRUE, forces re-download and rebuild of annotations
#'   even if a valid cache exists. Defaults to FALSE.
#'
#' @return A dataframe with the following columns:
#'   \item{ensembl_id}{Ensembl gene identifier (primary key)}
#'   \item{gene_symbol}{Official gene symbol (prioritizes RGD nomenclature)}
#'   \item{gene_biotype}{Gene biotype classification from BioMart}
#'   \item{chromosome}{Chromosome location}
#'   \item{gene_start}{Genomic start position from BioMart/Ensembl}
#'   \item{gene_end}{Genomic end position from BioMart/Ensembl}
#'   \item{strand}{Strand orientation (+ or -)}
#'   \item{metadata_source}{Metadata source: Ensembl or NCBI}
#'   \item{gene_synonyms}{Pipe-separated merged gene synonyms from all sources}
#'   \item{ncbi_id}{NCBI Gene ID (validated against NCBI reference)}
#'   \item{rgd_id}{RGD identifier}
#'   \item{gene_description}{Gene description}
#'
#' @details
#' The merging strategy uses multiple tiers:
#' 1. BioMart -> RGD (perfect match via Ensembl ID AND Symbol)
#' 2. BioMart -> RGD (Ensembl ID match with updated symbol from RGD)
#' 3. BioMart -> NCBI (rescue for RGD orphans via Ensembl ID)
#' 4. Validate multiple NCBI IDs against the NCBI reference database
#'
#' RGD provides the official rat nomenclature and is preferred for gene symbols
#' when available. Genomic coordinates are sourced from BioMart/Ensembl for consistency.
#' NCBI IDs are validated and cleaned when multiple values are present (separated
#' by semicolons in the source data). Results are cached locally as TSV with
#' accompanying JSON metadata log.
#'
#' @examples
#' \dontrun{
#' rat_annotation <- buildRatAnnotation()
#' rat_annotation <- buildRatAnnotation(force_update = TRUE)
#' }
#'
#' @import dplyr
#' @export
#'
buildRatAnnotation <- function(annotation_dir = .get_data_path(), force_update = FALSE, offline = getOption("geneRosetta.offline", FALSE)) {
  annotation_file <- file.path(annotation_dir, "rat_merged_annotation.tsv")
  summary_log_file <- file.path(annotation_dir, "rat_merged_annotation_log.json")

  # Check cache
  cache_complete <- file.exists(annotation_file) && file.exists(summary_log_file)

  if (cache_complete && !force_update) {
    cli::cli_alert_info("Loading cached rat annotation from: {.file {annotation_file}}")
    combined_final <- readr::read_tsv(
      annotation_file,
      # The builders emit an all-character table and use "" for "no symbol", so
      # the cache must be read back the same way: na = "NA" keeps "" intact, and
      # the explicit col_types stops readr guessing ncbi_id/gene_start numeric.
      # Without this, buildAnnotation() returned different types and values on a
      # cache hit than it did on the initial build.
      na = "NA",
      col_types = readr::cols(.default = readr::col_character())
    )
    return(combined_final)
  }

  if (!cache_complete && !force_update) {
    cli::cli_alert_warning("No cache found. Creating merged rat annotation...")
  }

  cli::cli_alert_info("Fetching rat annotations...")

  rat_rgd_annot <- getAnnotationRGDrat(annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  rat_biomart_annot <- getAnnotationBiomart("rat", annotation_dir = annotation_dir, force_update = force_update, offline = offline)
  rat_ncbi_annot <- getAnnotationNCBI("rat", annotation_dir = annotation_dir, force_update = force_update, offline = offline)

  # Standardize column types and rename to prevent collisions
  clean_na_str <- function(df) {
    df %>% mutate(across(where(is.character), ~ if_else(is.na(.x) | .x == "" | tolower(.x) %in% c("na", "null", "nan"), NA_character_, .x)))
  }

  # Standardize IDs and rename synonyms to prevent collision
  df_bm <- rat_biomart_annot %>%
    clean_na_str() %>%
    mutate(across(c(ensembl_id), as.character)) %>%
    mutate(across(any_of(c("gene_start", "gene_end", "strand")), as.integer)) %>%
    rename(symbol_bm = gene_symbol, syn_bm = gene_synonyms) %>%
    rename(any_of(c(
      chromosome_bm = "chromosome",
      gene_start_bm = "gene_start",
      gene_end_bm = "gene_end",
      strand_bm = "strand"
    )))

  # RGD: Exclude duplicate genomic columns to keep BioMart/Ensembl as coordinate authority
  ensembl_dup_cols <- c(
    "gene_biotype", "chromosome", "gene_start", "gene_end",
    "strand", "gene_description"
  )

  df_rgd <- rat_rgd_annot %>%
    clean_na_str() %>%
    dplyr::select(-any_of(ensembl_dup_cols)) %>%
    mutate(across(c(rgd_id, ensembl_id, ncbi_id), as.character)) %>%
    rename(symbol_rgd = gene_symbol, syn_rgd = gene_synonyms)

  df_ncbi <- rat_ncbi_annot %>%
    clean_na_str() %>%
    mutate(across(c(ncbi_id, ensembl_id), as.character)) %>%
    mutate(across(any_of(c("gene_start", "gene_end", "strand")), as.integer)) %>%
    rename(any_of(c(
      symbol_ncbi = "gene_symbol",
      syn_ncbi = "gene_synonyms",
      gene_biotype_ncbi = "gene_biotype",
      gene_description_ncbi = "gene_description",
      chromosome_ncbi = "chromosome",
      gene_start_ncbi = "gene_start",
      gene_end_ncbi = "gene_end",
      strand_ncbi = "strand"
    )))

  if (!"symbol_ncbi" %in% colnames(df_ncbi)) df_ncbi$symbol_ncbi <- NA_character_
  if (!"syn_ncbi" %in% colnames(df_ncbi)) df_ncbi$syn_ncbi <- NA_character_
  if (!"gene_biotype_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_biotype_ncbi <- NA_character_
  if (!"gene_description_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_description_ncbi <- NA_character_

  if (!"chromosome_bm" %in% colnames(df_bm)) df_bm$chromosome_bm <- NA_character_
  if (!"gene_start_bm" %in% colnames(df_bm)) df_bm$gene_start_bm <- NA_integer_
  if (!"gene_end_bm" %in% colnames(df_bm)) df_bm$gene_end_bm <- NA_integer_
  if (!"strand_bm" %in% colnames(df_bm)) df_bm$strand_bm <- NA_integer_

  if (!"chromosome_ncbi" %in% colnames(df_ncbi)) df_ncbi$chromosome_ncbi <- NA_character_
  if (!"gene_start_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_start_ncbi <- NA_integer_
  if (!"gene_end_ncbi" %in% colnames(df_ncbi)) df_ncbi$gene_end_ncbi <- NA_integer_
  if (!"strand_ncbi" %in% colnames(df_ncbi)) df_ncbi$strand_ncbi <- NA_integer_

  cli::cli_alert_info("Merging rat annotations...")

  # Merge database records using Ensembl IDs
  bm_rgd_match <- inner_join(df_bm %>% filter(!is.na(ensembl_id)), df_rgd %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_unmatched <- anti_join(df_bm, df_rgd %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_ncbi_match <- inner_join(bm_unmatched %>% filter(!is.na(ensembl_id)), df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  bm_orphans <- anti_join(bm_unmatched, df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")

  rgd_orphans <- df_rgd %>% filter(!ensembl_id %in% df_bm$ensembl_id | is.na(ensembl_id))
  rgd_ncbi_match <- inner_join(rgd_orphans %>% filter(!is.na(ensembl_id)), df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")
  rgd_orphans_unmatched <- anti_join(rgd_orphans, df_ncbi %>% filter(!is.na(ensembl_id)), by = "ensembl_id", na_matches = "never")

  ncbi_orphans <- df_ncbi %>%
    filter(!ensembl_id %in% df_bm$ensembl_id | is.na(ensembl_id)) %>%
    filter(!ensembl_id %in% df_rgd$ensembl_id | is.na(ensembl_id))

  combined_raw <- bind_rows(
    bm_rgd_match,
    bm_ncbi_match,
    bm_orphans,
    rgd_ncbi_match,
    rgd_orphans_unmatched,
    ncbi_orphans
  )

  # Prioritize BioMart coordinates, fallback to NCBI coordinates if BioMart is missing
  combined_raw <- combined_raw %>%
    mutate(
      chromosome = coalesce(chromosome_bm, chromosome_ncbi),
      gene_start = coalesce(gene_start_bm, gene_start_ncbi),
      gene_end = coalesce(gene_end_bm, gene_end_ncbi),
      strand = coalesce(strand_bm, strand_ncbi),
      metadata_source = case_when(
        !is.na(chromosome_bm) | !is.na(gene_start_bm) |
          !is.na(gene_end_bm) | !is.na(strand_bm) ~ "Ensembl",
        !is.na(chromosome_ncbi) | !is.na(gene_start_ncbi) |
          !is.na(gene_end_ncbi) | !is.na(strand_ncbi) ~ "NCBI",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::select(-any_of(c("chromosome_bm", "chromosome_ncbi",
                            "gene_start_bm", "gene_start_ncbi",
                            "gene_end_bm", "gene_end_ncbi",
                            "strand_bm", "strand_ncbi")))

  # Validate semicolon-separated NCBI IDs from RGD against NCBI reference
  valid_ncbi_ids <- unique(na.omit(df_ncbi$ncbi_id))

  combined_raw <- combined_raw %>%
    mutate(
      ncbi_id = purrr::map_chr(ncbi_id, function(x) {
        if (is.na(x) || !stringr::str_detect(x, ";")) {
          return(x)
        }

        ids <- stringi::stri_split_fixed(x, ";")[[1]] %>%
          stringi::stri_trim_both()

        valid_ids <- ids[ids %in% valid_ncbi_ids]

        if (length(valid_ids) == 0) {
          NA_character_
        } else {
          stringi::stri_join(unique(valid_ids), collapse = "|")
        }
      })
    )

  # Merge Ensembl-unlinked records into existing matches via NCBI ID
  ens_present <- combined_raw %>% filter(!is.na(ensembl_id))
  ens_missing_ncbi_present <- combined_raw %>% filter(is.na(ensembl_id) & !is.na(ncbi_id))
  ens_missing_ncbi_missing <- combined_raw %>% filter(is.na(ensembl_id) & is.na(ncbi_id))

  authority_orphan_synonyms <- .collect_authority_orphan_synonyms(
    ens_missing_ncbi_present, "symbol_rgd", "syn_rgd"
  )

  ens_missing_ncbi_present_collapsed <- ens_missing_ncbi_present %>%
    group_by(ncbi_id) %>%
    summarise(across(everything(), ~ {
      if (all(is.na(.x))) .x[1] else na.omit(.x)[1]
    }), .groups = "drop")

  overlap_ncbi <- intersect(ens_present$ncbi_id, ens_missing_ncbi_present_collapsed$ncbi_id)
  overlap_ncbi <- overlap_ncbi[!is.na(overlap_ncbi)]

  if (length(overlap_ncbi) > 0) {
    to_merge_into <- ens_present %>% filter(ncbi_id %in% overlap_ncbi)
    to_merge_from <- ens_missing_ncbi_present_collapsed %>% filter(ncbi_id %in% overlap_ncbi)

    merged <- to_merge_into %>%
      rows_patch(to_merge_from, by = "ncbi_id", unmatched = "ignore") %>%
      left_join(authority_orphan_synonyms, by = "ncbi_id")

    ens_present <- bind_rows(
      ens_present %>% filter(!ncbi_id %in% overlap_ncbi | is.na(ncbi_id)),
      merged
    )

    ens_missing_ncbi_present_collapsed <- ens_missing_ncbi_present_collapsed %>%
      filter(!ncbi_id %in% overlap_ncbi)
  }

  combined_raw <- bind_rows(
    ens_present,
    ens_missing_ncbi_present_collapsed,
    ens_missing_ncbi_missing
  )

  # Retain conflicting NCBI symbols and synonyms for Ensembl-linked records
  ncbi_all_syns <- df_ncbi %>%
    filter(!is.na(ensembl_id)) %>%
    tidyr::unite("all_syns", symbol_ncbi, syn_ncbi, sep = "|", remove = FALSE, na.rm = TRUE) %>%
    group_by(ensembl_id) %>%
    summarise(syn_ncbi_extra = paste(unique(all_syns[all_syns != ""]), collapse = "|"), .groups = "drop")

  combined_raw <- combined_raw %>%
    left_join(ncbi_all_syns, by = "ensembl_id") %>%
    mutate(
      gene_biotype = coalesce(gene_biotype, gene_biotype_ncbi),
      gene_description = coalesce(gene_description, gene_description_ncbi)
    )

  # Define primary gene symbol (Priority: RGD -> NCBI -> BioMart)
  combined_symbol <- combined_raw %>%
    mutate(gene_symbol = coalesce(symbol_rgd, symbol_ncbi, symbol_bm))

  # Match additional NCBI IDs by symbol for records lacking NCBI ID
  combined_symbol_ncbi_is_NA <- combined_symbol %>%
    filter(is.na(ncbi_id)) %>%
    dplyr::select(-ncbi_id)
  combined_symbol_1 <- combined_symbol %>%
    filter(!is.na(ncbi_id))
  df_ncbi_noEnsembl <- df_ncbi %>%
    dplyr::select(symbol_ncbi, syn_ncbi, ncbi_id)
  combined_symbol_extraNCBI <- combined_symbol_ncbi_is_NA %>%
    left_join(df_ncbi_noEnsembl, join_by("gene_symbol" == "symbol_ncbi")) %>%
    tidyr::unite("syn_ncbi", syn_ncbi.x, syn_ncbi.y, sep = "|", na.rm = TRUE)
  combined_symbol_2 <- bind_rows(combined_symbol_1, combined_symbol_extraNCBI)

  # Match additional RGD IDs by symbol for records lacking RGD ID
  combined_symbol_2_rgd_is_NA <- combined_symbol_2 %>%
    filter(is.na(rgd_id)) %>%
    dplyr::select(-rgd_id)
  combined_symbol_3 <- combined_symbol_2 %>%
    filter(!is.na(rgd_id))
  df_rgd_noEnsembl <- df_rgd %>%
    dplyr::select(symbol_rgd, syn_rgd, rgd_id)
  combined_symbol_extraRGD <- combined_symbol_2_rgd_is_NA %>%
    left_join(df_rgd_noEnsembl, join_by("gene_symbol" == "symbol_rgd")) %>%
    tidyr::unite("syn_rgd", syn_rgd.x, syn_rgd.y, sep = "|", na.rm = TRUE)
  combined_raw <- bind_rows(combined_symbol_3, combined_symbol_extraRGD)

  # Normalize missing values
  combined_raw <- combined_raw %>%
    mutate(across(everything(), ~ ifelse(.x %in% c("null", "", "NA"), NA_character_, .x))) %>%
    mutate(gene_symbol = if_else(!is.na(ensembl_id) & gene_symbol == ensembl_id | is.na(gene_symbol), "", gene_symbol))

  # Concatenate, clean, and deduplicate gene synonyms across sources
  cli::cli_alert_info("Cleaning synonyms and finalizing...")

  combined_final <- combined_raw %>%
    tidyr::unite("temp_synonyms", any_of(c("symbol_bm", "symbol_ncbi", "symbol_rgd", "syn_bm", "syn_ncbi", "syn_rgd", "syn_ncbi_extra", "syn_authority_orphan")), sep = "|", na.rm = TRUE) %>%
    mutate(gene_synonyms = purrr::map2_chr(stringi::stri_split_fixed(temp_synonyms, "|"), gene_symbol, function(x, sym) {
      x <- stringi::stri_trim_both(x)
      x <- x[x != "" & !is.na(x) & x != "NA" & x != "-" & x != sym]
      unique_syns <- stringi::stri_sort(unique(x))
      if (length(unique_syns) == 0) {
        return(NA_character_)
      }
      stringi::stri_join(unique_syns, collapse = "|")
    })) %>%
    dplyr::select(
      any_of(c("ensembl_id", "gene_symbol", "gene_biotype", "chromosome",
      "gene_start", "gene_end", "strand", "metadata_source", "gene_synonyms",
      "ncbi_id", "rgd_id", "gene_description"))
    ) %>%
    .backfill_ncbi_metadata_by_id(df_ncbi) %>%
    .resolve_deprecated_records(
      df_bm$ensembl_id,
      authority_id = "rgd_id"
    ) %>%
    dplyr::filter(!(is.na(ensembl_id) & is.na(ncbi_id))) %>%
    dplyr::mutate(gene_biotype = .harmonize_gene_biotype(gene_biotype)) %>%
    dplyr::distinct()

  if (!dir.exists(annotation_dir)) dir.create(annotation_dir, recursive = TRUE)
  readr::write_tsv(combined_final, annotation_file)

  # Save metadata log
  source_logs <- list(
    Ensembl = file.path(annotation_dir, "Ensembl_Genes", "rat_biomart_mapping_log.json"),
    NCBI    = file.path(annotation_dir, "NCBI_Genes", "rat_ncbi_mapping_log.json"),
    RGD     = file.path(annotation_dir, "RGD_rat_genes", "rat_rgd_mapping_log.json")
  )

  log_contents <- purrr::map(source_logs, function(path) {
    if (file.exists(path)) {
      return(jsonlite::read_json(path))
    } else {
      return("Log file not found")
    }
  })

  summary_log <- list(
    project_name = "geneRosetta rat merged annotation",
    gene_count = nrow(combined_final),
    source_details = log_contents
  )

  jsonlite::write_json(summary_log, summary_log_file, auto_unbox = TRUE, pretty = TRUE)

  cli::cli_rule()
  cli::cli_alert_success("Build Complete")
  cli::cli_rule()

  return(combined_final)
}


#' Summarize and Export Annotation Build Logs
#'
#' @description
#' Provides a human-readable summary of the gene annotation build process for
#' a specific species. It can print a detailed report to the console, export
#' the data as a formatted text file, generate manuscript-ready Methods text,
#' or copy the raw JSON log.
#'
#' @param species Character. One of \code{"human"}, \code{"mouse"}, or \code{"rat"}.
#'   Specifies which log file to retrieve. Defaults to \code{"human"}.
#' @param ensembl_annot Logical. If \code{TRUE}, summarises the Ensembl BioMart
#'   source log for the species instead of the merged annotation log.
#'   Defaults to \code{FALSE}.
#' @param export_format Character. The format to export the log:
#'   \itemize{
#'     \item \code{"none"}: Just print the summary to the R console (default).
#'     \item \code{"txt"}: Create a formatted, human-readable text file.
#'     \item \code{"methods"}: Create a prose summary suitable for a scientific
#'       Methods section.
#'     \item \code{"json"}: Copy the raw JSON log file to a local destination.
#'   }
#' @param export_path Character (optional). The file path where the log should be saved.
#'   If the directory does not exist, it will be created recursively.
#' @param annotation_dir Character string specifying the directory holding the
#'   annotation logs. Defaults to the persistent user data directory, which is
#'   where \code{buildAnnotation()} writes unless it was given a custom directory.
#'
#' @return The log data as a nested list, returned invisibly.
#' @export
#' @examples
#' \dontrun{
#' # Print human log to console
#' getAnnotationLogSummary(species = "human")
#'
#' # Export mouse log to a specific folder
#' getAnnotationLogSummary(
#'   species = "mouse",
#'   export_format = "methods",
#'   export_path = "logs/mouse_annotation_methods.txt"
#' )
#' }
#'
getAnnotationLogSummary <- function(
  species = c("human", "mouse", "rat"),
  ensembl_annot = FALSE,
  export_format = c("none", "txt", "methods", "json"),
  export_path = NULL,
  annotation_dir = .get_data_path()
) {
  species <- match.arg(species)
  export_format <- match.arg(export_format)

  # 1. Path setup
  if (ensembl_annot) {
    log_filename <- sprintf("%s_biomart_mapping_log.json", species)
    source_log <- file.path(annotation_dir, "Ensembl_Genes", log_filename)
  } else {
    log_filename <- sprintf("%s_merged_annotation_log.json", species)
    source_log <- file.path(annotation_dir, log_filename)
  }


  if (!file.exists(source_log)) {
    cli::cli_abort("No log file found for {.val {species}} at {.file {source_log}}. Run buildAnnotation(species={.val {species}}) first")
  }

  # Load data using jsonlite::
  log_data <- jsonlite::fromJSON(source_log, simplifyVector = FALSE)

  # If it's the Ensembl flat file, restructure it to match the 'Merged' format
  if (is.null(log_data$source_details)) {
    standardized <- list(
      project_name = paste("Ensembl BioMart Mapping -", stringr::str_to_title(species)),
      gene_count = NULL,
      source_details = list("BioMart" = log_data)
    )
    log_data <- standardized
  }

  # --- Helper: Ensure Directory Exists Recursively ---
  .ensure_dir <- function(path) {
    if (!is.null(path)) {
      dir_path <- dirname(path)
      if (!dir.exists(dir_path)) {
        cli::cli_alert_info("Creating directory: {.file {dir_path}}")
        dir.create(dir_path, recursive = TRUE)
      }
    }
  }

  # --- Case 1: JSON Export ---
  if (export_format == "json") {
    if (is.null(export_path)) export_path <- sprintf("%s_annotation_log.json", species)
    .ensure_dir(export_path)
    file.copy(from = source_log, to = export_path, overwrite = TRUE)
    cli::cli_alert_success("Log (json) exported to {.file {export_path}}")
    return(invisible(log_data))
  }

  # --- Case 2: TXT Export ---
  if (export_format == "txt") {
    if (is.null(export_path)) export_path <- sprintf("%s_annotation_log_readable.txt", species)
    .ensure_dir(export_path)

    txt <- c(
      paste("Project:", log_data$project_name),
      if (!is.null(log_data$gene_count)) paste("Gene count:", log_data$gene_count) else NULL,
      "",
      "Source Details:",
      ""
    )

    for (s_name in names(log_data$source_details)) {
      det <- log_data$source_details[[s_name]]
      txt <- c(txt, paste0(s_name, ":"))
      for (key in names(det)) {
        label <- stringr::str_to_title(gsub("_", " ", key))
        txt <- c(txt, paste0("    ", label, ": ", det[[key]]))
      }
      txt <- c(txt, "")
    }

    writeLines(txt, export_path)
    cli::cli_alert_success("Log (txt) exported to {.file {export_path}}")
    return(invisible(log_data))
  }

  # --- Case 3: Manuscript-ready Methods text ---
  if (export_format == "methods") {
    if (is.null(export_path)) export_path <- sprintf("%s_annotation_methods.txt", species)
    .ensure_dir(export_path)

    format_access_date <- function(details) {
      timestamp <- details$data_access_timestamp
      if (is.null(timestamp)) return(NULL)
      date <- as.Date(sub("T.*$", "", timestamp))
      if (is.na(date)) return(sub("T.*$", "", timestamp))
      format(date, "%B %Y")
    }

    package_description <- tryCatch(
      suppressWarnings(utils::packageDescription("geneRosetta")),
      error = function(e) NULL
    )
    if ((!is.list(package_description) || is.null(package_description$Version)) &&
        file.exists("DESCRIPTION")) {
      package_description <- as.list(read.dcf("DESCRIPTION")[1, ])
    }
    package_version <- if (is.null(package_description$Version)) {
      "unknown version"
    } else {
      paste0("v", package_description$Version)
    }
    package_url <- if (is.null(package_description$URL)) {
      "https://github.com/alexankc/geneRosetta"
    } else {
      sub(",.*$", "", package_description$URL)
    }

    details <- log_data$source_details
    ensembl <- details$Ensembl
    ensembl_release <- if (is.null(ensembl$ensembl_version)) {
      "Ensembl"
    } else {
      paste("Ensembl release", sub(".*_", "", ensembl$ensembl_version))
    }
    ensembl_build <- if (is.null(ensembl$genome_build)) "" else {
      paste0(" (", ensembl$genome_build, ")")
    }

    authority <- switch(species, human = "HGNC", mouse = "MGI", rat = "RGD")
    authority_name <- switch(
      species,
      human = "the HGNC complete gene set",
      mouse = "the MGI gene set",
      rat = "the RGD gene set"
    )
    ncbi_date <- format_access_date(details$NCBI)
    authority_date <- format_access_date(details[[authority]])
    ncbi_text <- paste0(
      "the NCBI Gene database",
      if (is.null(ncbi_date)) "" else paste0(" (accessed ", ncbi_date, ")")
    )
    authority_text <- paste0(
      authority_name,
      if (is.null(authority_date)) "" else paste0(" (accessed ", authority_date, ")")
    )

    authority_symbol <- switch(
      species,
      human = "HGNC-approved",
      mouse = "MGI-approved",
      rat = "RGD-approved"
    )
    count <- if (is.null(log_data$gene_count)) "an unspecified number of" else {
      format(log_data$gene_count, big.mark = ",", scientific = FALSE)
    }

    methods_text <- paste0(
      stringr::str_to_title(species), " gene annotations were generated using ",
      "geneRosetta (", package_version, "; ", package_url, ") by integrating ",
      ensembl_release, ensembl_build, ", ", ncbi_text, ", and ", authority_text, ". ",
      "Records were first matched by Ensembl gene IDs, with unmatched records ",
      "cross-referenced against NCBI Gene IDs; ", authority_symbol,
      " gene symbols were preferred as primary symbols, while conflicting source symbols ",
      "were retained as synonyms. Coding genes, non-coding RNAs, and pseudogenes were ",
      "retained; records lacking both an Ensembl and an NCBI identifier were excluded, ",
      "resulting in a final merged reference set ",
      "containing ", count, " unique entries."
    )

    writeLines(methods_text, export_path)
    cli::cli_alert_success("Methods text exported to {.file {export_path}}")
    return(invisible(log_data))
  }

  # --- Case 4: Console Output (Full Information) ---
  cli::cli_h1(log_data$project_name)
  if (!is.null(log_data$gene_count)) {
    cli::cli_alert_success("Total unique genes: {.val {log_data$gene_count}}")
    cli::cli_text(" ")
  }

  for (source_name in names(log_data$source_details)) {
    details <- log_data$source_details[[source_name]]

    cli::cli_div(theme = list(span.strong = list(color = "cyan")))
    cli::cli_alert_info("{.strong {source_name}}")
    cli::cli_ul()

    for (key in names(details)) {
      label <- stringr::str_to_title(gsub("_", " ", key))
      val <- details[[key]]

      if (grepl("url", key, ignore.case = TRUE)) {
        cli::cli_li("{label}: {.url {val}}")
      } else if (grepl("version|build|pkg|dataset", key, ignore.case = TRUE)) {
        cli::cli_li("{label}: {.val {val}}")
      } else if (grepl("timestamp|modified", key, ignore.case = TRUE)) {
        cli::cli_li("{label}: {.field {val}}")
      } else {
        cli::cli_li("{label}: {val}")
      }
    }
    cli::cli_end()
    cli::cli_end()
    cli::cli_text(" ")
  }

  cli::cli_rule()
  return(invisible(log_data))
}


# =======================================
