#' Map gene symbols to Ensembl IDs
#'
#' Maps input gene symbols to Ensembl IDs using a provided merged annotation.
#' Supports gene synonym rescue, configurable handling of multi-mapped symbols, and
#' optional retention of unmapped rows.
#'
#' @param input_data Character vector or data frame of gene symbols to map.
#' @param col_id If input_data is a data frame, the column containing the gene symbols.
#' @param annotation The output of buildHumanAnnotation(), buildMouseAnnotation(), buildRatAnnotation() or getAnnotationBiomart()).
#' @param add_metadata Logical; if FALSE, return only query, ensembl_id, gene_symbol.
#' @param multi_handling Strategy for multi-mapped symbols: "keepAll", "collapse", or "unique".
#' @param keep_unmapped Logical; if TRUE, retain rows that could not be mapped. Default is FALSE.
#' @return A list with elements:
#'   \item{mapped}{Tibble of mapped (and optionally unmapped) rows.}
#'   \item{unmapped_primary}{Character vector of symbols unmapped after primary pass.}
#'   \item{unmapped_final}{Character vector of symbols still unmapped after synonym rescue.}
#'   \item{multi_mapped}{Tibble summarizing symbols with multiple matches.}
#' @export

mapSymbolToID <- function(
  input_data,
  col_id = NULL,
  annotation,
  add_metadata = TRUE,
  multi_handling = c("keepAll", "collapse", "unique"),
  keep_unmapped = FALSE
) {
  multi_handling <- match.arg(multi_handling)

  # Validate required annotation columns
  required_annot_cols <- c("gene_symbol", "ensembl_id")
  if (!all(required_annot_cols %in% names(annotation))) {
    missing_cols <- setdiff(required_annot_cols, names(annotation))
    cli::cli_abort("The {.arg annotation} table is missing required column{?s}: {.val {missing_cols}}.")
  }

  # --- 1. Handle Vector vs Data Frame Input ---
  if (is.vector(input_data) || is.factor(input_data)) {
    original_df <- data.frame(query_gene_symbol = as.character(input_data), stringsAsFactors = FALSE) %>%
      mutate(row_id_internal = row_number())
  } else if (is.data.frame(input_data)) {
    if (is.null(col_id)) cli::cli_abort("Argument {.arg col_id} must be provided when input is a data frame.")
    if (!col_id %in% names(input_data)) cli::cli_abort("Column {.val {col_id}} not found in input data frame.")
    original_df <- input_data %>%
      dplyr::rename(query_gene_symbol = !!sym(col_id)) %>%
      mutate(row_id_internal = row_number())
  } else {
    cli::cli_abort("Input must be either a vector or a data frame.")
  }

  # --- 2. Handle Name Collisions ---
  overlap_cols <- intersect(names(original_df), names(annotation))
  overlap_cols <- setdiff(overlap_cols, "row_id_internal")

  if (length(overlap_cols) > 0) {
    cli::cli_alert_info("Prefixing overlapping input columns with {.val query_}: {.field {overlap_cols}}")
    original_df <- original_df %>%
      rename_with(~ paste0("query_", .), all_of(overlap_cols))
  }

  # --- 3. Standardize Annotation ---
  annotation <- annotation %>%
    mutate(
      across(everything(), as.character),
      annotation_row_internal = dplyr::row_number()
    )

  # --- 4. Mapping Pass 1: Primary Symbols ---
  to_map <- original_df %>% filter(!is.na(query_gene_symbol) & query_gene_symbol != "" & query_gene_symbol != "null")
  no_id <- original_df %>% filter(is.na(query_gene_symbol) | query_gene_symbol == "" | query_gene_symbol == "null")

  mapped_main_raw <- to_map %>%
    mutate(join_symbol_tmp = tolower(trimws(query_gene_symbol))) %>%
    left_join(
      annotation %>% filter(!is.na(gene_symbol)) %>% mutate(join_symbol_tmp = tolower(trimws(gene_symbol))),
      by = "join_symbol_tmp",
      relationship = "many-to-many"
    ) %>%
    dplyr::select(-join_symbol_tmp)

  mapped_main <- mapped_main_raw %>% filter(!is.na(annotation_row_internal))
  unmapped_main <- mapped_main_raw %>% filter(is.na(annotation_row_internal))
  unmapped_genes_primary <- unique(unmapped_main$query_gene_symbol)

  # --- 5. Mapping Pass 2: Synonym Mapping ---
  synonym_matches <- mapped_main[0, ]

  if ("gene_synonyms" %in% names(annotation)) {
    cli::cli_alert_info("Checking annotation synonyms...")

    synonym_db <- annotation %>%
      filter(!is.na(gene_synonyms) & gene_synonyms != "") %>%
      tidyr::separate_rows(gene_synonyms, sep = "\\|")

    synonym_matches <- to_map %>%
      mutate(query_gene_symbol_clean = tolower(trimws(query_gene_symbol))) %>%
      left_join(
        synonym_db %>% mutate(gene_synonyms_clean = tolower(trimws(gene_synonyms))),
        by = c("query_gene_symbol_clean" = "gene_synonyms_clean"),
        relationship = "many-to-many"
      ) %>%
      dplyr::select(-query_gene_symbol_clean) %>%
      filter(!is.na(annotation_row_internal))
  }

  all_matches <- bind_rows(mapped_main, synonym_matches) %>%
    distinct(row_id_internal, annotation_row_internal, .keep_all = TRUE)
  matched_row_ids <- unique(all_matches$row_id_internal)
  unmapped_genes_final <- to_map %>%
    filter(!row_id_internal %in% matched_row_ids) %>%
    pull(query_gene_symbol) %>%
    unique()
  unmapped_final_rows <- unmapped_main %>%
    filter(!row_id_internal %in% matched_row_ids)

  final_df <- if (keep_unmapped) {
    bind_rows(all_matches, unmapped_final_rows, no_id)
  } else {
    all_matches
  }

  # --- 6. Detect and Handle Multiple Ensembl IDs ---
  required_cols <- c("ensembl_id", "gene_symbol")
  for (col in required_cols) {
    if (!col %in% names(final_df)) {
      final_df[[col]] <- rep(NA_character_, nrow(final_df))
    }
  }

  multi_mapped_report <- final_df %>%
    filter(!is.na(query_gene_symbol)) %>%
    group_by(row_id_internal) %>%
    filter(n() > 1) %>%
    summarise(
      query_gene_symbol = first(query_gene_symbol),
      ensembl_ids = paste(unique(na.omit(ensembl_id)), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::select(-row_id_internal)

  if (nrow(multi_mapped_report) > 0) {
    cli::cli_alert_warning("{nrow(multi_mapped_report)} symbols produced multiple matches.")
    cli::cli_alert_info("Handling strategy: {.val {multi_handling}}")

    if (multi_handling == "collapse") {
      final_df <- final_df %>%
        group_by(row_id_internal) %>%
        summarise(
          across(everything(), ~ {
            vals <- unique(na.omit(.x))
            if (length(vals) == 0) NA_character_ else paste(vals, collapse = "|")
          }),
          .groups = "drop"
        )
    } else if (multi_handling == "unique") {
      final_df <- final_df %>%
        group_by(row_id_internal) %>%
        arrange(is.na(ensembl_id), ensembl_id) %>%
        slice(1) %>%
        ungroup()
    }
  }

  # --- 7. Cleanup and Output Formatting ---
  final_df <- final_df %>%
    arrange(row_id_internal) %>%
    dplyr::select(-any_of(c("gene_synonyms", "row_id_internal", "annotation_row_internal"))) %>%
    mutate(gene_symbol = if_else(query_gene_symbol %in% unmapped_genes_final,
      NA_character_,
      coalesce(gene_symbol, query_gene_symbol)
    )) %>%
    relocate(any_of(c("query_gene_symbol", "gene_symbol", "ensembl_id")))

  if (!add_metadata) {
    input_cols <- names(original_df)[names(original_df) != "row_id_internal"]
    final_df <- final_df %>% dplyr::select(query_gene_symbol, gene_symbol, ensembl_id, any_of(input_cols))
  }

  # --- 8. Summary Statistics ---
  successfully_mapped_symbols <- final_df %>%
    filter(!is.na(gene_symbol)) %>%
    pull(query_gene_symbol) %>%
    unique()
  total_input_count <- length(unique(original_df$query_gene_symbol))
  mapped_input_count <- length(successfully_mapped_symbols)
  num_primary_empty <- length(unique(no_id$query_gene_symbol))
  total_unmapped <- length(unmapped_genes_final)

  perc_mapped <- if (total_input_count > 0) {
    round((mapped_input_count / total_input_count) * 100, 1)
  } else {
    0
  }

  rescued_count <- length(setdiff(successfully_mapped_symbols, mapped_main$query_gene_symbol))

  cli::cli_h1("Mapping Summary (Symbol to ID)")
  cli::cli_alert_info("Total Unique Input Genes: {.val {total_input_count}}")

  if (rescued_count > 0) {
    cli::cli_alert_success(paste0(
      "Successfully Mapped: {.val {mapped_input_count}} ({perc_mapped}%) ",
      "- of which {.val {rescued_count}} {?was/were} mapped via synonyms"
    ))
  } else {
    cli::cli_alert_success("Successfully Mapped: {.val {mapped_input_count}} ({perc_mapped}%)")
  }

  if (total_unmapped > 0 || num_primary_empty > 0) {
    if (total_unmapped > 0) {
      cli::cli_alert_danger(paste0(
        "Final Unmapped: {.val {total_unmapped}} gene{?s} ",
        "could not be mapped."
      ))
    }

    if (num_primary_empty > 0) {
      cli::cli_alert_warning(paste0(
        "{.val {num_primary_empty}} input gene{?s} ",
        "{?was/were} NA, empty, or 'null' and {?was/were} skipped."
      ))
    }
  }

  cli::cli_rule()

  return(list(
    mapped = final_df,
    unmapped_primary = unmapped_genes_primary,
    unmapped_final = unmapped_genes_final,
    multi_mapped = multi_mapped_report
  ))
}


#' Map Ensembl/NCBI IDs to gene symbols
#'
#' Maps input Ensembl or NCBI gene IDs to symbols using a provided merged annotation.
#' Auto-detects ID type, supports multi-mapping strategies, and optional unmapped retention.
#'
#' @param input_data Character vector or data frame of IDs to map.
#' @param col_id If input_data is a data frame, the column containing IDs.
#' @param annotation The output of buildHumanAnnotation(), buildMouseAnnotation(), buildRatAnnotation() or getAnnotationBiomart()).
#' @param add_metadata Logical; if FALSE, return only query ID and gene_symbol.
#' @param multi_handling Strategy for multi-mapped IDs: "keepAll", "collapse", or "unique".
#' @param keep_unmapped Logical; if TRUE, retain rows that could not be mapped. Default is FALSE.
#' @return A list with elements:
#'   \item{mapped}{Tibble of mapped (and optionally unmapped) rows.}
#'   \item{unmapped}{Character vector of IDs that could not be mapped.}
#'   \item{multi_mapped}{Tibble summarizing IDs with multiple symbol matches.}
#' @export

mapIDToSymbol <- function(
  input_data,
  col_id = NULL,
  annotation,
  add_metadata = TRUE,
  multi_handling = c("keepAll", "collapse", "unique"),
  keep_unmapped = FALSE
) {
  multi_handling <- match.arg(multi_handling)

  # Validate required annotation columns
  if (!"gene_symbol" %in% names(annotation)) {
    cli::cli_abort("The {.arg annotation} table must contain column {.val gene_symbol}.")
  }

  # --- 1. Handle Vector vs Data Frame Input ---
  if (is.vector(input_data) || is.factor(input_data)) {
    original_df <- data.frame(temp_id = as.character(input_data), stringsAsFactors = FALSE) %>%
      mutate(row_id_internal = row_number())
    input_is_vector <- TRUE
  } else if (is.data.frame(input_data)) {
    if (is.null(col_id)) cli::cli_abort("Argument {.arg col_id} must be provided when input is a data frame.")
    if (!col_id %in% names(input_data)) cli::cli_abort("Column {.val {col_id}} not found in input data frame.")

    original_df <- input_data %>%
      dplyr::rename(temp_id = !!sym(col_id)) %>%
      mutate(row_id_internal = row_number())
    input_is_vector <- FALSE
  } else {
    cli::cli_abort("Input must be either a vector or a data frame.")
  }

  # --- 2. Auto-detect ID Type ---
  id_vector <- as.character(original_df$temp_id)

  valid_index <- which(!is.na(id_vector) &
    id_vector != "" &
    !tolower(trimws(id_vector)) %in% c("null", "na", "nan"))[1]

  if (is.na(valid_index)) {
    cli::cli_abort("All input IDs are missing (NA, empty, or 'null'). Cannot detect ID type.")
  }

  sample_id <- trimws(id_vector[valid_index])

  if (grepl("^ENS", sample_id, ignore.case = TRUE)) {
    id_type <- "ensembl_id"
    new_col_name <- "query_ensembl_id"
    cli::cli_alert_info("Detected ID type: {.val Ensembl} (found at row {valid_index})")
  } else if (grepl("^[0-9]+$", sample_id)) {
    id_type <- "ncbi_id"
    new_col_name <- "query_ncbi_id"
    cli::cli_alert_info("Detected ID type: {.val NCBI / Entrez} (found at row {valid_index})")
  } else {
    cli::cli_abort(paste0(
      "Could not determine ID type from sample: {.val ", sample_id, "}. ",
      "Ensure IDs are Ensembl (ENSG...) or NCBI (Numeric)."
    ))
  }

  if (!id_type %in% names(annotation)) {
    cli::cli_abort("The detected ID column {.val {id_type}} was not found in {.arg annotation}.")
  }

  # --- 3. Handle Name Collisions for Metadata/Symbols ---
  overlap_cols <- intersect(names(original_df), names(annotation))
  overlap_cols <- setdiff(overlap_cols, c("temp_id", "row_id_internal", id_type))

  if (length(overlap_cols) > 0) {
    original_df <- original_df %>%
      rename_with(~ paste0("query_", .), all_of(overlap_cols))
  }

  # --- 4. Standardize Annotation & Join Keys ---
  annotation <- annotation %>%
    mutate(across(everything(), as.character))

  to_map <- original_df %>%
    filter(!is.na(temp_id) &
      temp_id != "" &
      !tolower(trimws(temp_id)) %in% c("null", "na", "nan")) %>%
    mutate(join_id_tmp = if (id_type == "ensembl_id") toupper(sub("\\.\\d+$", "", trimws(temp_id))) else trimws(temp_id))

  no_id <- original_df %>% filter(!row_id_internal %in% to_map$row_id_internal)

  annot_subset <- annotation %>%
    filter(!is.na(!!sym(id_type))) %>%
    mutate(join_id_tmp = if (id_type == "ensembl_id") toupper(sub("\\.\\d+$", "", trimws(!!sym(id_type)))) else trimws(!!sym(id_type)))

  mapped_df <- to_map %>%
    left_join(
      annot_subset,
      by = "join_id_tmp",
      relationship = "many-to-many"
    ) %>%
    dplyr::select(-join_id_tmp)

  if (keep_unmapped) {
    final_df <- bind_rows(mapped_df, no_id)
  } else {
    final_df <- mapped_df %>% filter(!is.na(gene_symbol))
  }

  # --- 5. Detect and Handle Multiple Matches ---
  multi_mapped_report <- final_df %>%
    filter(!is.na(temp_id)) %>%
    group_by(row_id_internal) %>%
    filter(n_distinct(gene_symbol, na.rm = TRUE) > 1) %>%
    summarise(
      temp_id = first(temp_id),
      n_symbols = n_distinct(gene_symbol, na.rm = TRUE),
      symbols = paste(unique(na.omit(gene_symbol)), collapse = ","),
      .groups = "drop"
    ) %>%
    dplyr::select(-row_id_internal)

  if (nrow(multi_mapped_report) > 0) {
    cli::cli_alert_warning("{nrow(multi_mapped_report)} input IDs produced multiple matches.")
    cli::cli_alert_info("Handling strategy: {.val {multi_handling}}")

    if (multi_handling == "collapse") {
      final_df <- final_df %>%
        group_by(row_id_internal) %>%
        summarise(
          across(everything(), ~ {
            vals <- unique(na.omit(.x))
            if (length(vals) == 0) NA_character_ else paste(vals, collapse = "|")
          }),
          .groups = "drop"
        )
    } else if (multi_handling == "unique") {
      final_df <- final_df %>%
        group_by(row_id_internal) %>%
        arrange(is.na(gene_symbol), gene_symbol) %>%
        slice(1) %>%
        ungroup()
    }
  }

  # --- 6. Cleanup ---
  final_col_id <- if (!is.null(col_id)) col_id else id_type

  # If final_col_id already exists in final_df (from annotation), drop it first so renaming temp_id is clean
  if (final_col_id %in% names(final_df)) {
    final_df <- final_df %>% dplyr::select(-all_of(final_col_id))
  }

  final_df <- final_df %>%
    arrange(row_id_internal) %>%
    dplyr::select(-any_of(c("gene_synonyms", "row_id_internal"))) %>%
    dplyr::rename(!!sym(final_col_id) := temp_id)

  if (!add_metadata) {
    orig_names <- if (input_is_vector) final_col_id else names(input_data)
    current_names <- names(final_df)
    cols_to_keep <- current_names[current_names %in% orig_names | grepl("^query_", current_names)]

    final_df <- final_df %>% dplyr::select(any_of(cols_to_keep), gene_symbol)
  }

  final_df <- final_df %>%
    relocate(!!sym(final_col_id), gene_symbol)

  # --- 7. Summary Statistics ---
  total_input_count <- length(unique(original_df$temp_id))
  mapped_ids <- unique(final_df[[final_col_id]][!is.na(final_df$gene_symbol)])
  mapped_input_count <- length(mapped_ids)

  primary_empty_ids <- unique(as.character(no_id$temp_id))
  num_primary_empty <- length(primary_empty_ids)

  all_unique_ids <- unique(original_df$temp_id)
  unmapped_ids <- setdiff(all_unique_ids, c(mapped_ids, primary_empty_ids))
  total_unmapped <- length(unmapped_ids)

  perc_mapped <- if (total_input_count > 0) {
    round((mapped_input_count / total_input_count) * 100, 1)
  } else {
    0
  }

  cli::cli_h1("Mapping Summary (ID to Symbol)")
  cli::cli_alert_info("Total Unique Input IDs: {.val {total_input_count}}")
  cli::cli_alert_success("Successfully Mapped: {.val {mapped_input_count}} ({perc_mapped}%)")

  if (total_unmapped > 0 || num_primary_empty > 0) {
    if (total_unmapped > 0) {
      cli::cli_alert_danger(paste0(
        "Final Unmapped: {.val {total_unmapped}} ID{?s} ",
        "could not be mapped."
      ))
    }

    if (num_primary_empty > 0) {
      cli::cli_alert_warning(paste0(
        "{.val {num_primary_empty}} input ID{?s} ",
        "{?was/were} NA, empty, or 'null' and {?was/were} skipped."
      ))
    }
  }

  cli::cli_rule()

  return(list(
    mapped = final_df,
    unmapped = unmapped_ids,
    multi_mapped = multi_mapped_report
  ))
}
