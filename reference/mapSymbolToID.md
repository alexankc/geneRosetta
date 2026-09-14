# Map gene symbols to Ensembl IDs

Maps input gene symbols to Ensembl IDs using a provided merged
annotation. Supports gene synonym rescue, configurable handling of
multi-mapped symbols, and optional retention of unmapped rows.

## Usage

``` r
mapSymbolToID(
  input_data,
  col_id = NULL,
  annotation,
  add_metadata = TRUE,
  multi_handling = c("keepAll", "collapse", "unique"),
  keep_unmapped = FALSE
)
```

## Arguments

- input_data:

  Character vector or data frame of gene symbols to map.

- col_id:

  If input_data is a data frame, the column containing the gene symbols.

- annotation:

  The output of buildHumanAnnotation(), buildMouseAnnotation(),
  buildRatAnnotation() or getAnnotationBiomart()).

- add_metadata:

  Logical; if FALSE, return only query, ensembl_id, gene_symbol.

- multi_handling:

  Strategy for multi-mapped symbols: "keepAll", "collapse", or "unique".

- keep_unmapped:

  Logical; if TRUE, retain rows that could not be mapped. Default is
  FALSE.

## Value

A list with elements:

- mapped:

  Tibble of mapped (and optionally unmapped) rows.

- unmapped_primary:

  Character vector of symbols unmapped after primary pass.

- unmapped_final:

  Character vector of symbols still unmapped after synonym rescue.

- multi_mapped:

  Tibble summarizing symbols with multiple matches.
