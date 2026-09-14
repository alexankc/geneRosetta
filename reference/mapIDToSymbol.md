# Map Ensembl/NCBI IDs to gene symbols

Maps input Ensembl or NCBI gene IDs to symbols using a provided merged
annotation. Auto-detects ID type, supports multi-mapping strategies, and
optional unmapped retention.

## Usage

``` r
mapIDToSymbol(
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

  Character vector or data frame of IDs to map.

- col_id:

  If input_data is a data frame, the column containing IDs.

- annotation:

  The output of buildHumanAnnotation(), buildMouseAnnotation(),
  buildRatAnnotation() or getAnnotationBiomart()).

- add_metadata:

  Logical; if FALSE, return only query ID and gene_symbol.

- multi_handling:

  Strategy for multi-mapped IDs: "keepAll", "collapse", or "unique".

- keep_unmapped:

  Logical; if TRUE, retain rows that could not be mapped. Default is
  FALSE.

## Value

A list with elements:

- mapped:

  Tibble of mapped (and optionally unmapped) rows.

- unmapped:

  Character vector of IDs that could not be mapped.

- multi_mapped:

  Tibble summarizing IDs with multiple symbol matches.
