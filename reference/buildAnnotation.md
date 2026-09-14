# Build Merged Gene Annotation for Multiple Species

A unified entry point to build or load gene annotations for Human,
Mouse, or Rat. This function automatically routes the request to the
species-specific builder while maintaining a consistent interface.

## Usage

``` r
buildAnnotation(
  species = c("human", "mouse", "rat"),
  annotation_dir = .get_data_path(),
  force_update = FALSE,
  offline = getOption("geneRosetta.offline", FALSE)
)
```

## Arguments

- species:

  Character string. One of "human", "mouse", or "rat".

- annotation_dir:

  Character string specifying the directory where annotation data files
  are located. Defaults to the persistent user data directory.

- force_update:

  Logical. If TRUE, forces re-download and rebuild of annotations even
  if a valid cache exists. Defaults to FALSE.

## Value

A dataframe containing the merged annotation for the selected species.

## Examples

``` r
if (FALSE) { # \dontrun{
# The simple way
human_annot <- buildAnnotation("human")

# Forcing an update for mouse
mouse_annot <- buildAnnotation("mouse", force_update = TRUE)
} # }
```
