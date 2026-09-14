# Getting Started with geneRosetta

## Overview

Comparing your gene lists across studies (especially older ones) can be
frustrating. Gene symbols change, IDs are deprecated, and different
databases sometimes disagree. This makes it hard to compare results
across studies, or even to reproduce your own analysis after a few
months.

`geneRosetta` solves this problem by building a multi-source reference
and uses it to map gene symbols and IDs. The package combines several
sources:

| Species | Sources                     |
|---------|-----------------------------|
| Human   | Ensembl, NCBI Gene and HGNC |
| Mouse   | Ensembl, NCBI Gene and MGI  |
| Rat     | Ensembl, NCBI Gene and RGD  |

Combining sources helps recover aliases and records missing from one
database.

Unlike packages with static, outdated reference data — or those that
need a live server connection for every lookup, which can fail in the
middle of an analysis — it downloads the references from every source in
real-time, merges them, and caches the merged reference locally.

Are you completely offline? No problem - use the built-in cached data.
Need to update your reference? Refresh it using the most recently
available data.

Furthermore, `geneRosetta` is designed to forgive real-world messy
inputs: it automatically trims surrounding whitespace, ignores case
differences, and strips version numbers from Ensembl IDs (e.g.,
converting `ENSG00000141510.5` to `ENSG00000141510`).

Species Note

Choose the species that matches your experiment. These functions map
identifiers within a species; they do not find equivalent orthologs
across species.

## Installation

``` r

if (!require("remotes")) install.packages("remotes")
remotes::install_github("alexankc/geneRosetta")

# or

if (!require("devtools")) install.packages("devtools")
devtools::install_github("alexankc/geneRosetta")
```

## Step 1: Build Merged Annotation

How geneRosetta builds Annotations

1.  **Multi-Tier Matching**

    First, Ensembl, NCBI, and official naming authorities (HGNC, MGI, or
    RGD) are linked using Ensembl IDs. If a gene isn’t found that way,
    the pipeline “rescues” it by looking for matches using NCBI IDs or
    the common gene symbol.

2.  **Smart Prioritization**

    When different databases disagree, the following rules are used:

    - **Official Symbols:** Species-specific authority symbols come
      first (such as HGNC, MGI, or RGD).
    - **Map Locations:** Ensembl coordinates are preferred (start, end,
      and chromosome).
    - **Synonym Cleaning:** All gene synonyms/aliases found across all
      three databases are gathered in a single list.

3.  **Extended Coverage**

    A record can be retained even when it lacks a current Ensembl ID
    (for example, if recognized by HGNC or NCBI). This broadens your
    lookup coverage, though some recognized genes may not have every ID
    type.

``` r

library(geneRosetta)

annot <- buildAnnotation(species = "human") # "human", "mouse", "rat" are supported
```

    ! No cache found (missing data or log). Creating merged annotation...

    Fetching annotations...

    ℹ Connecting to Ensembl...

    Ensembl site unresponsive, trying useast mirror

    Ensembl site unresponsive, trying asia mirror

    ℹ Downloading NCBI assembly report for human...

    ✔ Data and metadata log saved to disk.

    ℹ Downloading NCBI assembly report for human...
    ✔ Downloading NCBI assembly report for human... [473ms]

    ℹ Connecting to Ensembl...
    ✔ Connecting to Ensembl... [20.4s]

    ℹ Connecting to Ensembl and querying...
    ℹ Materializing bundled reference data into cache (may be outdated).
    ℹ Connecting to Ensembl and querying...
    ! Ensembl is currently unresponsive; using bundled package data.
    ℹ Connecting to Ensembl and querying...
    ✔ Connecting to Ensembl and querying... [11s]

    ℹ Downloading and parsing NCBI data...
    ℹ Downloading and parsing NCBI GFF3 (GCF_000001405.40_GRCh38.p14)...
    ✔ Cache is up-to-date (2024-09-28). Loading...
    ℹ Downloading and parsing NCBI GFF3 (GCF_000001405.40_GRCh38.p14)...
    ✔ Cache is up-to-date (2024-09-28). Loading...
    ℹ Downloading and parsing NCBI GFF3 (GCF_000001405.40_GRCh38.p14)...
    ✔ Downloading and parsing NCBI GFF3 (GCF_000001405.40_GRCh38.p14)... [1m 51.7s]

    ℹ Downloading and parsing NCBI data...
    ✔ Data and metadata log saved to disk.
    ℹ Downloading and parsing NCBI data...
    ✔ Downloading and parsing NCBI data... [1m 58.4s]

    ℹ Downloading and processing HGNC data...
    ✔ Data and metadata log saved to disk.
    ℹ Downloading and processing HGNC data...
    ✔ Downloading and processing HGNC data... [4.1s]

    Merging annotations...
    ────────────────────────────────────────────────────────────────────────────────
    ✔ Build Complete
    ────────────────────────────────────────────────────────────────────────────────

The annotation file and corresponding log file are both stored in the
user’s local cache directory by default (using
[`tools::R_user_dir()`](https://rdrr.io/r/tools/userdir.html)). This
ensures data is saved across R sessions while keeping the project folder
clean. You can optionally store files in a custom directory by setting
`annotation_dir`.

This cached annotation will be sufficiently updated for some time.
However, if you want to refresh it with the latest online data, run:

``` r

annot <- buildAnnotation(species = "human", force_update = TRUE)
```

#### Working completely offline

If you are working on a cluster node without internet access or need
strictly reproducible, air-gapped execution, use `offline = TRUE`:

``` r

annot <- buildAnnotation(species = "human", offline = TRUE)
```

Offline mode will use your locally cached data or bundled package data
where available, without making any external network requests.

Internally
[`buildAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildAnnotation.md)
calls:

- For species=human:
  [`buildHumanAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildHumanAnnotation.md)
  which in turn calls
  `getAnnotationBiomart(species = "human", version = "current")`,
  `getAnnotationNCBI(species = "human")` and
  [`getAnnotationHGNChuman()`](https://alexankc.github.io/geneRosetta/reference/getAnnotationHGNChuman.md).
- For species=mouse:
  [`buildMouseAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildMouseAnnotation.md)
  which in turn calls
  `getAnnotationBiomart(species = "mouse", version = "current")`,
  `getAnnotationNCBI(species = "mouse")` and
  [`getAnnotationMGImouse()`](https://alexankc.github.io/geneRosetta/reference/getAnnotationMGImouse.md).
- For species=rat:
  [`buildRatAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildRatAnnotation.md)
  which in turn calls
  `getAnnotationBiomart(species = "rat", version = "current")`,
  `getAnnotationNCBI(species = "rat")` and
  `getAnnotationRGDrat(genome_build = "GRCr8")`.

Important

Only the output of
[`buildAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildAnnotation.md)
or
[`getAnnotationBiomart()`](https://alexankc.github.io/geneRosetta/reference/getAnnotationBiomart.md)
should be used for downstream mapping with `geneRosetta`.

#### Reference Anatomy

The merged table produced by
[`buildAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildAnnotation.md)
contains standardized columns across databases:

| Column | Description |
|----|----|
| `gene_symbol` | Current official gene symbol (prioritizing naming authority) |
| `ensembl_id` | Ensembl gene identifier |
| `ncbi_id` | NCBI Entrez Gene identifier |
| `hgnc_id` / `mgi_id` / `rgd_id` | Species-specific authority accession |
| `gene_biotype` | Gene classification (e.g., protein_coding, lncRNA, pseudogene) |
| `chromosome`, `gene_start`, `gene_end`, `strand` | Genomic coordinates and orientation (stored as character strings) |
| `metadata_source` | Source database providing the coordinates and metadata |
| `gene_synonyms` | Aggregated aliases and past symbols, separated by `\|` |

#### Chromosome labels and `chromosome_scope`

Genes on standard chromosomes are reported with bare names (`"14"`,
`"X"`, `"MT"`). Everything else - unplaced and unlocalized scaffolds,
alternate haplotypes, and patch releases - is labelled
`chr<molecule>(<GenBank accession>)`, resolved against the NCBI assembly
report for the exact genome build Ensembl is serving.

So the unlocalized scaffold carrying `MAFIP` is reported as
`chr14(GL000194.1)`, and a scaffold with no assigned molecule becomes
`chrUn(<accession>)`. The label always carries the GenBank accession,
making it identifiable across Ensembl versions.

These sequences are included by default, but you can narrow them using
the `chromosome_scope` argument in
[`getAnnotationBiomart()`](https://alexankc.github.io/geneRosetta/reference/getAnnotationBiomart.md):

| Scope | Keeps |
|----|----|
| `"all"` (default) | All sequences |
| `"primary"` | Sequences in standard chromosomes + unlocalized/unplaced scaffolds, excluding alternate haplotype scaffolds (e.g., HLA regions) and patches |
| `"standard"` | Sequences in standard chromosomes only (1-22, X, Y, MT) |

``` r

# Everything Ensembl reports (default)
all_genes <- getAnnotationBiomart("human")

# Drop alternate haplotypes and patches, but keep scaffolds like chr14(GL000194.1)
primary_only <- getAnnotationBiomart("human", chromosome_scope = "primary")

# Standard chromosomes only
standard_only <- getAnnotationBiomart("human", chromosome_scope = "standard")
```

Filtering Merged Annotations

[`buildAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildAnnotation.md)
always merges with `chromosome_scope = "all"`. Because the `chromosome`
column is self-describing, you can easily filter the merged table
afterwards:

``` r

# Filter out scaffolds, alternate haplotypes, and patches:
annot_standard <- subset(annot, !is.na(chromosome) & !grepl("^chr", chromosome))
```

#### Version control & Manuscript Methods

Annotations change frequently, so keeping track of metadata matters. Log
files record what was downloaded, when, and from which source release.
***Run this after building your annotation***.

The log summary can be printed to the console or exported in various
formats:

| `export_format` | Description |
|----|----|
| `"none"` (default) | Interactive formatted summary in the R console |
| `"txt"` | Human-readable text report |
| `"json"` | Full structured JSON log with all source URLs, timestamps, and row counts |
| `"methods"` | Auto-generated Methods paragraph tailored for manuscripts |

``` r

# Text log for your records
getAnnotationLogSummary(species = "human", export_format = "txt", export_path = "logs/merged_annotation.txt")

# Export structured JSON
getAnnotationLogSummary(species = "human", export_format = "json", export_path = "logs/merged_annotation.json")

# Generate a publication-ready Methods paragraph for your manuscript:
getAnnotationLogSummary(species = "human", export_format = "methods", export_path = "logs/annotation_methods.txt")
```

The `"methods"` export produces a clean paragraph describing the
database releases, access dates, and build procedures, which you can
paste directly into your paper’s Methods section:

``` r

getAnnotationLogSummary(species = "human", export_format = "methods")
```

    ✔ Methods text exported to 'human_annotation_methods.txt'

``` r

cat(readLines("human_annotation_methods.txt", warn = FALSE), sep = "\n")
```

    Human gene annotations were generated using geneRosetta (v1.0.0; https://github.com/alexankc/geneRosetta) by integrating Ensembl release 116 (GRCh38.p14), the NCBI Gene database (accessed September 2026), and the HGNC complete gene set (accessed September 2026). Records were first matched by Ensembl gene IDs, with unmatched records cross-referenced against NCBI Gene IDs; HGNC-approved gene symbols were preferred as primary symbols, while conflicting source symbols were retained as synonyms. Coding genes, non-coding RNAs, and pseudogenes were retained; records lacking both an Ensembl and an NCBI identifier were excluded, resulting in a final merged reference set containing 234,128 unique entries.

For reproducibility reasons, you can pin a specific Ensembl archive
release via `biomaRt`:

``` r

annot_biomart <- getAnnotationBiomart(
  species = "human", # "human", "mouse", "rat" are supported
  version = "114"    # pins Ensembl release 114
)
getAnnotationLogSummary(species = "human", export_format = "txt", ensembl_annot = TRUE, export_path = "logs/biomart_annotation.txt")
```

If that specific archive is temporarily unavailable (biomaRt servers
down again…), the latest saved cached annotation is used.

## Step 2: Map gene symbols ↔︎ gene IDs (Ensembl/NCBI)

### Gene Symbols → Gene IDs

Here, we face the big problem of frequent gene symbol changes.
Nomenclature changes whenever the scientific community discovers a new
function - or - whenever Excel decides a gene name looks suspiciously
like a calendar date.

This is why mapping is performed in two phases (passes):

1.  **Exact Match:** Matching the input symbol as provided.
2.  **Alias Search:** Searching the input against known gene synonyms
    and historical aliases.

For example, although `p53` is a widely used gene symbol for Tumor
Protein P53, the official symbol is `TP53`. If `p53` is searched against
a strict official database, there would not be any match. With the
second pass, the symbol is “rescued”.

Another challenge is that one symbol can sometimes map to multiple gene
IDs. In that case, choose a strategy via `multi_handling`:

- `keepAll`: one row per match (all mappings preserved; may expand rows)
- `collapse`: one row per input, values joined with `|` (cleaner table;
  prevents row expansion)
- `unique`: keep one mapping only (first sorted match; prevents row
  expansion)

When `multi_handling = "unique"` is used for symbol mapping, it
prioritizes records with a non-missing `ensembl_id` and sorts by
`ensembl_id`.

Multi-mapping typically has two distinct biological causes:

1.  **Genuinely distinct loci** — for example `CD99`, which sits in the
    pseudoautosomal region (PAR) and has separate functional IDs on the
    X and Y chromosomes.
2.  **Alternate haplotype copies** — genes such as `HLA-A` appear once
    on the primary assembly (chromosome `6`) and again on each
    alternative haplotype scaffold (e.g., `chr6(GL000256.2)`). Their
    `chromosome` values carry the `chr6(...)` form, making them easy to
    spot and filter out.

#### Vectors as input

In this example, we include:

- `" p53 "`: gene whose name has been updated (`TP53`), with
  leading/trailing spaces to test whitespace trimming
- `INS`: standard gene symbol
- `CD99`: gene with multiple Ensembl IDs on X & Y
- `FAKEGENE`: intentionally invalid gene symbol

``` r

id_lookup <- mapSymbolToID(
    input_data = c(" p53 ", "INS", "CD99", "FAKEGENE"), # input as vector
    annotation = annot,
    add_metadata = FALSE
)
```

    ℹ Checking annotation synonyms...

    ! 1 symbols produced multiple matches.

    ℹ Handling strategy: "keepAll"

    ── Mapping Summary (Symbol to ID) ──────────────────────────────────────────────

    ℹ Total Unique Input Genes: 4

    ✔ Successfully Mapped: 3 (75%) - of which 1 was mapped via synonyms

    ✖ Final Unmapped: 1 gene could not be mapped.

    ────────────────────────────────────────────────────────────────────────────────

Console output summarizes what happened.

The returned list object contains:

##### Mapped symbols

``` r

id_mapped <- id_lookup$mapped # mapped genes

id_mapped
```

| query_gene_symbol | gene_symbol | ensembl_id      |
|:------------------|:------------|:----------------|
| p53               | TP53        | ENSG00000141510 |
| INS               | INS         | ENSG00000254647 |
| CD99              | CD99        | ENSG00000002586 |
| CD99              | CD99        | ENSG00000292348 |

Input symbols are stored in `query_gene_symbol`, and current official
symbols in `gene_symbol` (showing updates like `p53` → `TP53`).

##### Symbols with no match in 1st pass (before alias rescue)

``` r

id_unmapped_firstPass <- id_lookup$unmapped_primary # genes not mapped in the first passage
id_unmapped_firstPass
```

    [1] " p53 "    "FAKEGENE"

##### Symbols with no match in 1st + 2nd pass

``` r

id_unmapped_final <- id_lookup$unmapped_final # completely unmapped genes
id_unmapped_final
```

    [1] "FAKEGENE"

##### Symbols with multiple ID matches

``` r

id_multi_mapped <- id_lookup$multi_mapped # genes that mapped to multiple IDs
id_multi_mapped
```

| query_gene_symbol | ensembl_ids                     |
|:------------------|:--------------------------------|
| CD99              | ENSG00000002586,ENSG00000292348 |

You can include extra metadata (NCBI + naming authority IDs, biotype,
genomic coordinates, description) with `add_metadata = TRUE` (default):

``` r

id_lookup <- mapSymbolToID(
    input_data = c(" p53 ", "INS", "CD99", "FAKEGENE"),
    annotation = annot,
    add_metadata = TRUE 
)
```

    ℹ Checking annotation synonyms...

    ! 1 symbols produced multiple matches.

    ℹ Handling strategy: "keepAll"

    ── Mapping Summary (Symbol to ID) ──────────────────────────────────────────────

    ℹ Total Unique Input Genes: 4

    ✔ Successfully Mapped: 3 (75%) - of which 1 was mapped via synonyms

    ✖ Final Unmapped: 1 gene could not be mapped.

    ────────────────────────────────────────────────────────────────────────────────

``` r

id_lookup$mapped
```

| query_gene_symbol | gene_symbol | ensembl_id | gene_biotype | chromosome | gene_start | gene_end | strand | metadata_source | ncbi_id | hgnc_id | gene_description |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| p53 | TP53 | ENSG00000141510 | protein_coding | 17 | 7661779 | 7687546 | -1 | Ensembl | 7157 | 11998 | tumor protein p53 \[Source:HGNC Symbol;Acc:HGNC:11998\] |
| INS | INS | ENSG00000254647 | protein_coding | 11 | 2159779 | 2161221 | -1 | Ensembl | 3630 | 6081 | insulin \[Source:HGNC Symbol;Acc:HGNC:6081\] |
| CD99 | CD99 | ENSG00000002586 | protein_coding | X | 2690988 | 2741312 | 1 | Ensembl | 4267 | 7082 | CD99 molecule (Xg blood group) \[Source:HGNC Symbol;Acc:HGNC:7082\] |
| CD99 | CD99 | ENSG00000292348 | protein_coding | Y | 2690988 | 2741312 | 1 | Ensembl | 4267 | 7082 | CD99 molecule (Xg blood group) \[Source:HGNC Symbol;Acc:HGNC:7082\] |

Because `CD99` maps to multiple Ensembl IDs, `keepAll` prints one row
per match. Use `multi_handling = "collapse"` for one-row summaries:

``` r

id_lookup <- mapSymbolToID(
    input_data = c(" p53 ", "INS", "CD99", "FAKEGENE"),
    annotation = annot,
    add_metadata = TRUE,
  multi_handling = "collapse"
)
```

    ℹ Checking annotation synonyms...

    ! 1 symbols produced multiple matches.

    ℹ Handling strategy: "collapse"

    ── Mapping Summary (Symbol to ID) ──────────────────────────────────────────────

    ℹ Total Unique Input Genes: 4

    ✔ Successfully Mapped: 3 (75%) - of which 1 was mapped via synonyms

    ✖ Final Unmapped: 1 gene could not be mapped.

    ────────────────────────────────────────────────────────────────────────────────

``` r

id_lookup$mapped
```

| query_gene_symbol | gene_symbol | ensembl_id | gene_biotype | chromosome | gene_start | gene_end | strand | metadata_source | ncbi_id | hgnc_id | gene_description |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| p53 | TP53 | ENSG00000141510 | protein_coding | 17 | 7661779 | 7687546 | -1 | Ensembl | 7157 | 11998 | tumor protein p53 \[Source:HGNC Symbol;Acc:HGNC:11998\] |
| INS | INS | ENSG00000254647 | protein_coding | 11 | 2159779 | 2161221 | -1 | Ensembl | 3630 | 6081 | insulin \[Source:HGNC Symbol;Acc:HGNC:6081\] |
| CD99 | CD99 | ENSG00000002586\|ENSG00000292348 | protein_coding | X\|Y | 2690988 | 2741312 | 1 | Ensembl | 4267 | 7082 | CD99 molecule (Xg blood group) \[Source:HGNC Symbol;Acc:HGNC:7082\] |

Now multiple values are joined with `|` in a single row.

Note on Collapsed Columns

When `multi_handling = "collapse"` is used, each column is collapsed
independently. If strict row-by-row linkage between IDs and coordinates
is required, keep `multi_handling = "keepAll"`.

By default, unmapped genes are omitted from the `mapped` table (they are
always accessible in `id_lookup$unmapped_final`). Set
`keep_unmapped = TRUE` to keep them in the output table with `NA`
values:

``` r

id_lookup <- mapSymbolToID(
    input_data = c(" p53 ", "INS", "CD99", "FAKEGENE"),
    annotation = annot,
    add_metadata = TRUE,
  multi_handling = "collapse",
  keep_unmapped = TRUE
)
```

    ℹ Checking annotation synonyms...

    ! 1 symbols produced multiple matches.

    ℹ Handling strategy: "collapse"

    ── Mapping Summary (Symbol to ID) ──────────────────────────────────────────────

    ℹ Total Unique Input Genes: 4

    ✔ Successfully Mapped: 3 (75%) - of which 1 was mapped via synonyms

    ✖ Final Unmapped: 1 gene could not be mapped.

    ────────────────────────────────────────────────────────────────────────────────

``` r

id_lookup$mapped
```

| query_gene_symbol | gene_symbol | ensembl_id | gene_biotype | chromosome | gene_start | gene_end | strand | metadata_source | ncbi_id | hgnc_id | gene_description |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| p53 | TP53 | ENSG00000141510 | protein_coding | 17 | 7661779 | 7687546 | -1 | Ensembl | 7157 | 11998 | tumor protein p53 \[Source:HGNC Symbol;Acc:HGNC:11998\] |
| INS | INS | ENSG00000254647 | protein_coding | 11 | 2159779 | 2161221 | -1 | Ensembl | 3630 | 6081 | insulin \[Source:HGNC Symbol;Acc:HGNC:6081\] |
| CD99 | CD99 | ENSG00000002586\|ENSG00000292348 | protein_coding | X\|Y | 2690988 | 2741312 | 1 | Ensembl | 4267 | 7082 | CD99 molecule (Xg blood group) \[Source:HGNC Symbol;Acc:HGNC:7082\] |
| FAKEGENE | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

#### Dataframes as input

You can also pass a full dataframe (such as differential expression
results) and specify the symbol column using `col_id`. This is handy
when your table already contains statistics (for example, gene
expression results) and you want your results and the annotation data in
one step.

``` r

test_df_symbols <- data.frame(
  comparison = c("A_vs_B", "A_vs_C", "A_vs_B", "A_vs_B", "A_vs_B", "A_vs_B"),
  symbol_input = c("INS", "INS", "C1orf123", "POLR2J3", "FAKE", NA), 
  log2FC = c(2.5, -1.2, 0.5, 0.2, 0.5, 0.6),
  p_val = c(0.001, 0.05, 0.1, 0.65, 0.8, 1),
  stringsAsFactors = FALSE
)

test_df_annotated <- mapSymbolToID(
  input_data = test_df_symbols,
  annotation = annot,
  col_id = "symbol_input", 
  add_metadata = TRUE, 
  multi_handling = "collapse", 
  keep_unmapped = TRUE
)
```

    ℹ Checking annotation synonyms...

    ! 1 symbols produced multiple matches.

    ℹ Handling strategy: "collapse"

    ── Mapping Summary (Symbol to ID) ──────────────────────────────────────────────

    ℹ Total Unique Input Genes: 5

    ✔ Successfully Mapped: 3 (60%) - of which 1 was mapped via synonyms

    ✖ Final Unmapped: 1 gene could not be mapped.

    ! 1 input gene was NA, empty, or 'null' and was skipped.

    ────────────────────────────────────────────────────────────────────────────────

``` r

test_df_annotated$mapped
```

| query_gene_symbol | gene_symbol | ensembl_id | comparison | log2FC | p_val | gene_biotype | chromosome | gene_start | gene_end | strand | metadata_source | ncbi_id | hgnc_id | gene_description |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| INS | INS | ENSG00000254647 | A_vs_B | 2.5 | 0.001 | protein_coding | 11 | 2159779 | 2161221 | -1 | Ensembl | 3630 | 6081 | insulin \[Source:HGNC Symbol;Acc:HGNC:6081\] |
| INS | INS | ENSG00000254647 | A_vs_C | -1.2 | 0.05 | protein_coding | 11 | 2159779 | 2161221 | -1 | Ensembl | 3630 | 6081 | insulin \[Source:HGNC Symbol;Acc:HGNC:6081\] |
| C1orf123 | CZIB | ENSG00000162384 | A_vs_B | 0.5 | 0.1 | protein_coding | 1 | 53214099 | 53220657 | -1 | Ensembl | 54987 | 26059 | CXXC motif containing zinc binding protein \[Source:HGNC Symbol;Acc:HGNC:26059\] |
| POLR2J3 | POLR2J3 | ENSG00000168255\|ENSG00000285437 | A_vs_B | 0.2 | 0.65 | protein_coding | 7 | 102537918\|102552328 | 102572653\|102572620 | -1 | Ensembl | 548644 | 33853 | RNA polymerase II subunit J3 \[Source:HGNC Symbol;Acc:HGNC:33853\] |
| FAKE | NA | NA | A_vs_B | 0.5 | 0.8 | NA | NA | NA | NA | NA | NA | NA | NA | NA |
| NA | NA | NA | A_vs_B | 0.6 | 1 | NA | NA | NA | NA | NA | NA | NA | NA | NA |

Tips for Annotating Dataframes

- **Handling Name Collisions:** If your input dataframe already has
  columns with the same names as the annotation (`ensembl_id`,
  `gene_biotype`, `chromosome`, `gene_start`, `gene_end`, `strand`,
  `gene_description`), the `query_` prefix is automatically added to the
  original column names, so you can compare old and new values.
- **Tracking Expanded Rows:** If you choose `multi_handling = "keepAll"`
  on a dataframe, multi-mapped genes will expand into multiple rows.
  Having an optional row or observation ID in your table can make it
  easier to trace expanded rows back and avoid overcounting measurements
  in downstream enrichment tests.

### Gene IDs → Gene Symbols

[`mapIDToSymbol()`](https://alexankc.github.io/geneRosetta/reference/mapIDToSymbol.md)
accepts a vector or dataframe column of **Ensembl gene IDs** or **NCBI
Gene (Entrez) IDs** and auto-detects the ID type.

Key details to know about ID mapping:

- **Single-Pass:** Unlike symbol mapping (which has an alias rescue
  pass), ID mapping is a single-pass lookup.
- **Suffix & Space Trimming:** Version numbers on Ensembl IDs (e.g.,
  `.5` in `ENSG00000141510.5`) are automatically stripped, and
  leading/trailing whitespace is removed.
- **Default Output:** Like
  [`mapSymbolToID()`](https://alexankc.github.io/geneRosetta/reference/mapSymbolToID.md),
  `keep_unmapped = FALSE` by default. Set `keep_unmapped = TRUE` if you
  wish to retain unmapped IDs in the output table.

#### Mapping Ensembl IDs

In this example, we include a mix of versioned, whitespace,
multi-mapped, and unmapped IDs:

- `"ENSG00000141510.5"`: versioned Ensembl ID (`TP53`)
- `" ensg00000012048 "`: Ensembl ID with whitespace and lowercase
  formatting (`BRCA1`)
- `"ENSG00000226087"`: multi-mapped Ensembl ID (maps to both
  `LOC105374588` and `LOC124907763`)
- `"ENSG00000000000"`: non-existent/unmapped Ensembl ID

``` r

genes_df <- data.frame(
    gene_ids = c("ENSG00000141510.5", " ensg00000012048 ", "ENSG00000226087", "ENSG00000000000"),
    scores = c(5, 3, 2, 1)
)

symbol_lookup <- mapIDToSymbol(
    input_data = genes_df, # input as dataframe
    col_id = "gene_ids",
    annotation = annot,
    multi_handling = "keepAll",
    keep_unmapped = TRUE
)
```

    ℹ Detected ID type: "Ensembl" (found at row 1)

    ! 1 input IDs produced multiple matches.

    ℹ Handling strategy: "keepAll"

    ── Mapping Summary (ID to Symbol) ──────────────────────────────────────────────

    ℹ Total Unique Input IDs: 4

    ✔ Successfully Mapped: 3 (75%)

    ✖ Final Unmapped: 1 ID could not be mapped.

    ────────────────────────────────────────────────────────────────────────────────

``` r

symbol_lookup$mapped
```

| gene_ids | gene_symbol | scores | ensembl_id | gene_biotype | chromosome | gene_start | gene_end | strand | metadata_source | ncbi_id | hgnc_id | gene_description |
|:---|:---|---:|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| ENSG00000141510.5 | TP53 | 5 | ENSG00000141510 | protein_coding | 17 | 7661779 | 7687546 | -1 | Ensembl | 7157 | 11998 | tumor protein p53 \[Source:HGNC Symbol;Acc:HGNC:11998\] |
| ensg00000012048 | BRCA1 | 3 | ENSG00000012048 | protein_coding | 17 | 43044292 | 43170245 | -1 | Ensembl | 672 | 1100 | BRCA1 DNA repair associated \[Source:HGNC Symbol;Acc:HGNC:1100\] |
| ENSG00000226087 | LOC105374588 | 2 | ENSG00000226087 | lncRNA | 2 | 47175209 | 47240886 | 1 | Ensembl | 105374588 | NA | novel transcript |
| ENSG00000226087 | LOC124907763 | 2 | ENSG00000226087 | lncRNA | 2 | 47175209 | 47240886 | 1 | Ensembl | 124907763 | NA | novel transcript |
| ENSG00000000000 | NA | 1 | NA | NA | NA | NA | NA | NA | NA | NA | NA | NA |

Like
[`mapSymbolToID()`](https://alexankc.github.io/geneRosetta/reference/mapSymbolToID.md),
the returned object carries the mapped table and diagnostic diagnostics:

``` r

symbol_lookup$unmapped # IDs with no matching gene symbol
```

    [1] "ENSG00000000000"

``` r

symbol_lookup$multi_mapped # IDs matching more than one symbol
```

| temp_id         | n_symbols | symbols                   |
|:----------------|----------:|:--------------------------|
| ENSG00000226087 |         2 | LOC105374588,LOC124907763 |

#### Mapping NCBI / Entrez IDs

If input is a number/numeric vector, it is assumed to represent NCBI /
Entrez IDs.

``` r

ncbi_lookup <- mapIDToSymbol(
  input_data = c("7157", "672"), # TP53 and BRCA1
  annotation = annot,
  add_metadata = TRUE
)
```

    ℹ Detected ID type: "NCBI / Entrez" (found at row 1)

    ── Mapping Summary (ID to Symbol) ──────────────────────────────────────────────

    ℹ Total Unique Input IDs: 2

    ✔ Successfully Mapped: 2 (100%)

    ────────────────────────────────────────────────────────────────────────────────

``` r

ncbi_lookup$mapped
```

| ncbi_id | gene_symbol | ensembl_id | gene_biotype | chromosome | gene_start | gene_end | strand | metadata_source | hgnc_id | gene_description |
|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|:---|
| 7157 | TP53 | ENSG00000141510 | protein_coding | 17 | 7661779 | 7687546 | -1 | Ensembl | 11998 | tumor protein p53 \[Source:HGNC Symbol;Acc:HGNC:11998\] |
| 672 | BRCA1 | ENSG00000012048 | protein_coding | 17 | 43044292 | 43170245 | -1 | Ensembl | 1100 | BRCA1 DNA repair associated \[Source:HGNC Symbol;Acc:HGNC:1100\] |

## Troubleshooting & Tips

| What you see | Likely cause & Recommended action |
|----|----|
| **Many unmapped symbols** | Confirm species matches your data. Check if your inputs are IDs rather than symbols, or if a spreadsheet tool altered names to dates (e.g., `SEPT4` $`\rightarrow`$`Sept-04`). |
| **Row count increased after mapping** | When using `multi_handling = "keepAll"`, multi-mapped loci create multiple rows. Filter by standard chromosomes or track rows using an observation ID. |
| **Gene symbol found, but `ensembl_id` is blank** | HGNC or NCBI recognized the gene, but no Ensembl link currently exists. Inspect `ncbi_id` or authority IDs in the metadata. |
| **ID-type detection error** | Ensure you do not mix Ensembl IDs and NCBI IDs in the same vector. Map each identifier type in a separate call. |
| **Log summary finds nothing** | Make sure the `annotation_dir` parameter passed to [`getAnnotationLogSummary()`](https://alexankc.github.io/geneRosetta/reference/getAnnotationLogSummary.md) matches the directory used in [`buildAnnotation()`](https://alexankc.github.io/geneRosetta/reference/buildAnnotation.md). |

## Session Information

``` r

sessionInfo()
```

    R version 4.6.1 (2026-06-24)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.5 LTS

    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0

    locale:
     [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
     [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
     [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C

    time zone: UTC
    tzcode source: system (glibc)

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] geneRosetta_1.0.0

    loaded via a namespace (and not attached):
     [1] tidyr_1.3.2          generics_0.1.4       xml2_1.6.0
     [4] RSQLite_3.53.3       stringi_1.8.9        hms_1.1.4
     [7] digest_0.6.39        magrittr_2.0.5       evaluate_1.0.5
    [10] fastmap_1.2.0        blob_1.3.0           jsonlite_2.0.0
    [13] progress_1.2.3       AnnotationDbi_1.74.0 DBI_1.3.0
    [16] httr_1.4.9           purrr_1.2.2          Biostrings_2.80.2
    [19] httr2_1.3.0          cli_3.6.6            rlang_1.3.0
    [22] crayon_1.5.3         XVector_0.52.0       dbplyr_2.6.0
    [25] Biobase_2.72.0       bit64_4.8.6          withr_3.0.3
    [28] cachem_1.1.0         yaml_2.3.12          otel_0.2.0
    [31] parallel_4.6.1       tools_4.6.1          tzdb_0.5.0
    [34] memoise_2.0.1        dplyr_1.2.1          filelock_1.0.3
    [37] BiocGenerics_0.58.1  curl_8.0.0           vctrs_0.7.3
    [40] R6_2.6.1             png_0.1-9            stats4_4.6.1
    [43] BiocFileCache_3.2.0  lifecycle_1.0.5      Seqinfo_1.2.0
    [46] KEGGREST_1.52.2      stringr_1.6.0        S4Vectors_0.50.2
    [49] IRanges_2.46.0       bit_4.6.0            vroom_1.7.1
    [52] pkgconfig_2.0.3      pillar_1.11.1        glue_1.8.1
    [55] xfun_0.60            tibble_3.3.1         tidyselect_1.2.1
    [58] knitr_1.52           htmltools_0.5.9      rmarkdown_2.32
    [61] readr_2.2.0          compiler_4.6.1       prettyunits_1.2.0
    [64] biomaRt_2.68.0      
