# Build Scaffold Name Lookups From an Assembly Report

Ensembl names scaffolds inconsistently: unlocalized and unplaced
scaffolds are named by their GenBank accession ("GL000194.1"), while
alt-scaffolds and patches are named by their sequence name ("HG1_PATCH",
"HSCHR6_MHC_APD_CTG1"). Indexing every alias a report offers therefore
resolves all of them.

## Usage

``` r
.chromosome_lookup(report)
```

## Arguments

- report:

  A tibble from
  [`.get_assembly_report()`](https://alexankc.github.io/geneRosetta/reference/dot-get_assembly_report.md).

## Value

A list with two named character vectors:

- label:

  alias -\> `chr<molecule>(<accession>)` label

- acc_role:

  GenBank accession -\> sequence role
