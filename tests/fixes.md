# Including non-standard chromosomes in geneRosetta

Record of the change that brought scaffold, alt-haplotype and patch genes into the
annotation tables, and labelled them as e.g. `chr14(GL000194.1)`.

Reference test case: **MAFIP must map to `chr14(GL000194.1)`.**

## What was actually blocking MAFIP

The original report was that all five `getAnnotation*()` functions restricted
themselves to standard chromosomes. Checking each one, only a single function was
really filtering:

| Function | Filtered by chromosome? | Action taken |
| --- | --- | --- |
| `getAnnotationBiomart()` | Yes — `filter(chromosome %in% c(1:22, "X", "Y", "MT"))` | Filter removed, scaffolds labelled |
| `getAnnotationNCBI()` | No | None needed |
| `getAnnotationMGImouse()` | No | Labelling applied (no-op today) |
| `getAnnotationRGDrat()` | No | Labelling applied (no-op today) |
| `getAnnotationHGNChuman()` | No, but read `non_alt_loci_set.txt` | Switched to `hgnc_complete_set.txt` |

`getAnnotationBiomart()` was the sole gate, for two reasons:

1. All three merge builders are **BioMart-anchored** — every join in
   `mergeAnnotations.R` hangs off `df_bm`. A gene absent from BioMart can never
   reach the merged table.
2. MGI and RGD explicitly drop their own chromosome column (`ensembl_dup_cols`)
   in favour of Ensembl's, so the merged `chromosome` column **always** comes
   from BioMart.

MAFIP (`ENSG00000274847`) lives on `GL000194.1`, so the filter removed it from
BioMart and it disappeared from the merge — even though it was present in both
NCBI (`727764`) and HGNC (`31102`) all along.

HGNC's `non_alt_loci_set.txt` was not filtering by chromosome, but by
construction it omits genes on alternate reference loci. Switching to
`hgnc_complete_set.txt` adds 47 genes, including `HLA-DRB3` and `HLA-DRB4`.

## How the labelling works

New module: `R/chromosomes.R`.

It downloads the NCBI assembly report for whatever build Ensembl is currently
serving, and rewrites non-standard chromosome names to
`chr<assigned-molecule>(<GenBank accession>)`.

- Ensembl's `listDatasets()$version` string (`GRCh38.p14`, `GRCm39`, `GRCr8`)
  matches NCBI's assembly-directory suffix exactly, so the correct report is
  found by listing the species' GCA directory. A pinned accession is used as
  fallback if the listing is unreachable or the build is not yet published.
- Assembled molecules pass through bare, so `TP53` stays `17`.
- Scaffolds with no assigned molecule (`Assigned-Molecule` = `na`) become
  `chrUn(<accession>)`.
- Names absent from the report are returned unchanged.

### The subtlety that makes this work

Ensembl names scaffolds **inconsistently**:

- unlocalized and unplaced scaffolds by GenBank accession — `GL000194.1`
- alt-scaffolds and patches by sequence name — `HG1_PATCH`,
  `HSCHR6_MHC_APD_CTG1`

Indexing *every* alias the report offers (`Sequence-Name`, `GenBank-Accn`,
`RefSeq-Accn`, `UCSC-style-name`) resolves all **503** non-standard human names
with zero unmatched. Keying on accession alone would have missed 458 of them;
keying on sequence name alone would have missed 45 — including MAFIP's.

Examples:

| Ensembl name | Label | Role |
| --- | --- | --- |
| `GL000194.1` | `chr14(GL000194.1)` | unlocalized-scaffold |
| `GL000195.1` | `chrUn(GL000195.1)` | unplaced-scaffold |
| `HG1_PATCH` | `chr14(KZ208920.1)` | fix-patch |
| `HSCHR6_MHC_APD_CTG1` | `chr6(GL000250.2)` | alt-scaffold |

Note the label always carries the **GenBank accession**, never the Ensembl
sequence name, so it is stable and self-describing.

## Gene counts

Human: **78,733 -> 86,411** genes.

| Sequence role | Genes added |
| --- | --- |
| unlocalized-scaffold | 138 |
| unplaced-scaffold | 70 |
| alt-scaffold | 4,847 |
| fix-patch | 1,828 |
| novel-patch | 795 |
| **total** | **7,678** |

Mouse gains 59 scaffold genes (78,348 total), rat gains 29 (43,360 total).
The merged human annotation is 86,451 rows.

## `chromosome_scope`

Alt-scaffold and patch genes are haplotype/patch copies that share gene symbols
with their primary counterparts, so `mapSymbolToID("HLA-DRB4")` now legitimately
returns three Ensembl IDs. To keep that controllable, `getAnnotationBiomart()`
gained a `chromosome_scope` argument:

| Scope | Keeps | Human rows | MAFIP |
| --- | --- | --- | --- |
| `"all"` (default) | everything | 86,411 | yes |
| `"primary"` | assembled molecules + unlocalized + unplaced | 78,941 | yes |
| `"standard"` | assembled molecules only | 78,733 | no |

`"standard"` reproduces the pre-change behaviour exactly. Filtering is applied on
the way out and the cache always stores the full table, so one cache serves every
scope and narrowing costs no extra download.

The merge builders are not scope-aware: they always use `"all"`, which is the
requested behaviour. Downstream filtering is easy because the `chromosome` column
is now self-describing.

## Bug fixed along the way

`buildHumanAnnotation()`, `buildMouseAnnotation()` and `buildRatAnnotation()`
accepted `annotation_dir` but never passed it to their source functions, which
silently fell back to the default user data directory. The visible symptom was
the merged log recording:

```json
"source_details": { "NCBI": "Log file not found", "HGNC": "Log file not found" }
```

Fixed in `R/mergeAnnotations.R`; the source logs now resolve.

## Files touched

- `R/chromosomes.R` — new: assembly report fetch/parse/cache, scaffold lookup,
  `.format_chromosome()`, `.chromosome_role()`, `.filter_chromosome_scope()`
- `R/getAnnotation.R` — BioMart filter removed + `chromosome_scope`; HGNC source
  switched; MGI and RGD chromosome columns labelled
- `R/mergeAnnotations.R` — `annotation_dir` passed through to source functions
- `inst/extdata/Assembly_Reports/` — new bundled offline fallback (~100 KB)
- `inst/extdata/Ensembl_Genes/`, `inst/extdata/HGNC_human_Genes/` — regenerated
  so offline users get the scaffold-aware tables
- `DESCRIPTION` — `Config/testthat/edition: 3`; `curl` and `withr` suggested

## Verification

43 tests pass, no failures, no skips.

- `tests/testthat/test-chromosomes.R` — offline unit tests against a checked-in
  assembly-report excerpt (`tests/testthat/fixtures/`). Covers parsing, alias
  resolution, `chrUn`, passthrough of assembled molecules and unknown names,
  the NULL-lookup no-op, role recovery and scope filtering.
- `tests/testthat/test-scaffold-genes.R` — live tests against Ensembl, NCBI and
  HGNC (`skip_on_cran()`, `skip_if_offline()`), sharing one download directory.

The headline case, through the full merge:

```
ensembl_id       gene_symbol  chromosome         ncbi_id  hgnc_id
ENSG00000274847  MAFIP        chr14(GL000194.1)  727764   31102
```

And through `mapSymbolToID(..., multi_handling = "collapse")`:

```
query_gene_symbol  gene_symbol  ensembl_id                                       chromosome
MAFIP              MAFIP        ENSG00000274847                                  chr14(GL000194.1)
TP53               TP53         ENSG00000141510                                  17
HLA-DRB4           HLA-DRB4     ENSG00000227357|ENSG00000227826|ENSG00000231021  chr6(GL000256.2)|chr6(GL000253.2)|chr6(GL000254.2)
```

Also confirmed with the network stubbed out and an empty cache:
`getAnnotationBiomart("human")` serves 86,411 rows from the bundled data with
labels intact, and `chromosome_scope = "primary"` still filters correctly using
the bundled assembly reports.

The merged output schema is unchanged. The 79 duplicate `ensembl_id` rows in the
merged human table are pre-existing NCBI symbol-rescue artifacts on standard
chromosomes, not introduced by this change.
