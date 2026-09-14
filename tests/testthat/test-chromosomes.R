fixture_report <- function() {
  .parse_assembly_report(test_path("fixtures", "grch38_excerpt_assembly_report.txt"))
}

test_that(".parse_assembly_report reads the ten report columns as character", {
  report <- fixture_report()

  expect_equal(nrow(report), 8L)
  expect_named(report, .ASSEMBLY_REPORT_COLS)
  expect_true(all(vapply(report, is.character, logical(1))))
  expect_equal(report$genbank_accn[report$sequence_name == "HG1_PATCH"], "KZ208920.1")
})

test_that(".parse_assembly_report rejects an all-comment file", {
  empty <- tempfile(fileext = ".txt")
  on.exit(unlink(empty), add = TRUE)
  writeLines(c("# header only", "# nothing else"), empty)

  expect_error(.parse_assembly_report(empty), "no records")
})

test_that(".format_chromosome labels an unlocalized scaffold with its molecule", {
  lookup <- .chromosome_lookup(fixture_report())

  # MAFIP lives on GL000194.1, an unlocalized scaffold assigned to chromosome 14.
  expect_equal(.format_chromosome("GL000194.1", lookup), "chr14(GL000194.1)")
})

test_that(".format_chromosome resolves scaffolds by every alias the report offers", {
  lookup <- .chromosome_lookup(fixture_report())

  # Ensembl names alt-scaffolds and patches by sequence name, not accession.
  expect_equal(.format_chromosome("HG1_PATCH", lookup), "chr14(KZ208920.1)")
  expect_equal(.format_chromosome("HSCHR6_MHC_APD_CTG1", lookup), "chr6(GL000250.2)")
  # ...but names unlocalized/unplaced scaffolds by accession.
  expect_equal(.format_chromosome("GL000250.2", lookup), "chr6(GL000250.2)")
  expect_equal(.format_chromosome("chr14_GL000194v1_random", lookup), "chr14(GL000194.1)")
})

test_that(".format_chromosome sends scaffolds with no assigned molecule to chrUn", {
  lookup <- .chromosome_lookup(fixture_report())

  expect_equal(.format_chromosome("GL000195.1", lookup), "chrUn(GL000195.1)")
})

test_that(".format_chromosome passes through assembled molecules and unknowns", {
  lookup <- .chromosome_lookup(fixture_report())

  input <- c("14", "X", "MT", "NOT_A_SCAFFOLD", NA, "")
  expect_equal(.format_chromosome(input, lookup), input)
})

test_that(".format_chromosome is a no-op when the assembly report is unavailable", {
  input <- c("14", "GL000194.1", NA)

  expect_equal(.format_chromosome(input, NULL), input)
})

test_that(".chromosome_role recovers the sequence role from a formatted label", {
  lookup <- .chromosome_lookup(fixture_report())

  labels <- c("14", "chr14(GL000194.1)", "chrUn(GL000195.1)",
              "chr14(KZ208920.1)", "chr6(GL000250.2)", "NOT_A_SCAFFOLD")

  expect_equal(
    .chromosome_role(labels, lookup),
    c("assembled-molecule", "unlocalized-scaffold", "unplaced-scaffold",
      "fix-patch", "alt-scaffold", "unknown")
  )
})

test_that(".filter_chromosome_scope keeps the roles each scope admits", {
  lookup <- .chromosome_lookup(fixture_report())

  df <- data.frame(
    gene = c("primary", "unlocalized", "unplaced", "patch", "alt"),
    chromosome = c("14", "chr14(GL000194.1)", "chrUn(GL000195.1)",
                   "chr14(KZ208920.1)", "chr6(GL000250.2)"),
    stringsAsFactors = FALSE
  )

  expect_equal(.filter_chromosome_scope(df, "all", lookup)$gene, df$gene)
  expect_equal(
    .filter_chromosome_scope(df, "primary", lookup)$gene,
    c("primary", "unlocalized", "unplaced")
  )
  expect_equal(.filter_chromosome_scope(df, "standard", lookup)$gene, "primary")
})

test_that(".assembly_report_url falls back to the pinned assembly without a build", {
  url <- .assembly_report_url("human", genome_build = NULL)

  expect_match(url, "GCA_000001405\\.29_GRCh38\\.p14_assembly_report\\.txt$")
})

test_that(".assembly_report_url rejects unsupported species", {
  expect_error(.assembly_report_url("zebrafish"), "Unsupported species")
})
