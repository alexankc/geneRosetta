test_that("authority symbols and synonyms survive consolidation for every builder", {
  authorities <- c(hgnc = "hgnc", mgi = "mgi", rgd = "rgd")

  for (authority in authorities) {
    symbol_col <- paste0("symbol_", authority)
    synonym_col <- paste0("syn_", authority)
    orphans <- data.frame(
      ncbi_id = c("54302", "54302", NA_character_),
      symbol = c("Slc14a2_v3", "Slc14a2_v3", "unused"),
      synonyms = c("Slc14a1_v3|UT-A3", "UT-A3", "unused"),
      stringsAsFactors = FALSE
    )
    names(orphans)[2:3] <- c(symbol_col, synonym_col)

    result <- .collect_authority_orphan_synonyms(orphans, symbol_col, synonym_col)

    expect_equal(result$ncbi_id, "54302", info = authority)
    expect_equal(
      result$syn_authority_orphan,
      "Slc14a2_v3|Slc14a1_v3|UT-A3|Slc14a2_v3|UT-A3",
      info = authority
    )
  }
})

test_that("orphan synonym collection changes no primary merge fields", {
  retained <- data.frame(
    ensembl_id = c("ENSRNOG1", "ENSRNOG2"),
    ncbi_id = c("54302", "54302"),
    symbol_rgd = c("Slc14a2", "Slc14a2"),
    rgd_id = c("3689", "3689"),
    stringsAsFactors = FALSE
  )
  orphans <- data.frame(
    ncbi_id = "54302",
    symbol_rgd = "Slc14a2_v3",
    syn_rgd = "UT-A3",
    stringsAsFactors = FALSE
  )

  extras <- .collect_authority_orphan_synonyms(orphans, "symbol_rgd", "syn_rgd")
  result <- dplyr::left_join(retained, extras, by = "ncbi_id")

  expect_equal(nrow(result), nrow(retained))
  expect_identical(result$ensembl_id, retained$ensembl_id)
  expect_identical(result$symbol_rgd, retained$symbol_rgd)
  expect_identical(result$rgd_id, retained$rgd_id)
  expect_equal(result$syn_authority_orphan, rep("Slc14a2_v3|UT-A3", 2))
})

test_that("deprecated records are resolved by stable identity, not synonyms", {
  annotation <- data.frame(
    ensembl_id = c("ENSG_OLD", "ENSG_CURRENT", "ENSG_ALIAS_OLD", "ENSG_ONLY"),
    gene_symbol = c("OLD42", "CURRENT42", "OLD_ALIAS", "ONLY99"),
    gene_synonyms = c(NA_character_, "OLD42|OLD_ALIAS", NA_character_, NA_character_),
    ncbi_id = c("42", "42", "77", "99"),
    hgnc_id = c("10", "10", "20", "30"),
    gene_biotype = c(NA_character_, "protein-coding", "ncRNA", NA_character_),
    gene_description = c(NA_character_, "complete NCBI record", "distinct record", NA_character_),
    stringsAsFactors = FALSE
  )

  result <- .resolve_deprecated_records(annotation, "ENSG_CURRENT")

  expect_equal(result$gene_symbol, c("CURRENT42", "OLD_ALIAS", "ONLY99"))
  expect_equal(result$ensembl_id, c("ENSG_CURRENT", NA_character_, NA_character_))
})

test_that("mouse deprecated IDs use MGI identity and retain distinct NCBI genes", {
  annotation <- data.frame(
    ensembl_id = c(
      "ENSMUSG_ARHGAP_CURRENT", "ENSMUSG_ARHGAP_OLD",
      "ENSMUSG_KISS_CURRENT", "ENSMUSG_KISS_OLD", NA_character_
    ),
    gene_symbol = c("Arhgap26", "Arhgap26", "Gm28040", "Kiss1", "Gm12863"),
    gene_biotype = c("protein_coding", NA, "protein_coding", NA, NA),
    gene_description = c("Rho GTPase activating protein 26", NA,
                         "predicted gene, 28040", NA, NA),
    chromosome = c("18", NA, "1", NA, NA),
    gene_start = c(39126097L, NA, 133237565L, NA, NA),
    gene_end = c(39509337L, NA, 133257465L, NA, NA),
    strand = c(1L, NA, 1L, NA, NA),
    metadata_source = c("Ensembl", NA, "Ensembl", NA, NA),
    ncbi_id = c("71302", "71302", NA, "280287", "100038571"),
    mgi_id = c("1918552", "1918552", "5547776", "2663985", "3651212"),
    stringsAsFactors = FALSE
  )
  ncbi <- data.frame(
    ncbi_id = c("71302", "280287", "100038571"),
    gene_biotype_ncbi = c("protein-coding", "protein-coding", "ncRNA"),
    gene_description_ncbi = c(
      "Rho GTPase activating protein 26",
      "KiSS-1 metastasis-suppressor",
      "predicted gene 12863"
    ),
    chromosome_ncbi = c("18", "1", "4"),
    gene_start_ncbi = c(38734531L, 133249625L, 118811470L),
    gene_end_ncbi = c(39509338L, 133257460L, 118812116L),
    strand_ncbi = c(1L, 1L, -1L),
    stringsAsFactors = FALSE
  )

  result <- annotation %>%
    .backfill_ncbi_metadata_by_id(ncbi) %>%
    .resolve_deprecated_records(
      c("ENSMUSG_ARHGAP_CURRENT", "ENSMUSG_KISS_CURRENT"),
      authority_id = "mgi_id"
    )

  expect_equal(sum(result$gene_symbol == "Arhgap26"), 1L)
  expect_equal(result$ensembl_id[result$gene_symbol == "Arhgap26"],
               "ENSMUSG_ARHGAP_CURRENT")

  kiss <- result[result$gene_symbol == "Kiss1", ]
  expect_equal(kiss$ensembl_id, NA_character_)
  expect_equal(kiss$ncbi_id, "280287")
  expect_equal(kiss$mgi_id, "2663985")
  expect_equal(kiss$chromosome, "1")
  expect_equal(kiss$gene_description, "KiSS-1 metastasis-suppressor")
  expect_equal(kiss$metadata_source, "NCBI")

  gm <- result[!is.na(result$ncbi_id) & result$ncbi_id == "100038571", ]
  expect_equal(gm$ensembl_id, NA_character_)
  expect_equal(gm$chromosome, "4")
  expect_equal(gm$gene_start, 118811470L)
  expect_equal(gm$gene_description, "predicted gene 12863")
  expect_equal(gm$metadata_source, "NCBI")
})

test_that("rat obsolete Ensembl rows sharing the NCBI-RGD pair are discarded", {
  annotation <- data.frame(
    ensembl_id = c(
      "ENSRNOG_A2M_CURRENT", "ENSRNOG_A2M_OLD",
      "ENSRNOG_AANAT_CURRENT", "ENSRNOG_AANAT_OLD_1",
      "ENSRNOG_AANAT_OLD_2", "ENSRNOG_DISTINCT_OLD"
    ),
    gene_symbol = c("A2m", "A2m", "Aanat", "Aanat", "Aanat", "Distinct"),
    gene_biotype = c("protein_coding", NA, "protein_coding", NA, NA, NA),
    gene_description = c("alpha-2-macroglobulin", NA,
                         "aralkylamine N-acetyltransferase", NA, NA, NA),
    chromosome = c("4", NA, "10", NA, NA, NA),
    gene_start = c(156569860L, NA, 102326135L, NA, NA, NA),
    gene_end = c(156619868L, NA, 102330635L, NA, NA, NA),
    strand = c(1L, NA, 1L, NA, NA, NA),
    metadata_source = c("Ensembl", NA, "Ensembl", NA, NA, NA),
    ncbi_id = c("24153", "24153", "25120", "25120", "25120", "99999"),
    rgd_id = c("2004", "2004", "2006", "2006", "2006", "9999"),
    stringsAsFactors = FALSE
  )
  ncbi <- data.frame(
    ncbi_id = c("24153", "25120", "99999"),
    gene_biotype_ncbi = c("protein-coding", "protein-coding", "ncRNA"),
    gene_description_ncbi = c(
      "alpha-2-macroglobulin",
      "aralkylamine N-acetyltransferase",
      "distinct rat record"
    ),
    chromosome_ncbi = c("4", "10", "7"),
    gene_start_ncbi = c(156570163L, 102323647L, 700L),
    gene_end_ncbi = c(156619870L, 102330639L, 900L),
    strand_ncbi = c(1L, 1L, -1L),
    stringsAsFactors = FALSE
  )

  result <- annotation %>%
    .backfill_ncbi_metadata_by_id(ncbi) %>%
    .resolve_deprecated_records(
      c("ENSRNOG_A2M_CURRENT", "ENSRNOG_AANAT_CURRENT"),
      authority_id = "rgd_id"
    )

  expect_equal(sum(result$gene_symbol == "A2m"), 1L)
  expect_equal(sum(result$gene_symbol == "Aanat"), 1L)
  expect_equal(result$ensembl_id[result$gene_symbol == "A2m"],
               "ENSRNOG_A2M_CURRENT")
  expect_equal(result$ensembl_id[result$gene_symbol == "Aanat"],
               "ENSRNOG_AANAT_CURRENT")

  distinct <- result[result$gene_symbol == "Distinct", ]
  expect_equal(distinct$ensembl_id, NA_character_)
  expect_equal(distinct$ncbi_id, "99999")
  expect_equal(distinct$rgd_id, "9999")
  expect_equal(distinct$chromosome, "7")
  expect_equal(distinct$metadata_source, "NCBI")
})

test_that("a partial stable-ID match does not discard a distinct authority record", {
  annotation <- data.frame(
    ensembl_id = c("ENS_CURRENT", "ENS_OLD_NCBI", "ENS_OLD_AUTHORITY"),
    gene_symbol = c("CURRENT", "DISTINCT_NCBI", "DISTINCT_AUTHORITY"),
    gene_biotype = c("protein_coding", "ncRNA", "ncRNA"),
    gene_description = c("current", "different authority", "different NCBI"),
    ncbi_id = c("1", "1", "2"),
    rgd_id = c("10", "20", "10"),
    stringsAsFactors = FALSE
  )

  result <- .resolve_deprecated_records(
    annotation, "ENS_CURRENT", authority_id = "rgd_id"
  )

  expect_equal(nrow(result), 3L)
  expect_equal(
    result$ensembl_id,
    c("ENS_CURRENT", NA_character_, NA_character_)
  )
})

test_that("an Ensembl-less duplicate of a complete mouse record is discarded", {
  annotation <- data.frame(
    ensembl_id = c("ENSMUSG00000097205", NA_character_),
    gene_symbol = c("Gm17540", "Gm17540"),
    gene_biotype = c("lncRNA", "ncRNA"),
    gene_description = c(
      "predicted gene, 17540 [Source:MGI Symbol;Acc:MGI:4937174]",
      "predicted gene, 17540"
    ),
    chromosome = c("9", NA_character_),
    gene_start = c(40442451L, NA_integer_),
    gene_end = c(40444916L, NA_integer_),
    strand = c(1L, NA_integer_),
    metadata_source = c("Ensembl", NA_character_),
    ncbi_id = c("100503925", "100503925"),
    mgi_id = c("4937174", "4937174"),
    stringsAsFactors = FALSE
  )

  result <- .resolve_deprecated_records(
    annotation, "ENSMUSG00000097205", authority_id = "mgi_id"
  )

  expect_equal(nrow(result), 1L)
  expect_equal(result$ensembl_id, "ENSMUSG00000097205")
  expect_equal(result$metadata_source, "Ensembl")
})

test_that("final annotations exclude rows lacking Ensembl and NCBI IDs", {
  annotation <- data.frame(
    ensembl_id = c("ENSRNOG1", NA_character_, NA_character_),
    ncbi_id = c(NA_character_, "123", NA_character_),
    rgd_id = c("10", "20", "30"),
    gene_symbol = c("ENSEMBL_ONLY", "NCBI_ONLY", "RGD_ONLY"),
    stringsAsFactors = FALSE
  )

  result <- annotation %>%
    dplyr::filter(!(is.na(ensembl_id) & is.na(ncbi_id)))

  expect_equal(result$gene_symbol, c("ENSEMBL_ONLY", "NCBI_ONLY"))
})

test_that("NCBI gene types are harmonized to Ensembl-style biotypes", {
  ncbi_types <- c(
    "protein-coding", "ncRNA", "pseudo", "biological-region", "other",
    "unknown", "tRNA", "rRNA", "snRNA", "snoRNA", "scRNA", NA_character_
  )

  expect_equal(
    .harmonize_gene_biotype(ncbi_types),
    c(
      "protein_coding", "lncRNA", "pseudogene", "biological_region",
      "misc_RNA", NA_character_, "tRNA", "rRNA", "snRNA", "snoRNA",
      "scRNA", NA_character_
    )
  )
})

test_that("final annotations retain one instance of each exact row", {
  annotation <- data.frame(
    ensembl_id = c(NA_character_, NA_character_, "ENS_DISTINCT"),
    gene_symbol = c("Or6c202d", "Or6c202d", "Or6c202d"),
    gene_synonyms = c("Olr1845", "Olr1845", "Olr1845"),
    ncbi_id = c("404902", "404902", "404902"),
    rgd_id = c("1333847", "1333847", "1333847"),
    stringsAsFactors = FALSE
  )

  result <- dplyr::distinct(annotation)

  expect_equal(nrow(result), 2L)
  expect_equal(sum(is.na(result$ensembl_id)), 1L)
  expect_true("ENS_DISTINCT" %in% result$ensembl_id)
})

test_that("NCBI metadata is backfilled by exact NCBI ID", {
  annotation <- data.frame(
    ncbi_id = "145438",
    gene_biotype = NA_character_,
    gene_description = NA_character_,
    chromosome = NA_character_,
    gene_start = NA_integer_,
    gene_end = NA_integer_,
    strand = NA_integer_,
    metadata_source = NA_character_,
    stringsAsFactors = FALSE
  )
  ncbi <- data.frame(
    ncbi_id = "145438",
    gene_biotype_ncbi = "ncRNA",
    gene_description_ncbi = "FRMD6 antisense RNA 1",
    chromosome_ncbi = "14",
    gene_start_ncbi = 51649516L,
    gene_end_ncbi = 51651744L,
    strand_ncbi = -1L,
    stringsAsFactors = FALSE
  )

  result <- .backfill_ncbi_metadata_by_id(annotation, ncbi)

  expect_equal(result$gene_biotype, "ncRNA")
  expect_equal(result$gene_description, "FRMD6 antisense RNA 1")
  expect_equal(result$chromosome, "14")
  expect_equal(result$gene_start, 51649516L)
  expect_equal(result$metadata_source, "NCBI")
})

test_that("records missing both biotype and description are removable", {
  annotation <- data.frame(
    gene_symbol = c("DROP", "KEEP_TYPE", "KEEP_DESCRIPTION"),
    gene_biotype = c(NA_character_, "ncRNA", NA_character_),
    gene_description = c(NA_character_, NA_character_, "described"),
    stringsAsFactors = FALSE
  )

  result <- annotation %>%
    dplyr::filter(!(is.na(gene_biotype) & is.na(gene_description)))

  expect_equal(result$gene_symbol, c("KEEP_TYPE", "KEEP_DESCRIPTION"))
})
