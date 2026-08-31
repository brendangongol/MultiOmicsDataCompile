make_condition_fixture <- function(directory) {
  raw <- data.table::data.table(
    ENTREZID = c("1", "2"),
    SYMBOL = c("A", "B"),
    GSE1_S1 = c(5, 3), GSE1_S2 = c(6, 4),
    GSE1_S3 = c(9, 1), GSE1_S4 = c(8, 2)
  )
  saveRDS(raw, file.path(directory, "GSE1_Raw.rds"))
  data.frame(
    ID = "GSE1",
    GEO2R = TRUE,
    ComparisonVector = "0011",
    FCColumnNames = "NASH-Ctrl_logFC,NASH-Ctrl_Pvalue,NASH-Ctrl_AdjPValue",
    stringsAsFactors = FALSE
  )
}

test_that("conditions are read from the manifest comparison code", {
  directory <- withr::local_tempdir()
  manifest <- make_condition_fixture(directory)
  conditions <- SampleConditions(directory, manifest)
  expect_equal(nrow(conditions), 4L)
  expect_equal(
    conditions$condition[match(paste0("GSE1_S", 1:4), conditions$rawcolumnnames)],
    c("Ctrl", "Ctrl", "NASH", "NASH")
  )
})

test_that("samples excluded from every comparison are labelled, not guessed", {
  directory <- withr::local_tempdir()
  manifest <- make_condition_fixture(directory)
  manifest$ComparisonVector <- "0X1X"
  conditions <- SampleConditions(directory, manifest)
  labelled <- conditions$condition[match(paste0("GSE1_S", 1:4), conditions$rawcolumnnames)]
  expect_equal(labelled, c("Ctrl", "Unassigned", "NASH", "Unassigned"))
})

test_that("the first comparison that includes a sample names its condition", {
  directory <- withr::local_tempdir()
  manifest <- make_condition_fixture(directory)
  manifest <- rbind(manifest, manifest)
  manifest$ComparisonVector <- c("00XX", "XX01")
  manifest$FCColumnNames <- c(
    "NASH-Ctrl_logFC,NASH-Ctrl_Pvalue,NASH-Ctrl_AdjPValue",
    "Fibrosis-Mild_logFC,Fibrosis-Mild_Pvalue,Fibrosis-Mild_AdjPValue"
  )
  conditions <- SampleConditions(directory, manifest)
  expect_equal(
    conditions$condition[match(paste0("GSE1_S", 1:4), conditions$rawcolumnnames)],
    c("Ctrl", "Ctrl", "Mild", "Fibrosis")
  )
})

test_that("a comparison code that disagrees with the stored samples is rejected", {
  directory <- withr::local_tempdir()
  manifest <- make_condition_fixture(directory)
  manifest$ComparisonVector <- "001"
  expect_error(SampleConditions(directory, manifest), "3 entries but")
})

test_that("MetDataCompile attaches supplied conditions", {
  rnaseq_dir <- withr::local_tempdir()
  array_dir <- withr::local_tempdir()
  data.table::fwrite(
    data.table::data.table(
      Run = c("SRR1", "SRR2"),
      `Sample Name` = c("GSM1", "GSM2"),
      `Assay Type` = "RNA-Seq"
    ),
    file.path(rnaseq_dir, "GSE1_SraRunTable.txt")
  )
  overview <- data.frame(
    ID = "GSE1", Tissue = "Liver", Disease = "NAFLD",
    Species = "Human", Technology = "RNAseq",
    stringsAsFactors = FALSE
  )
  conditions <- data.frame(
    rawcolumnnames = c("GSE1_GSM1", "GSE1_GSM2"),
    condition = c("Ctrl", "NASH"),
    stringsAsFactors = FALSE
  )
  result <- MetDataCompile(rnaseq_dir, array_dir, overview, conditions = conditions)
  expect_true("condition" %in% names(result))
  expect_equal(result[c("GSE1_GSM1", "GSE1_GSM2"), "condition"], c("Ctrl", "NASH"))
})

test_that("MetDataCompile without conditions leaves the column absent", {
  rnaseq_dir <- withr::local_tempdir()
  array_dir <- withr::local_tempdir()
  data.table::fwrite(
    data.table::data.table(
      Run = c("SRR1", "SRR2"),
      `Sample Name` = c("GSM1", "GSM2"),
      `Assay Type` = "RNA-Seq"
    ),
    file.path(rnaseq_dir, "GSE1_SraRunTable.txt")
  )
  overview <- data.frame(
    ID = "GSE1", Tissue = "Liver", Disease = "NAFLD",
    stringsAsFactors = FALSE
  )
  result <- expect_message(
    MetDataCompile(rnaseq_dir, array_dir, overview),
    "SampleConditions"
  )
  expect_false("condition" %in% names(result))
})
