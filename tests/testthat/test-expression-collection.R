make_expression_fixture <- function(missing_condition = FALSE) {
  samples <- c("GSE1_A", "GSE1_B", "GSE1_C", "GSE1_D", "GSE2_A", "GSE2_B")
  compiled <- data.table::data.table(
    FEATURE_ID = paste0("Human:", 1:6),
    SPECIES = "Human",
    ENTREZID = as.character(1:6),
    SYMBOL = paste0("GENE", 1:6),
    GSE1_A = c(10, 20, 30, 40, 50, 60),
    GSE1_B = c(12, 18, 35, 44, 48, 61),
    GSE1_C = c(30, 40, 20, 55, 70, 80),
    GSE1_D = c(28, 44, 19, 51, 72, 84),
    GSE2_A = c(5.1, 6.2, 7.3, 8.4, 9.5, 10.6),
    GSE2_B = c(5.4, 6.0, 7.8, 8.1, 9.8, 10.2)
  )
  metadata <- data.frame(
    dataset = c(rep("GSE1", 4), rep("GSE2", 2)),
    assay.type = c(rep("RNA-Seq", 4), rep("Array", 2)),
    rawcolumnnames = samples,
    condition = c("control", "control", "treated", "treated", "A", "B"),
    row.names = samples,
    stringsAsFactors = FALSE
  )
  if (missing_condition) metadata$condition[1] <- NA_character_
  overview <- data.frame(
    ID = c("GSE1", "GSE2"), Species = "Human",
    Technology = c("RNAseq", "Array"),
    stringsAsFactors = FALSE
  )
  list(compiled = compiled, metadata = metadata, overview = overview)
}

test_that("RNA-seq and array measurements remain separate", {
  fixture <- make_expression_fixture()
  result <- GenerateRawSE(fixture$metadata, fixture$compiled, fixture$overview)
  expect_named(result$RNAseqSE, "GSE1")
  expect_named(result$ArraySE, "GSE2")
  expect_equal(
    SummarizedExperiment::assayNames(result$RNAseqSE$GSE1),
    c("counts", "normalized_counts", "transformed")
  )
  expect_equal(
    SummarizedExperiment::assayNames(result$ArraySE$GSE2),
    "expression"
  )
  expect_equal(
    SummarizedExperiment::assayNames(result$ExplorationSE),
    "Expression"
  )
  expect_identical(result$RPKMSE, result$ExplorationSE)
})

test_that("missing RNA-seq conditions are never fabricated", {
  fixture <- make_expression_fixture(missing_condition = TRUE)
  expect_error(
    GenerateRawSE(fixture$metadata, fixture$compiled, fixture$overview),
    "never inferred"
  )
})

test_that("samples without metadata name the dataset responsible", {
  fixture <- make_expression_fixture()
  # An externally harmonized dataset now reaches this stage, so it needs
  # metadata like any other.
  fixture$compiled[, PID1_S1 := c(3, 4, 5, 6, 7, 8)]
  expect_error(
    GenerateRawSE(fixture$metadata, fixture$compiled, fixture$overview),
    "Affected dataset\\(s\\): PID1"
  )
})

test_that("RNA-seq QC returns sample summaries", {
  fixture <- make_expression_fixture()
  result <- GenerateRawSE(fixture$metadata, fixture$compiled, fixture$overview)
  qc <- RNAseqQC(result$RNAseqSE$GSE1)
  expect_equal(nrow(qc$sample_metrics), 4L)
  expect_equal(dim(qc$correlation), c(4L, 4L))
  expect_s3_class(qc$library_plot, "ggplot")
  expect_s3_class(qc$pca_plot, "ggplot")
})
