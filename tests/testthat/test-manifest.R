test_that("dataset manifests are normalized and validated", {
  manifest <- data.frame(
    ID = c("GSE1", "GSE2"),
    Technology = c("RNA-Seq", "microarray"),
    Species = c("Homo sapiens", "mouse"),
    GEO2R = c(TRUE, TRUE),
    ComparisonVector = c("0011", "01"),
    stringsAsFactors = FALSE
  )
  result <- ValidateDatasetManifest(manifest)
  expect_equal(result$Technology, c("RNAseq", "Array"))
  expect_equal(result$Species, c("Human", "Mouse"))
})

test_that("dataset manifests reject missing comparisons", {
  manifest <- data.frame(
    ID = "GSE1", Technology = "RNAseq", Species = "Human",
    GEO2R = TRUE, ComparisonVector = "",
    stringsAsFactors = FALSE
  )
  expect_error(ValidateDatasetManifest(manifest), "require a ComparisonVector")
})

test_that("comparison codes must align with samples", {
  expect_error(
    .parse_comparison_code("01", c("a", "b", "c"), "GSE1"),
    "3 samples"
  )
  expect_error(
    .parse_comparison_code("00X", c("a", "b", "c"), "GSE1"),
    "control.*treatment"
  )
})
