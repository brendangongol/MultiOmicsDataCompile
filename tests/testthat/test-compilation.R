test_that("raw RNA-seq duplicate identifiers are summed", {
  input <- withr::local_tempdir()
  output <- tempfile(fileext = ".txt")
  raw <- data.table::data.table(
    ENTREZID = c("1", "1", "2"),
    SYMBOL = c("A", "A", "B"),
    GSE1_S1 = c(5, 7, 3),
    GSE1_S2 = c(2, 4, 9)
  )
  saveRDS(raw, file.path(input, "GSE1_Raw.rds"))
  overview <- data.frame(
    ID = "GSE1", Species = "Human", Technology = "RNAseq",
    stringsAsFactors = FALSE
  )
  result <- RawDataCompile(input, output, overview = overview)
  expect_equal(result[ENTREZID == "1", GSE1_S1], 12)
  expect_true(file.exists(output))
})

test_that("externally harmonized count files are compiled", {
  input <- withr::local_tempdir()
  output <- tempfile(fileext = ".txt")
  saveRDS(
    data.table::data.table(
      ENTREZID = c("1", "2"), SYMBOL = c("A", "B"),
      PID1_S1 = c(4, 8), PID1_S2 = c(6, 2)
    ),
    file.path(input, "PID1_CountCodingGenes_CountRaw.rds")
  )
  overview <- data.frame(
    ID = "PID1", Species = "Human", Technology = "RNAseq",
    stringsAsFactors = FALSE
  )
  result <- suppressMessages(RawDataCompile(input, output, overview = overview))
  expect_true(all(c("PID1_S1", "PID1_S2") %in% names(result)))
  expect_equal(result[ENTREZID == "2", PID1_S1], 8)
})

test_that("RPKM files are skipped with an explicit message", {
  input <- withr::local_tempdir()
  output <- tempfile(fileext = ".txt")
  saveRDS(
    data.table::data.table(
      ENTREZID = "1", SYMBOL = "A", GSE1_S1 = 5, GSE1_S2 = 6
    ),
    file.path(input, "GSE1_Raw.rds")
  )
  saveRDS(
    data.table::data.table(
      ENTREZID = "1", SYMBOL = "A", PID1_S1 = 1.25, PID1_S2 = 2.5
    ),
    file.path(input, "PID1_CodingGenes_RPKMRaw.rds")
  )
  overview <- data.frame(
    ID = c("GSE1", "PID1"), Species = "Human", Technology = "RNAseq",
    stringsAsFactors = FALSE
  )
  expect_message(
    result <- RawDataCompile(input, output, overview = overview),
    "RPKM"
  )
  expect_false("PID1_S1" %in% names(result))
})

test_that("differential-expression compilation uses stable feature keys", {
  input <- withr::local_tempdir()
  output <- tempfile(fileext = ".rds")
  make_result <- function(offset) {
    x <- data.table::data.table(
      SYMBOL = c("A", "B"), ID = c("1", "2"), ENTREZID = c("1", "2"),
      GENENAME = c("alpha", "beta"),
      test_logFC = c(1, -1) + offset,
      test_Pvalue = c(0.01, 0.02),
      test_AdjPValue = c(0.02, 0.04)
    )
    attr(x, "provenance") <- list(species = "Human")
    x
  }
  saveRDS(make_result(0), file.path(input, "GSE1_A-B.rds"))
  saveRDS(make_result(0.5), file.path(input, "GSE2_A-B.rds"))
  se <- DESEGenerate(input, output)
  expect_true(all(grepl("^Human:", rownames(se))))
  expect_equal(length(SummarizedExperiment::assays(se)), 2L)
  expect_true(file.exists(output))
})
