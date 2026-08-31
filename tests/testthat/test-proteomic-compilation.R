make_proteomic_experiment <- function(symbols, diff, contrast = "A_vs_B") {
  values <- matrix(
    seq_along(symbols), nrow = length(symbols), ncol = 2L,
    dimnames = list(symbols, c("s1", "s2"))
  )
  row_data <- data.frame(
    name = symbols,
    ID = symbols,
    GeneSymbol = symbols,
    stringsAsFactors = FALSE
  )
  row_data[[paste0(contrast, "_diff")]] <- diff
  row_data[[paste0(contrast, "_p.val")]] <- rep(0.01, length(symbols))
  row_data[[paste0(contrast, "_p.adj")]] <- rep(0.02, length(symbols))
  row_data[[paste0(contrast, "_BHCorrection")]] <- rep(0.03, length(symbols))
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = values), rowData = row_data
  )
}

test_that("a proteomics database is built from scratch", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  input <- withr::local_tempdir()
  output <- tempfile(fileext = ".rds")
  saveRDS(
    make_proteomic_experiment(c("GAPDH", "ACTB"), c(1.5, -2)),
    file.path(input, "PXD000001-A_vs_B.rds")
  )
  se <- suppressMessages(DESEProtGenerate(input, output))
  expect_named(SummarizedExperiment::assays(se), "PXD000001-A_vs_B")
  expect_equal(
    colnames(SummarizedExperiment::assays(se)[[1]]),
    c("logFC", "AdjPValue", "Pvalue", "BHCorrection")
  )
  assay_values <- SummarizedExperiment::assays(se)[[1]]
  expect_equal(assay_values["GAPDH", "logFC"], 1.5)
  expect_equal(assay_values["ACTB", "logFC"], -2)
  expect_true(file.exists(output))
})

test_that("an existing proteomics database is extended, not rebuilt", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  input <- withr::local_tempdir()
  output <- tempfile(fileext = ".rds")
  saveRDS(
    make_proteomic_experiment(c("GAPDH", "ACTB"), c(1.5, -2)),
    file.path(input, "PXD000001-A_vs_B.rds")
  )
  first <- suppressMessages(DESEProtGenerate(input, output))

  saveRDS(
    make_proteomic_experiment(c("GAPDH", "ACTB"), c(-1, 4), contrast = "C_vs_D"),
    file.path(input, "PXD000002-C_vs_D.rds")
  )
  second <- suppressMessages(DESEProtGenerate(input, output))

  expect_equal(
    names(SummarizedExperiment::assays(second)),
    c("PXD000001-A_vs_B", "PXD000002-C_vs_D")
  )
  expect_equal(nrow(second), nrow(first))
  expect_equal(
    SummarizedExperiment::assays(second)[["PXD000001-A_vs_B"]]["GAPDH", "logFC"],
    1.5
  )
  expect_equal(
    SummarizedExperiment::assays(second)[["PXD000002-C_vs_D"]]["ACTB", "logFC"],
    4
  )
})

test_that("compilation does not leak a Scratch variable into the caller", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  input <- withr::local_tempdir()
  saveRDS(
    make_proteomic_experiment("GAPDH", 1),
    file.path(input, "PXD000003-A_vs_B.rds")
  )
  expect_false(exists("Scratch", envir = globalenv(), inherits = FALSE))
  suppressMessages(DESEProtGenerate(input, tempfile(fileext = ".rds")))
  expect_false(exists("Scratch", envir = globalenv(), inherits = FALSE))
})

test_that("semicolon-separated protein groups are expanded", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  input <- withr::local_tempdir()
  saveRDS(
    make_proteomic_experiment("GAPDH;ACTB", 2.5),
    file.path(input, "PXD000004-A_vs_B.rds")
  )
  se <- suppressMessages(DESEProtGenerate(input, tempfile(fileext = ".rds")))
  assay_values <- SummarizedExperiment::assays(se)[[1]]
  expect_equal(assay_values["GAPDH", "logFC"], 2.5)
  expect_equal(assay_values["ACTB", "logFC"], 2.5)
})

test_that("files holding several comparisons are rejected with a clear message", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if_not_installed("org.Mm.eg.db")
  input <- withr::local_tempdir()
  experiment <- make_proteomic_experiment("GAPDH", 1)
  SummarizedExperiment::rowData(experiment)$C_vs_D_diff <- 2
  saveRDS(experiment, file.path(input, "PXD000005-A_vs_B.rds"))
  expect_error(
    suppressMessages(DESEProtGenerate(input, tempfile(fileext = ".rds"))),
    "Save one comparison per file"
  )
})
