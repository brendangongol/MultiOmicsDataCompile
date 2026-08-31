make_direction_fixture <- function(directory) {
  # Half the expressed features go up and half go down, so median-of-ratios
  # normalization cannot absorb the change as a library-size shift.
  control <- c(rep(1000, 6), 2000, 3000, 4000, 2000, 3000, 4000)
  treated <- c(2000, 3000, 4000, 2000, 3000, 4000, rep(1000, 6))
  expressed <- cbind(control, control, treated, treated)

  # Features that fail the count filter are placed on both sides of the
  # expressed block: any code that pairs the filtered fold changes with
  # unfiltered annotations by position picks up the wrong gene symbols.
  quiet <- matrix(0, nrow = 20L, ncol = 4L)
  counts <- rbind(quiet, expressed, quiet)
  colnames(counts) <- paste0("GSE1_S", 1:4)
  symbols <- c(
    paste0("QUIET_A", seq_len(20)),
    paste0("EXPRESSED", seq_len(12)),
    paste0("QUIET_B", seq_len(20))
  )
  raw <- data.table::data.table(
    ENTREZID = as.character(seq_along(symbols)),
    SYMBOL = symbols
  )
  raw <- cbind(raw, data.table::as.data.table(counts))
  saveRDS(raw, file.path(directory, "GSE1_Raw.rds"))

  de <- data.table::data.table(
    SYMBOL = symbols,
    ID = as.character(seq_along(symbols)),
    ENTREZID = as.character(seq_along(symbols)),
    GENENAME = symbols,
    `Tx-Ctrl_logFC` = c(rep(0, 20), log2(treated / control), rep(0, 20)),
    `Tx-Ctrl_Pvalue` = 0.01,
    `Tx-Ctrl_AdjPValue` = 0.02
  )
  saveRDS(de, file.path(directory, "GSE1_Tx_Ctrl.rds"))
  invisible(list(symbols = symbols, logfc = log2(treated / control)))
}

test_that("RNA-seq direction checks stay aligned after count filtering", {
  skip_if_not_installed("DESeq2")
  db <- withr::local_tempdir()
  graphs <- withr::local_tempdir()
  qc <- withr::local_tempdir()
  make_direction_fixture(db)

  result <- suppressWarnings(GEO2RDirectionCheck(
    DBPath = db,
    DS = "GSE1",
    namestr = "Tx-Ctrl_logFC,Tx-Ctrl_Pvalue,Tx-Ctrl_AdjPValue",
    gsm = "0011",
    Technology = "RNAseq",
    GraphPath = graphs,
    subsetRaw = FALSE,
    writeRaw = TRUE,
    RawQCPath = qc
  ))

  expect_equal(nrow(result$checkFile), 1L)
  # Only the expressed features survive the filter, and each one agrees in
  # direction with the stored fold change.
  expect_equal(result$checkFile$correct, 12L)
  expect_equal(result$checkFile$incorrect, 0L)
  expect_gt(result$checkFile$correlation, 0.99)
})

test_that("array direction checks still use every stored feature", {
  db <- withr::local_tempdir()
  graphs <- withr::local_tempdir()
  qc <- withr::local_tempdir()
  symbols <- paste0("GENE", 1:6)
  control <- c(1, 2, 3, 4, 5, 6)
  treated <- c(2, 6, 12, 8, 10, 24)
  raw <- data.table::data.table(
    ENTREZID = as.character(1:6),
    SYMBOL = symbols,
    GSE2_S1 = control, GSE2_S2 = control,
    GSE2_S3 = treated, GSE2_S4 = treated
  )
  saveRDS(raw, file.path(db, "GSE2_Raw.rds"))
  # These values stay below the automatic log2 threshold, so the recomputed
  # summary is a difference on the stored scale.
  de <- data.table::data.table(
    SYMBOL = symbols,
    ID = as.character(1:6),
    `Tx-Ctrl_logFC` = treated - control,
    `Tx-Ctrl_Pvalue` = 0.01,
    `Tx-Ctrl_AdjPValue` = 0.02
  )
  saveRDS(de, file.path(db, "GSE2_Tx_Ctrl.rds"))

  result <- GEO2RDirectionCheck(
    DBPath = db, DS = "GSE2",
    namestr = "Tx-Ctrl_logFC,Tx-Ctrl_Pvalue,Tx-Ctrl_AdjPValue",
    gsm = "0011", Technology = "Array", GraphPath = graphs,
    subsetRaw = FALSE, writeRaw = TRUE, RawQCPath = qc
  )
  expect_equal(result$checkFile$correct, 6L)
  expect_equal(result$checkFile$incorrect, 0L)
})

test_that("an unsupported technology is rejected rather than silently reused", {
  db <- withr::local_tempdir()
  make_direction_fixture(db)
  expect_error(
    GEO2RDirectionCheck(
      DBPath = db, DS = "GSE1",
      namestr = "Tx-Ctrl_logFC,Tx-Ctrl_Pvalue,Tx-Ctrl_AdjPValue",
      gsm = "0011", Technology = "Proteomics",
      GraphPath = withr::local_tempdir(), subsetRaw = FALSE,
      writeRaw = TRUE, RawQCPath = withr::local_tempdir()
    ),
    "Unsupported Technology"
  )
})
