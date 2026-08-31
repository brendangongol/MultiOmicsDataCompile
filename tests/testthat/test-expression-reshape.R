make_reshape_fixture <- function() {
  values <- matrix(
    c(1, 2, 3, 4, 5, 6), nrow = 3,
    dimnames = list(paste0("Human:", 1:3), c("GSE1_A", "GSE1_B"))
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(Expression = values),
    rowData = data.frame(
      SYMBOL = c("A", "B", "C"), row.names = rownames(values),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      dataset = c("GSE1", "GSE1"),
      condition = c("control", "treated"),
      row.names = colnames(values),
      stringsAsFactors = FALSE
    )
  )
}

test_that("assays are reshaped to one row per feature and sample", {
  se <- make_reshape_fixture()
  long <- MeltExpressionAssay(se, "Expression", "SYMBOL", c("dataset", "condition"))
  expect_equal(nrow(long), 6L)
  expect_true(all(c("feature", "sample", "SYMBOL", "dataset", "condition",
                    "Expression") %in% names(long)))
  expect_equal(
    long$Expression[long$feature == "Human:2" & long$sample == "GSE1_B"],
    5
  )
  expect_equal(
    long$condition[long$sample == "GSE1_B"][1],
    "treated"
  )
})

test_that("absent annotation columns are ignored rather than failing", {
  se <- make_reshape_fixture()
  long <- MeltExpressionAssay(
    se, "Expression", c("SYMBOL", "NOT_A_COLUMN"), c("condition", "cell_line")
  )
  expect_false("NOT_A_COLUMN" %in% names(long))
  expect_false("cell_line" %in% names(long))
  expect_true("condition" %in% names(long))
})

test_that("an unknown assay name is reported", {
  se <- make_reshape_fixture()
  expect_error(MeltExpressionAssay(se, "counts"), "does not contain an assay")
})
