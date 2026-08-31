# The explorer's plotting helpers are defined in the packaged server script.
# They are loaded in isolation here so that the app's reshaping path is covered
# without starting a Shiny session or building a full database.
load_app_helpers <- function(names) {
  script <- system.file(
    "extdata", "Scripts", "server.R", package = "MultiOmicsDataCompile"
  )
  skip_if(!nzchar(script), "packaged server script is unavailable")
  env <- new.env(parent = globalenv())
  for (expression in as.list(parse(script))) {
    if (!is.call(expression)) next
    operator <- as.character(expression[[1]])
    if (length(operator) != 1L || !operator %in% c("<-", "=")) next
    target <- expression[[2]]
    if (!is.name(target)) next
    if (as.character(target) %in% names) eval(expression, env)
  }
  missing <- setdiff(names, ls(env, all.names = TRUE))
  if (length(missing)) {
    stop("Not found in server.R: ", paste(missing, collapse = ", "))
  }
  env
}

make_explorer_fixture <- function() {
  samples <- c("GSE1_A", "GSE1_B", "GSE2_A", "GSE2_B")
  values <- matrix(
    c(-1, 1, -0.5, 0.5, 1, -1, 0.5, -0.5),
    nrow = 2, dimnames = list(c("Human:1", "Human:2"), samples)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(Expression = values),
    rowData = data.frame(
      SYMBOL = c("GAPDH", "ACTB"), row.names = rownames(values),
      stringsAsFactors = FALSE
    ),
    colData = data.frame(
      dataset = c("GSE1", "GSE1", "GSE2", "GSE2"),
      rawcolumnnames = samples,
      condition = c("Ctrl", "NASH", "Ctrl", "NASH"),
      row.names = samples,
      stringsAsFactors = FALSE
    )
  )
}

test_that("the single-dataset violin plot builds without a mia dependency", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")
  withr::local_package("ggplot2")
  withr::local_package("RColorBrewer")
  withr::local_package("SummarizedExperiment")
  helpers <- load_app_helpers(c(
    "ViolinSampleColumns", "ExpressionLongFormat", "ViolinPalette",
    "box_violin_plot_v2.single"
  ))
  plot <- helpers$box_violin_plot_v2.single(
    make_explorer_fixture(), "GAPDH", "GSE1", FALSE
  )
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 2L)
  expect_setequal(plot$data$condition, c("Ctrl", "NASH"))
})

test_that("the multi-dataset violin plot facets by dataset", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("RColorBrewer")
  withr::local_package("ggplot2")
  withr::local_package("RColorBrewer")
  withr::local_package("SummarizedExperiment")
  helpers <- load_app_helpers(c(
    "ViolinSampleColumns", "ExpressionLongFormat", "ViolinPalette",
    "box_violin_plot_v2.multi"
  ))
  overview <- data.frame(
    ID = c("GSE1", "GSE2"),
    Disease = "NAFLD", Technology = "RNAseq", Tissue = "Liver",
    stringsAsFactors = FALSE
  )
  plot <- helpers$box_violin_plot_v2.multi(
    make_explorer_fixture(), overview, "GAPDH",
    "NAFLD", "RNAseq", "Liver", FALSE
  )
  expect_s3_class(plot, "ggplot")
  expect_true("facet_name" %in% names(plot$data))
  expect_setequal(plot$data$facet_name, c("GSE1", "GSE2"))
})

test_that("reshaping tolerates sample annotations the metadata lack", {
  skip_if_not_installed("SummarizedExperiment")
  withr::local_package("SummarizedExperiment")
  helpers <- load_app_helpers(c("ViolinSampleColumns", "ExpressionLongFormat"))
  # cell_type, tissue_type and cell_line are absent from this fixture.
  long <- helpers$ExpressionLongFormat(
    make_explorer_fixture(), "ACTB", c("GSE1_A", "GSE1_B")
  )
  expect_equal(nrow(long), 2L)
  expect_true(all(c("SYMBOL", "dataset", "condition") %in% names(long)))
})

test_that("check_na reports genes that are absent from a dataset", {
  skip_if_not_installed("SummarizedExperiment")
  withr::local_package("SummarizedExperiment")
  helpers <- load_app_helpers("check_na")
  fixture <- make_explorer_fixture()
  expect_false(helpers$check_na(fixture, "GAPDH", "GSE1"))
  expect_true(helpers$check_na(fixture, "NOT_A_GENE", "GSE1"))
})
