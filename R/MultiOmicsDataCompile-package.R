#' MultiOmicsDataCompile: compile and explore multi-omics data
#'
#' MultiOmicsDataCompile provides data-ingestion and harmonization workflows for
#' public
#' transcriptomic and proteomic datasets, utilities for constructing integrated
#' `SummarizedExperiment` objects, and lightweight exploratory visualizations.
#'
#' @keywords internal
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom DEP theme_DEP1
#' @importFrom graphics boxplot par
#' @importFrom grDevices dev.off png
#' @importFrom grid unit
#' @importFrom readxl read_excel read_xlsx
#' @importFrom stats as.formula complete.cases cor lm p.adjust quantile rnorm sd
#' @importFrom textclean mgsub
"_PACKAGE"
