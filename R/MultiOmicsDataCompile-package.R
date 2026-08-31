#' MultiOmicsDataCompile: compile and explore multi-omics data
#'
#' MultiOmicsDataCompile provides data-ingestion and harmonization workflows for
#' public
#' transcriptomic and proteomic datasets, utilities for constructing integrated
#' `SummarizedExperiment` objects, and lightweight exploratory visualizations.
#'
#' @keywords internal
#' @importFrom AnnotationDbi mappedkeys
#' @importFrom Biobase exprs "exprs<-" fvarLabels "fvarLabels<-" pData
#' @importFrom dplyr %>% arrange desc ends_with filter group_by left_join n select summarize
#' @importFrom ggplot2 aes element_text facet_wrap geom_boxplot geom_col
#'   geom_density geom_point geom_text geom_vline ggplot ggtitle labs
#'   scale_color_manual theme theme_bw xlab xlim ylim
#' @importFrom graphics boxplot par
#' @importFrom grDevices dev.off png
#' @importFrom grid unit
#' @importFrom methods as
#' @importFrom readxl read_excel read_xlsx
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom org.Mm.eg.db org.Mm.eg.db
#' @importFrom org.Rn.eg.db org.Rn.eg.db
#' @importFrom stringr str_split
#' @importFrom SummarizedExperiment SummarizedExperiment assay "assay<-" assays rowData "rowData<-" colData
#' @importFrom stats as.formula complete.cases cor lm p.adjust quantile rnorm sd
#' @importFrom textclean mgsub
"_PACKAGE"
