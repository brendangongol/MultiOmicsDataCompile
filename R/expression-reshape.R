#' Reshape an experiment assay to long format
#'
#' Converts one assay of a [SummarizedExperiment::SummarizedExperiment()] to a
#' long data frame with one row per feature and sample, optionally carrying
#' selected `rowData` and `colData` columns. This replaces the equivalent
#' helpers in single-cell packages so that the explorer does not need an extra
#' dependency, and so that the reshaping used by the application can be tested.
#'
#' @param se A `SummarizedExperiment` object.
#' @param assay_name The name or index of the assay to reshape. When a name is
#'   supplied it is also used as the value column name.
#' @param row_columns An optional character vector of `rowData` columns to
#'   attach. Columns that are absent are ignored.
#' @param col_columns An optional character vector of `colData` columns to
#'   attach. Columns that are absent are ignored.
#'
#' @return A data frame with `feature`, `sample`, the requested annotation
#'   columns, and one value column.
#' @examples
#' counts <- matrix(
#'   1:6, nrow = 3,
#'   dimnames = list(paste0("gene", 1:3), c("s1", "s2"))
#' )
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(Expression = counts),
#'   rowData = data.frame(SYMBOL = c("A", "B", "C"), row.names = rownames(counts)),
#'   colData = data.frame(condition = c("ctrl", "tx"), row.names = colnames(counts))
#' )
#' head(MeltExpressionAssay(se, "Expression", "SYMBOL", "condition"))
#' @family exploration functions
#' @export
MeltExpressionAssay <- function(se, assay_name = 1L, row_columns = NULL,
                                col_columns = NULL) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must inherit from SummarizedExperiment.", call. = FALSE)
  }
  if (is.character(assay_name) &&
      !assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("`se` does not contain an assay named `", assay_name, "`.", call. = FALSE)
  }
  values <- as.matrix(SummarizedExperiment::assay(se, assay_name))
  if (is.null(rownames(values))) rownames(values) <- as.character(seq_len(nrow(values)))
  if (is.null(colnames(values))) colnames(values) <- as.character(seq_len(ncol(values)))
  value_name <- if (is.character(assay_name)) assay_name else "value"
  long <- data.frame(
    feature = rep(rownames(values), times = ncol(values)),
    sample = rep(colnames(values), each = nrow(values)),
    stringsAsFactors = FALSE
  )
  row_index <- rep(seq_len(nrow(values)), times = ncol(values))
  col_index <- rep(seq_len(ncol(values)), each = nrow(values))
  if (length(row_columns)) {
    row_data <- as.data.frame(SummarizedExperiment::rowData(se))
    for (column in intersect(row_columns, names(row_data))) {
      long[[column]] <- row_data[[column]][row_index]
    }
  }
  if (length(col_columns)) {
    col_data <- as.data.frame(SummarizedExperiment::colData(se))
    for (column in intersect(col_columns, names(col_data))) {
      long[[column]] <- col_data[[column]][col_index]
    }
  }
  long[[value_name]] <- as.vector(values)
  long
}
