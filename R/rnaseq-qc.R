#' Summarize RNA-seq sample quality
#'
#' Calculates library-size, detected-feature, sample-correlation,
#' sample-distance, and principal-component summaries from a per-dataset RNA-seq
#' experiment returned by [GenerateRawSE()]. This is an exploratory QC function;
#' it does not change counts or fit a differential-expression model.
#'
#' @param se A `SummarizedExperiment` containing a `counts` assay. A
#'   `transformed` assay is used for sample relationships when present;
#'   otherwise `log2(normalized_counts + 1)` or `log2(counts + 1)` is used.
#' @param min_count The minimum count used to identify a detected feature.
#'
#' @return A list containing `sample_metrics`, `correlation`, `distance`,
#'   `pca`, `library_plot`, and `pca_plot`.
#' @examples
#' counts <- matrix(
#'   c(10, 12, 30, 28, 50, 48, 8, 9, 20, 22, 60, 58),
#'   nrow = 3,
#'   dimnames = list(paste0("gene", 1:3), paste0("sample", 1:4))
#' )
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   assays = list(counts = counts),
#'   colData = data.frame(
#'     condition = rep(c("control", "treated"), each = 2),
#'     row.names = colnames(counts)
#'   )
#' )
#' qc <- RNAseqQC(se)
#' qc$sample_metrics
#' @family transcriptomics functions
#' @importFrom ggplot2 aes geom_col geom_point ggplot labs theme_bw
#' @export
RNAseqQC <- function(se, min_count = 10) {
  if (!inherits(se, "SummarizedExperiment")) {
    stop("`se` must inherit from SummarizedExperiment.", call. = FALSE)
  }
  if (!"counts" %in% SummarizedExperiment::assayNames(se)) {
    stop("`se` must contain a `counts` assay.", call. = FALSE)
  }
  if (!is.numeric(min_count) || length(min_count) != 1L ||
      !is.finite(min_count) || min_count < 0) {
    stop("`min_count` must be one finite, non-negative number.", call. = FALSE)
  }
  counts_matrix <- as.matrix(SummarizedExperiment::assay(se, "counts"))
  if (ncol(counts_matrix) < 2L) {
    stop("At least two samples are required for RNA-seq QC.", call. = FALSE)
  }
  transformed <- if ("transformed" %in% SummarizedExperiment::assayNames(se)) {
    as.matrix(SummarizedExperiment::assay(se, "transformed"))
  } else if ("normalized_counts" %in% SummarizedExperiment::assayNames(se)) {
    log2(as.matrix(SummarizedExperiment::assay(se, "normalized_counts")) + 1)
  } else {
    log2(counts_matrix + 1)
  }
  sample_metrics <- data.frame(
    sample = colnames(counts_matrix),
    library_size = colSums(counts_matrix, na.rm = TRUE),
    detected_features = colSums(counts_matrix >= min_count, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  sample_data <- as.data.frame(SummarizedExperiment::colData(se))
  if (nrow(sample_data)) {
    sample_data$sample <- rownames(sample_data)
    sample_metrics <- merge(sample_metrics, sample_data, by = "sample", all.x = TRUE, sort = FALSE)
    sample_metrics <- sample_metrics[match(colnames(counts_matrix), sample_metrics$sample), , drop = FALSE]
  }
  correlation <- stats::cor(transformed, use = "pairwise.complete.obs", method = "spearman")
  distance <- as.matrix(stats::dist(t(transformed)))
  variable <- apply(transformed, 1L, stats::var, na.rm = TRUE)
  usable <- is.finite(variable) & variable > 0
  if (sum(usable) < 2L) {
    stop("At least two variable features are required for PCA.", call. = FALSE)
  }
  pca_fit <- stats::prcomp(t(transformed[usable, , drop = FALSE]), scale. = FALSE)
  pca <- as.data.frame(pca_fit$x[, seq_len(min(3L, ncol(pca_fit$x))), drop = FALSE])
  pca$sample <- rownames(pca)
  if (nrow(sample_data)) {
    pca <- merge(pca, sample_data, by = "sample", all.x = TRUE, sort = FALSE)
  }
  pca <- pca[match(colnames(counts_matrix), pca$sample), , drop = FALSE]
  group_column <- intersect(c("condition", "dataset"), names(pca))[1]
  if (is.na(group_column)) {
    pca$group <- "all samples"
  } else {
    pca$group <- as.character(pca[[group_column]])
  }
  variance_explained <- pca_fit$sdev^2 / sum(pca_fit$sdev^2)
  library_plot <- ggplot2::ggplot(
    sample_metrics,
    ggplot2::aes(x = stats::reorder(sample, library_size), y = library_size)
  ) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Sample", y = "Library size") +
    ggplot2::theme_bw()
  pca_plot <- ggplot2::ggplot(
    pca,
    ggplot2::aes(x = PC1, y = PC2, color = group)
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::labs(
      x = sprintf("PC1 (%.1f%%)", 100 * variance_explained[1]),
      y = sprintf("PC2 (%.1f%%)", 100 * variance_explained[2]),
      color = if (is.na(group_column)) "Group" else group_column
    ) +
    ggplot2::theme_bw()
  list(
    sample_metrics = sample_metrics,
    correlation = correlation,
    distance = distance,
    pca = pca,
    library_plot = library_plot,
    pca_plot = pca_plot
  )
}
