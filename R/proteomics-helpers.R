# Internal helpers for compiling proteomics differential-abundance databases.

.PROTEOMIC_ASSAY_COLUMNS <- c("logFC", "AdjPValue", "Pvalue", "BHCorrection")

# Build the uppercase gene-symbol annotation shared by every proteomics assay.
.proteomic_symbol_annotation <- function() {
  human <- data.table::as.data.table(suppressMessages(
    AnnotationDbi::select(
      org.Hs.eg.db,
      mappedkeys(org.Hs.eg.db::org.Hs.egSYMBOL),
      "SYMBOL", "ENTREZID"
    )
  ))
  data.table::setnames(human, "ENTREZID", "ENTREZID_Human")
  human <- human[!is.na(human$SYMBOL), ]

  mouse <- data.table::as.data.table(suppressMessages(
    AnnotationDbi::select(
      org.Mm.eg.db,
      mappedkeys(org.Mm.eg.db::org.Mm.egSYMBOL),
      "SYMBOL", "ENTREZID"
    )
  ))
  data.table::setnames(mouse, "ENTREZID", "ENTREZID_Mouse")
  mouse <- mouse[!is.na(mouse$SYMBOL), ]

  human[, SYMBOL := toupper(SYMBOL)]
  mouse[, SYMBOL := toupper(SYMBOL)]
  merged <- merge(human, mouse, by = "SYMBOL", all = TRUE)
  merged <- merged[!grepl("RIK$|RIK[0-9]$|---|1-DEC|1-MAR", merged$SYMBOL), ]

  # Keep one row per symbol, preferring the row that carries the most
  # cross-species identifiers.
  merged[, annotation_completeness :=
           (!is.na(ENTREZID_Human)) + (!is.na(ENTREZID_Mouse))]
  data.table::setorderv(merged, c("SYMBOL", "annotation_completeness"), c(1L, -1L))
  merged <- unique(merged, by = "SYMBOL")
  merged[, annotation_completeness := NULL]
  data.table::setorderv(merged, "SYMBOL")
  merged[]
}

.proteomic_statistic_column <- function(columns, pattern, path, description) {
  found <- grep(pattern, columns, value = TRUE)
  if (length(found) == 0L) {
    stop(
      basename(path), " does not contain a ", description, " column. ",
      "Run DEAnalysis() before compiling the database.",
      call. = FALSE
    )
  }
  if (length(found) > 1L) {
    stop(
      basename(path), " contains ", length(found), " ", description,
      " columns (", paste(found, collapse = ", "),
      "). Save one comparison per file.",
      call. = FALSE
    )
  }
  found
}

# Read one saved differential-abundance experiment and return its statistics as
# a data frame aligned to `feature_symbols`.
.proteomic_assay <- function(path, feature_symbols) {
  row_data <- data.table::as.data.table(
    SummarizedExperiment::rowData(readRDS(path))
  )
  columns <- colnames(row_data)
  symbol_column <- .proteomic_statistic_column(
    columns, "GeneSymbol", path, "gene-symbol"
  )
  selected <- c(
    symbol_column,
    .proteomic_statistic_column(columns, "_diff$", path, "log fold-change"),
    .proteomic_statistic_column(columns, "_p[.]adj$", path, "adjusted p-value"),
    .proteomic_statistic_column(columns, "_p[.]val$", path, "p-value"),
    .proteomic_statistic_column(columns, "_BHCorrection$", path, "BH-corrected p-value")
  )
  values <- row_data[, selected, with = FALSE]
  data.table::setnames(values, selected, c("SYMBOL", .PROTEOMIC_ASSAY_COLUMNS))
  values <- values[!is.na(values$SYMBOL) & nzchar(values$SYMBOL), ]
  values[, SYMBOL := toupper(SYMBOL)]

  # Protein groups list several symbols; expand them without growing a list in
  # a loop.
  parts <- strsplit(values$SYMBOL, ";", fixed = TRUE)
  counts <- lengths(parts)
  expanded <- data.table::data.table(
    SYMBOL = unlist(parts, use.names = FALSE),
    logFC = rep(values$logFC, counts),
    AdjPValue = rep(values$AdjPValue, counts),
    Pvalue = rep(values$Pvalue, counts),
    BHCorrection = rep(values$BHCorrection, counts)
  )
  expanded <- expanded[expanded$SYMBOL %in% feature_symbols, ]
  # Where a symbol appears more than once, keep the most significant record.
  expanded <- expanded[, .SD[which.min(Pvalue)], by = "SYMBOL"]

  assay_values <- matrix(
    NA_real_,
    nrow = length(feature_symbols),
    ncol = length(.PROTEOMIC_ASSAY_COLUMNS),
    dimnames = list(feature_symbols, .PROTEOMIC_ASSAY_COLUMNS)
  )
  matched <- match(expanded$SYMBOL, feature_symbols)
  assay_values[matched, ] <- as.matrix(
    expanded[, .PROTEOMIC_ASSAY_COLUMNS, with = FALSE]
  )
  as.data.frame(assay_values)
}

.proteomic_row_data <- function(annotation) {
  row_data <- as.data.frame(annotation, stringsAsFactors = FALSE)
  if (!"SYMBOL" %in% names(row_data)) {
    stop("Proteomics row annotations must contain a `SYMBOL` column.", call. = FALSE)
  }
  rownames(row_data) <- row_data$SYMBOL
  row_data
}
