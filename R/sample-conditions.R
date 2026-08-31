# Helpers that turn manifest comparison codes into sample-level conditions.

.comparison_labels <- function(code, dataset) {
  if (length(code) != 1L || is.na(code) || !nzchar(trimws(as.character(code)))) {
    stop("ComparisonVector is missing for ", dataset, ".", call. = FALSE)
  }
  code <- gsub("[[:space:],;]", "", as.character(code))
  labels <- strsplit(toupper(code), "", fixed = TRUE)[[1]]
  invalid <- setdiff(unique(labels), c("0", "1", "X"))
  if (length(invalid)) {
    stop(
      "ComparisonVector for ", dataset,
      " contains invalid values: ", paste(invalid, collapse = ", "),
      ". Use only 0, 1, and X.",
      call. = FALSE
    )
  }
  labels
}

.comparison_tag <- function(namestr) {
  nam <- strsplit(gsub("[[:space:]]", "", as.character(namestr)), ",", fixed = TRUE)[[1]]
  if (!length(nam) || is.na(nam[1])) {
    stop("`namestr` must contain at least one output column name.", call. = FALSE)
  }
  gsub("-", "_", sub("_.+$", "", nam[1]))
}

.raw_file_name <- function(dataset, namestr, subsetRaw = FALSE) {
  if (isTRUE(subsetRaw)) {
    paste0(dataset, "_", .comparison_tag(namestr), "_Raw.rds")
  } else {
    paste0(dataset, "_Raw.rds")
  }
}

.de_file_name <- function(dataset, namestr) {
  paste0(dataset, "_", .comparison_tag(namestr), ".rds")
}

# Derive the control and treatment group labels from an FCColumnNames entry
# such as "NASH-Ctrl_logFC,NASH-Ctrl_Pvalue,NASH-Ctrl_AdjPValue".
.comparison_group_names <- function(namestr) {
  nam <- strsplit(gsub("[[:space:]]", "", as.character(namestr)), ",", fixed = TRUE)[[1]]
  if (!length(nam) || is.na(nam[1]) || !nzchar(nam[1])) {
    return(c(control = "Control", treatment = "Treatment"))
  }
  label <- sub("_[^_]*$", "", nam[1])
  if (!grepl("-", label, fixed = TRUE)) {
    return(c(control = "Control", treatment = if (nzchar(label)) label else "Treatment"))
  }
  c(control = sub("^[^-]*-", "", label), treatment = sub("-.*$", "", label))
}

.expression_annotation_columns <- c(
  "FEATURE_ID", "SPECIES", "ID", "ENTREZID", "SYMBOL", "GENENAME",
  "ENTREZID_Human", "ENTREZID_Mouse", "Length_Human", "Length_Mouse"
)

#' Derive sample conditions from a dataset manifest
#'
#' Reads the raw-expression files written by [GEOCompile()] and assigns each
#' sample the group name recorded in the manifest's `ComparisonVector` and
#' `FCColumnNames` columns. Conditions are therefore taken from an explicit,
#' user-authored comparison code; they are never inferred from sample order.
#'
#' A comparison code positions one character per sample of the stored raw
#' matrix, so the assignment is a direct lookup rather than a guess. When a
#' dataset appears in several manifest rows, the first comparison that includes
#' a sample determines its condition. Samples excluded (`X`) from every
#' comparison are labelled `"Unassigned"` so that downstream objects still carry
#' a condition for each column without implying group membership.
#'
#' @param DBPath A path containing the raw-expression RDS files written by
#'   [GEOCompile()].
#' @param overview A dataset manifest containing `ID`, `ComparisonVector`, and
#'   `FCColumnNames` columns, and optionally a logical `GEO2R` column.
#' @param subsetRaw A logical value indicating whether the stored raw data were
#'   subset to the samples used in each comparison.
#'
#' @return A data table with `dataset`, `rawcolumnnames`, `condition`,
#'   `comparison`, and `comparison_label` columns, containing one row per
#'   sample.
#' @examples
#' directory <- tempfile()
#' dir.create(directory)
#' raw <- data.frame(
#'   ENTREZID = c("1", "2"), SYMBOL = c("A", "B"),
#'   GSE1_S1 = c(5, 3), GSE1_S2 = c(6, 4),
#'   GSE1_S3 = c(9, 1), GSE1_S4 = c(8, 2)
#' )
#' saveRDS(raw, file.path(directory, "GSE1_Raw.rds"))
#' manifest <- data.frame(
#'   ID = "GSE1",
#'   ComparisonVector = "0011",
#'   FCColumnNames = "NASH-Ctrl_logFC,NASH-Ctrl_Pvalue,NASH-Ctrl_AdjPValue",
#'   stringsAsFactors = FALSE
#' )
#' SampleConditions(directory, manifest)
#' @family transcriptomics functions
#' @import data.table
#' @export
SampleConditions <- function(DBPath, overview, subsetRaw = FALSE) {
  overview <- as.data.frame(overview, stringsAsFactors = FALSE)
  required <- c("ID", "ComparisonVector", "FCColumnNames")
  missing <- setdiff(required, names(overview))
  if (length(missing)) {
    stop(
      "`overview` is missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  selected <- !is.na(overview$ComparisonVector) &
    nzchar(trimws(as.character(overview$ComparisonVector)))
  if ("GEO2R" %in% names(overview)) {
    selected <- selected & !is.na(overview$GEO2R) & as.logical(overview$GEO2R)
  }
  overview <- overview[selected, , drop = FALSE]
  if (!nrow(overview)) {
    stop("`overview` contains no rows with a ComparisonVector.", call. = FALSE)
  }
  assignments <- vector("list", nrow(overview))
  missing_files <- character()
  for (i in seq_len(nrow(overview))) {
    dataset <- as.character(overview$ID[i])
    namestr <- overview$FCColumnNames[i]
    path <- file.path(DBPath, .raw_file_name(dataset, namestr, subsetRaw))
    if (!file.exists(path)) {
      missing_files <- c(missing_files, basename(path))
      next
    }
    stored <- readRDS(path)
    sample_columns <- setdiff(names(stored), .expression_annotation_columns)
    labels <- .comparison_labels(overview$ComparisonVector[i], dataset)
    if (isTRUE(subsetRaw)) {
      labels <- labels[labels != "X"]
    }
    if (length(labels) != length(sample_columns)) {
      stop(
        "ComparisonVector for ", dataset, " has ", length(labels),
        " entries but ", basename(path), " stores ", length(sample_columns),
        " samples.",
        call. = FALSE
      )
    }
    groups <- .comparison_group_names(namestr)
    condition <- rep(NA_character_, length(labels))
    condition[labels == "1"] <- groups[["treatment"]]
    condition[labels == "0"] <- groups[["control"]]
    assignments[[i]] <- data.table::data.table(
      dataset = dataset,
      rawcolumnnames = sample_columns,
      comparison = .comparison_tag(namestr),
      comparison_label = labels,
      condition = condition
    )
  }
  if (length(missing_files)) {
    message(
      "No raw-expression file was found for ", length(missing_files),
      " comparison(s); skipping: ",
      paste(utils::head(missing_files, 5L), collapse = ", "),
      if (length(missing_files) > 5L) " ..." else ""
    )
  }
  assignments <- data.table::rbindlist(assignments, use.names = TRUE)
  if (!nrow(assignments)) {
    stop(
      "No raw-expression files matching the manifest were found in `DBPath`.",
      call. = FALSE
    )
  }
  samples <- unique(assignments[, c("dataset", "rawcolumnnames"), with = FALSE])
  assigned <- unique(
    assignments[!is.na(assignments$condition)],
    by = c("dataset", "rawcolumnnames")
  )
  resolved <- merge(samples, assigned, by = c("dataset", "rawcolumnnames"),
                    all.x = TRUE, sort = FALSE)
  unassigned <- is.na(resolved$condition)
  if (any(unassigned)) {
    message(
      sum(unassigned),
      " sample(s) are excluded from every comparison and are labelled 'Unassigned'."
    )
    resolved[unassigned, condition := "Unassigned"]
  }
  data.table::setcolorder(
    resolved,
    c("dataset", "rawcolumnnames", "condition", "comparison", "comparison_label")
  )
  resolved[]
}
