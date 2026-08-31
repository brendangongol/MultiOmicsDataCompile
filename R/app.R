#' Launch the MultiOmicsDataCompile explorer
#'
#' Starts the packaged Shiny application against an existing processing
#' directory. The directory is passed through an R option, so packaged source
#' files do not need to be rewritten before every launch.
#'
#' @param data_dir A directory created by [makeDirectory()] and populated by the
#'   processing workflow.
#' @param ... Additional arguments passed to [shiny::runApp()].
#'
#' @return The value returned by [shiny::runApp()], invisibly.
#' @examples
#' \dontrun{
#' launchMultiOmicsExplorer("path/to/MultiOmicsData")
#' }
#' @family exploration functions
#' @export
launchMultiOmicsExplorer <- function(data_dir, ...) {
  app_packages <- c(
    "shiny", "plotly", "shinyWidgets", "DT", "RColorBrewer", "DEP"
  )
  missing_packages <- app_packages[
    !vapply(app_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop(
      "Install the following optional packages to launch the explorer: ",
      paste(missing_packages, collapse = ", "),
      call. = FALSE
    )
  }
  data_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  required <- c(
    "ProcessFiles", "OverviewFiles", "Proteomic_3", "Metabolomics",
    "Methylation"
  )
  missing <- required[!dir.exists(file.path(data_dir, required))]
  if (length(missing)) {
    stop(
      "`data_dir` is missing required directories: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  required_files <- c(
    file.path("OverviewFiles", "GEODataOverview3.csv"),
    file.path("ProcessFiles", "SumarizedExp_DB.rds"),
    file.path("ProcessFiles", "expression_norm.v2.RDS"),
    file.path("ProcessFiles", "SumarizedProtExp_DB.rds"),
    file.path("OverviewFiles", "ProteomicProteins.RDS"),
    file.path("Metabolomics", "mtbls298.de.RDS")
  )
  missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]
  if (length(missing_files)) {
    stop(
      "The explorer cannot start because processed files are missing: ",
      paste(missing_files, collapse = ", "),
      ". Run the relevant workflow stages first.",
      call. = FALSE
    )
  }
  app_dir <- system.file(
    "extdata", "Scripts", package = "MultiOmicsDataCompile", mustWork = TRUE
  )
  old_option <- getOption("MultiOmicsDataCompile.data_dir")
  options(MultiOmicsDataCompile.data_dir = data_dir)
  on.exit(options(MultiOmicsDataCompile.data_dir = old_option), add = TRUE)
  invisible(shiny::runApp(appDir = app_dir, ...))
}
