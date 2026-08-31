.require_optional_package <- function(package, feature) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "`", feature, "` requires the optional package `", package,
      "`. Install it before using this feature.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}
