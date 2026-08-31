library(targets)
library(MultiOmicsDataCompile)

data_dir <- Sys.getenv("MULTIOMICS_DATA_DIR", unset = "")
if (!nzchar(data_dir) || !dir.exists(data_dir)) {
  stop("Set MULTIOMICS_DATA_DIR to the processing directory before running tar_make().")
}

tar_option_set(
  packages = c("MultiOmicsDataCompile", "data.table", "SummarizedExperiment"),
  format = "rds"
)

list(
  tar_target(
    overview_file,
    file.path(data_dir, "OverviewFiles", "GEODataOverview3.csv"),
    format = "file"
  ),
  tar_target(
    overview_all,
    data.table::as.data.table(
      ValidateDatasetManifest(data.table::fread(overview_file))
    )
  ),
  tar_target(overview_geo, overview_all[GEO2R == TRUE]),
  tar_target(
    platform_info_file,
    file.path(data_dir, "OverviewFiles", "GEOPlatformInfo.csv"),
    format = "file"
  ),
  tar_target(
    platform_annotation,
    PlatformAnnotationLoad(data.table::fread(platform_info_file))
  ),
  tar_target(
    geo_outputs,
    {
      GEOCompile(
        DS = overview_geo$ID,
        gpl = overview_geo$GPLNumber,
        gsm = overview_geo$ComparisonVector,
        namestr = overview_geo$FCColumnNames,
        nameraw = overview_geo$RawColumnNames,
        PlatAnnotInfo = platform_annotation,
        destdir = file.path(data_dir, "GEOcache"),
        writeDB = TRUE,
        writeRaw = overview_geo$DownloadRawData,
        GenerateMetaData = overview_geo$DownloadMetaData,
        MetaDataPath = file.path(data_dir, "ArrayMetaData"),
        writeMetaData = TRUE,
        DBPath = file.path(data_dir, "AppData"),
        Technology = overview_geo$Technology,
        Species = overview_geo$Species,
        subsetRaw = FALSE
      )
      list.files(file.path(data_dir, "AppData"), full.names = TRUE)
    },
    format = "file"
  ),
  tar_target(
    de_database,
    {
      geo_outputs
      output <- file.path(data_dir, "ProcessFiles", "SumarizedExp_DB.rds")
      DESEGenerate(file.path(data_dir, "AppData"), output)
      output
    },
    format = "file"
  ),
  tar_target(
    compiled_expression_file,
    {
      geo_outputs
      output <- file.path(data_dir, "ProcessFiles", "RawData.txt")
      RawDataCompile(
        file.path(data_dir, "AppData"), output, overview = overview_all
      )
      output
    },
    format = "file"
  ),
  tar_target(
    sample_conditions,
    {
      geo_outputs
      SampleConditions(
        DBPath = file.path(data_dir, "AppData"), overview = overview_all
      )
    }
  ),
  tar_target(
    sample_metadata,
    MetDataCompile(
      RNAseqFilePath = file.path(data_dir, "RunInfo"),
      ArrayFilePath = file.path(data_dir, "ArrayMetaData"),
      overview = overview_all,
      conditions = sample_conditions
    )
  ),
  tar_target(
    expression_collection,
    GenerateRawSE(
      df = sample_metadata,
      ArrayDT = data.table::fread(compiled_expression_file),
      overview = overview_all
    )
  ),
  tar_target(
    expression_collection_file,
    {
      output <- file.path(data_dir, "ProcessFiles", "ExpressionCollection.rds")
      saveRDS(expression_collection, output)
      output
    },
    format = "file"
  ),
  tar_target(
    exploration_file,
    {
      output <- file.path(data_dir, "ProcessFiles", "expression_norm.v2.RDS")
      saveRDS(expression_collection$ExplorationSE, output)
      output
    },
    format = "file"
  )
)
