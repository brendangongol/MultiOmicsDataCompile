# MultiOmicsDataCompile

`MultiOmicsDataCompile` provides batch-processing workflows for public
transcriptomic and proteomic data, tools for integrating differential and raw
expression data, and a Shiny application for brief multi-omics exploration.
The bundled example database focuses on cardiometabolic disease and includes
transcriptomic, proteomic, methylation, and metabolomic datasets.

## Installation

Most required dependencies (`DESeq2`, `GEOquery`, `limma`, `Biobase`,
`BiocFileCache`, `S4Vectors`, `SummarizedExperiment`, `rpx`, and the
`org.*.eg.db` annotation packages) come from Bioconductor rather than CRAN, so
the Bioconductor repositories must be on the search path before installing.
Without this step `install_github()` cannot resolve them:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

options(repos = BiocManager::repositories())

remotes::install_github(
  "brendangongol/MultiOmicsDataCompile",
  dependencies = TRUE
)
```

The transcriptomics workflow uses the required dependencies installed above.
Proteomics normalization, imputation, differential analysis, and the explorer
also use the suggested Bioconductor packages `DEP`, `MSnbase`, and
`NormalyzerDE`; install those suggestions when you need the proteomics stages:

```r
BiocManager::install(c("DEP", "MSnbase", "NormalyzerDE"))
```

Microarray platform annotation additionally needs the platform-specific
`.db` packages listed under `Suggests`; `PlatformAnnotationLoad()` reports the
name of any package it needs but cannot find.

Load the package with its new package name:

```r
library(MultiOmicsDataCompile)
```

## Set up the processing workspace

Choose a writable directory, copy the packaged resources, and configure the
copied scripts:

```r
library(MultiOmicsDataCompile)

home_dir <- "path/to/MultiOmicsDataCompileAnalysis"
makeDirectory(home_dir)
AppSetup(home_dir)
```

The configured processing driver is
`Scripts/WorkflowScript.R` beneath `home_dir`. Review its dataset selections and
then run it from that directory. Several stages download large public datasets
and should be run deliberately rather than during package installation.

For resumable execution, the copied `Scripts/_targets.R` file provides the same
transcriptomics stages as a dependency-tracked pipeline:

```r
Sys.setenv(MULTIOMICS_DATA_DIR = normalizePath(home_dir))
targets::tar_make(
  script = file.path(home_dir, "Scripts", "_targets.R"),
  store = file.path(home_dir, "Scripts", "_targets")
)
```

## Sample conditions

`GenerateRawSE()` requires a `condition` for every RNA-seq sample and will stop
rather than guess one. Neither SRA run tables nor GEO sample metadata record a
usable condition, so it is derived from the manifest instead:
`SampleConditions()` reads the stored raw-expression files and applies each
dataset's `ComparisonVector`, which positions one character per sample of that
matrix. Both processing drivers run this step before `MetDataCompile()`:

```r
conditions <- SampleConditions(
  DBPath = file.path(home_dir, "AppData"),
  overview = overview_all
)

metaDF <- MetDataCompile(
  RNAseqFilePath = file.path(home_dir, "RunInfo"),
  ArrayFilePath = file.path(home_dir, "ArrayMetaData"),
  overview = overview_all,
  conditions = conditions
)
```

Where a dataset contributes several comparisons, the first comparison that
includes a sample names its condition. Samples that every comparison excludes
(`X`) are labelled `"Unassigned"`, which records their exclusion without
implying group membership. Review the returned table, or the copied
`OverviewFiles/SampleConditions.xls`, before building the expression objects; a
hand-edited table with `rawcolumnnames` and `condition` columns can be passed to
`MetDataCompile()` instead.

The extended workflow is available in the bundled vignette source:

```r
file.show(system.file(
  "extdata",
  "vignette",
  "MultiOmicsDataCompileVignette.Rmd",
  package = "MultiOmicsDataCompile"
))
```

## Scope

The workflow covers:

- GEO microarray and RNA-seq download, harmonization, and two-group analysis;
- ProteomeXchange download and dataset-specific MaxQuant formatting;
- construction of integrated `SummarizedExperiment` databases;
- missing-value filtering, normalization, imputation, and differential protein
  abundance analysis; and
- lightweight exploration with tables, violin plots, volcano plots, heatmaps,
  and cross-dataset summaries.

RNA-seq counts and microarray expression are retained in separate per-dataset
experiments. The combined app-facing expression object contains within-dataset
z-scores for qualitative exploration only; it is not suitable for differential
testing or absolute abundance comparisons across platforms. RNA-seq sample
conditions must be supplied explicitly and are never inferred from sample
order.

Launch the explorer without editing packaged scripts:

```r
launchMultiOmicsExplorer(home_dir)
```

![Multi-omics workflow overview](images/Figure1.png)

![Exploration features](images/Figure2.png)

## Development status

The package is under active development. Proteomics source files vary widely in
format, so `FormatMaxQuant()` currently contains dataset-specific parsers and
new repositories may require an additional parser. Review generated quality
control outputs before adding processed data to a production database.
