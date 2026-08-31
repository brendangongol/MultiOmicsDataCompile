# MultiOmicsDataCompile 0.2.1

## Pipeline fixes

- Added `SampleConditions()`, which derives each sample's condition from the
  manifest `ComparisonVector` and the stored raw-expression columns.
  `MetDataCompile()` gains a `conditions` argument that attaches the result, so
  the transcriptomics workflow now reaches `GenerateRawSE()` instead of stopping
  for want of a `condition` column. Conditions remain explicit; they are still
  never inferred from sample order.
- `RawDataCompile()` now also reads the `*_CountRaw.rds` files written by
  `ExternalDataHarmonize()`, which were previously skipped without warning.
  `*_RPKMRaw.rds` files are still excluded, because they are not on the count
  scale that `GenerateRawSE()` requires, but the exclusion is now reported.
  Because those samples now reach `GenerateRawSE()`, they need sample metadata
  and conditions like any other dataset; `GenerateRawSE()` names the datasets
  whose metadata are missing instead of reporting sample identifiers alone.
- Fixed `GEO2RDirectionCheck()` for RNA-seq datasets. The direction-only fold
  change was computed from count-filtered rows but paired with unfiltered gene
  symbols, which raised a length error or silently misaligned genes. Row
  identities are now carried through the filter. The check no longer runs a full
  Wald fit, since only size factors are needed.
- Rewrote `DESEProtGenerate()`. Updating an existing database previously always
  failed: a local `rowData` variable shadowed the `SummarizedExperiment`
  generic, an annotation table from the other code path was referenced, and the
  two paths produced assays with different column sets. The scratch and update
  paths now share one implementation, and the "compile from scratch" test no
  longer relies on a warning handler or leaks a `Scratch` variable into the
  global environment.

## Explorer fixes

- Replaced the undefined `meltAssay()` calls with a new exported
  `MeltExpressionAssay()`, so both violin-plot tabs work without an
  undeclared single-cell dependency.
- The multi-dataset violin plot no longer requires a `facet_name` column that
  nothing created; it facets by dataset. Both violin tabs now attach only the
  sample annotations that are actually present.
- Gave the single- and multi-dataset submit buttons distinct input identifiers.
- The application now calls the packaged `add_rejections2()`, `plot_volcano2()`,
  `ProteomicSELoad()`, and `SigDEAnnotate()` instead of carrying copies that had
  drifted; the copied `plot_volcano2()` had disabled its own argument checks.

## Packaging

- Declared `tools` in `Imports`, and added `CONTRIBUTING.md` to
  `.Rbuildignore`.
- Documented that the Bioconductor repositories must be registered before
  installing, and documented the sample-condition stage.

# MultiOmicsDataCompile 0.2.0

- Renamed the package from CMDMultiOmics to MultiOmicsDataCompile.
- Added strict dataset-manifest, comparison-vector, and sample-metadata validation.
- Updated GEO RNA-seq downloads to use the current GEOquery RNA-seq API and retain source provenance.
- Added optional shrunken RNA-seq log2 fold changes while retaining maximum-likelihood estimates.
- Separated RNA-seq count, normalized-count, transformed, array-expression, and cross-dataset exploration representations.
- Added RNA-seq library-size, detection, correlation, distance, and PCA quality-control summaries.
- Rebuilt differential-expression compilation around species-qualified stable identifiers.
- Added a supported Shiny launcher, a `targets` workflow, tests, and continuous integration.
- Deprecated `DetermineMising()` and `LowSampleCountRmove()` in favor of their correctly spelled equivalents.

