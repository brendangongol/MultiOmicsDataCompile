# Contributing

Thank you for improving MultiOmicsDataCompile. Please open an issue before a large behavioral change so that its data-model and compatibility effects can be discussed.

For a code contribution:

1. Create a focused branch and add or update tests under `tests/testthat/`.
2. Document exported functions with roxygen2, including parameters, return values, and a runnable example when practical.
3. Run `roxygen2::roxygenise()`, `testthat::test_local()`, and `rcmdcheck::rcmdcheck()`.
4. Do not commit downloaded GEO data, credentials, caches, or machine-specific absolute paths.
5. Describe user-visible changes in `NEWS.md`.

Documentation is generated with roxygen2 8.1.0, which records its version in the
`Config/roxygen2/version` field of `DESCRIPTION`. Regenerating with roxygen2
7.x replaces that field with `RoxygenNote` and rewrites `NAMESPACE` into the
older one-import-per-line format, so check your roxygen2 version before
committing generated files.

New RNA-seq changes should preserve raw integer counts, record organism and genome provenance, validate sample-to-condition mappings, and avoid combining count-scale and expression-scale assays in one matrix.

Sample conditions are derived from the manifest `ComparisonVector` by
`SampleConditions()`. Do not add code paths that infer group membership from
sample order, column position, or file order.
