# Internal helpers for transcriptomics workflows.

.normalize_species <- function(x) {
  value <- tolower(trimws(as.character(x)))
  value[value %in% c("human", "homo sapiens", "9606")] <- "Human"
  value[value %in% c("mouse", "mus musculus", "10090")] <- "Mouse"
  value[value %in% c("rat", "rattus norvegicus", "10116")] <- "Rat"
  value[!value %in% c("Human", "Mouse", "Rat")] <- NA_character_
  value
}

.annotation_object <- function(package, object) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(
      "Platform annotation requires the optional package `", package,
      "`. Install it with BiocManager::install(\"", package, "\").",
      call. = FALSE
    )
  }
  get(object, envir = asNamespace(package), inherits = FALSE)
}

.validate_parallel_lengths <- function(reference, ..., names = NULL) {
  values <- list(...)
  expected <- length(reference)
  if (is.null(names)) {
    names <- paste0("argument ", seq_along(values))
  }
  bad <- vapply(
    values,
    function(x) !(length(x) %in% c(1L, expected)),
    logical(1)
  )
  if (any(bad)) {
    stop(
      "The following arguments must have length 1 or length(DS): ",
      paste(names[bad], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.recycle_to_length <- function(x, n) {
  if (length(x) == n) x else rep(x, n)
}

.parse_comparison_code <- function(code, sample_names, dataset) {
  labels <- .comparison_labels(code, dataset)
  if (length(labels) != length(sample_names)) {
    stop(
      "ComparisonVector for ", dataset, " has ", length(labels),
      " entries but the count matrix has ", length(sample_names), " samples.",
      call. = FALSE
    )
  }
  selected <- labels != "X"
  selected_labels <- labels[selected]
  if (!identical(sort(unique(selected_labels)), c("0", "1"))) {
    stop(
      "ComparisonVector for ", dataset,
      " must select at least one control (0) and one treatment (1) sample.",
      call. = FALSE
    )
  }
  list(labels = labels, selected = selected)
}

.find_column <- function(x, candidates, required = TRUE) {
  normalized <- tolower(gsub("[^a-z0-9]", "", names(x)))
  candidate_names <- tolower(gsub("[^a-z0-9]", "", candidates))
  index <- match(candidate_names, normalized, nomatch = 0L)
  index <- index[index > 0L]
  if (!length(index)) {
    if (required) {
      stop(
        "Required column not found. Expected one of: ",
        paste(candidates, collapse = ", "),
        call. = FALSE
      )
    }
    return(NULL)
  }
  names(x)[index[1]]
}

.standardize_geo_rnaseq_annotation <- function(annotation, gene_ids) {
  annotation <- as.data.frame(annotation, stringsAsFactors = FALSE)
  if (!nrow(annotation)) {
    annotation <- data.frame(row.names = as.character(gene_ids))
  }
  if (is.null(rownames(annotation)) || any(!nzchar(rownames(annotation)))) {
    id_column <- .find_column(annotation, c("GeneID", "gene_id", "ENTREZID"))
    rownames(annotation) <- as.character(annotation[[id_column]])
  }
  annotation <- annotation[match(as.character(gene_ids), rownames(annotation)), , drop = FALSE]
  symbol_column <- .find_column(
    annotation,
    c("Symbol", "gene_symbol", "GeneSymbol", "SYMBOL"),
    required = FALSE
  )
  description_column <- .find_column(
    annotation,
    c("Description", "gene_description", "GENENAME", "description"),
    required = FALSE
  )
  data.frame(
    ENTREZID = as.character(gene_ids),
    SYMBOL = if (is.null(symbol_column)) NA_character_ else as.character(annotation[[symbol_column]]),
    GENENAME = if (is.null(description_column)) NA_character_ else as.character(annotation[[description_column]]),
    stringsAsFactors = FALSE,
    row.names = as.character(gene_ids)
  )
}

.expression_provenance <- function(dataset, species = NA_character_, genome = NA_character_,
                                   platform = NA_character_, source = NA_character_) {
  list(
    dataset = dataset,
    species = species,
    genome = genome,
    platform = platform,
    source = source,
    retrieved_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    package_version = tryCatch(
      as.character(utils::packageVersion("MultiOmicsDataCompile")),
      error = function(e) NA_character_
    ),
    r_version = R.version.string
  )
}

.validate_expression_metadata <- function(df, sample_names) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  if (is.null(rownames(df)) || any(!nzchar(rownames(df)))) {
    stop("`df` must have non-empty sample identifiers as row names.", call. = FALSE)
  }
  if (anyDuplicated(rownames(df))) {
    stop("Sample identifiers in `df` must be unique.", call. = FALSE)
  }
  required <- c("dataset", "assay.type", "rawcolumnnames")
  missing <- setdiff(required, names(df))
  if (length(missing)) {
    stop(
      "`df` is missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  missing_metadata <- setdiff(sample_names, rownames(df))
  if (length(missing_metadata)) {
    datasets <- unique(sub("_.*$", "", missing_metadata))
    stop(
      "Metadata are missing for expression samples: ",
      paste(utils::head(missing_metadata, 10L), collapse = ", "),
      if (length(missing_metadata) > 10L) " ..." else "",
      ". Affected dataset(s): ", paste(utils::head(datasets, 10L), collapse = ", "),
      ". Every compiled expression column needs a row in `df`; add a sample ",
      "metadata table for these datasets, or remove their raw-data files from ",
      "the compilation directory.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.row_zscores <- function(x) {
  x <- as.matrix(x)
  means <- rowMeans(x, na.rm = TRUE)
  sds <- apply(x, 1L, stats::sd, na.rm = TRUE)
  result <- sweep(x, 1L, means, "-")
  usable <- is.finite(sds) & sds > 0
  result[usable, ] <- sweep(result[usable, , drop = FALSE], 1L, sds[usable], "/")
  result[!usable, ] <- 0
  result
}

.bind_feature_matrices <- function(matrices) {
  if (!length(matrices)) {
    return(NULL)
  }
  feature_ids <- Reduce(union, lapply(matrices, rownames))
  aligned <- lapply(matrices, function(x) {
    output <- matrix(
      NA_real_, nrow = length(feature_ids), ncol = ncol(x),
      dimnames = list(feature_ids, colnames(x))
    )
    output[rownames(x), ] <- x
    output
  })
  do.call(cbind, aligned)
}

.infer_species_from_entrez <- function(ids) {
  ids <- unique(as.character(ids[!is.na(ids)]))
  if (!length(ids)) return(NA_character_)
  human <- suppressMessages(AnnotationDbi::keys(org.Hs.eg.db, keytype = "ENTREZID"))
  mouse <- suppressMessages(AnnotationDbi::keys(org.Mm.eg.db, keytype = "ENTREZID"))
  overlap <- c(Human = sum(ids %in% human), Mouse = sum(ids %in% mouse))
  if (max(overlap) == 0L) NA_character_ else names(which.max(overlap))
}

.write_table_once <- function(x, path) {
  directory <- dirname(path)
  if (!dir.exists(directory)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  }
  temporary <- tempfile("multiomics-", tmpdir = directory, fileext = ".tmp")
  on.exit(unlink(temporary), add = TRUE)
  data.table::fwrite(x, temporary, row.names = FALSE, quote = FALSE, sep = "\t")
  if (!file.copy(temporary, path, overwrite = TRUE)) {
    stop("Unable to write compiled expression table to `", path, "`.", call. = FALSE)
  }
  invisible(path)
}

.raw_data_compile_impl <- function(Fpath, outPath, StartAt = 1, overview = NULL) {
  existing <- if (file.exists(outPath)) {
    data.table::fread(outPath)
  } else {
    NULL
  }
  # GEOCompile() writes "<dataset>_Raw.rds" and ExternalDataHarmonize() writes
  # "<dataset>_CountRaw.rds"; both are count- or expression-scale inputs.
  files <- list.files(Fpath, pattern = "_(Count)?Raw[.]rds$", full.names = TRUE)
  rpkm_files <- list.files(Fpath, pattern = "_RPKMRaw[.]rds$", full.names = TRUE)
  if (length(rpkm_files)) {
    message(
      "Skipping ", length(rpkm_files),
      " RPKM/FPKM file(s), which are not on the count scale expected by ",
      "GenerateRawSE(): ",
      paste(basename(rpkm_files), collapse = ", ")
    )
  }
  if (!length(files) || StartAt > length(files)) {
    message("No raw-data files remain to be processed.")
    return(invisible(existing))
  }
  files <- files[seq.int(StartAt, length(files))]
  for (index in seq_along(files)) {
    path <- files[index]
    dataset <- sub("_.*$", "", basename(path))
    source_object <- readRDS(path)
    provenance <- attr(source_object, "provenance", exact = TRUE)
    one <- data.table::as.data.table(source_object)
    required <- c("ENTREZID", "SYMBOL")
    missing <- setdiff(required, names(one))
    if (length(missing)) {
      stop(
        basename(path), " is missing required columns: ",
        paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    one <- one[!is.na(ENTREZID) & !is.na(SYMBOL) & nzchar(SYMBOL)]
    sample_columns <- setdiff(names(one), c("ID", "ENTREZID", "SYMBOL", "GENENAME"))
    if (!length(sample_columns)) {
      warning("No sample columns found in ", basename(path), ".", call. = FALSE)
      next
    }
    for (column in sample_columns) {
      one[[column]] <- suppressWarnings(as.numeric(one[[column]]))
    }
    dataset_row <- if (!is.null(overview) && "ID" %in% names(overview)) {
      overview_df <- as.data.frame(overview, stringsAsFactors = FALSE)
      overview_df[overview_df$ID == dataset, , drop = FALSE]
    } else {
      data.frame()
    }
    species <- if (nrow(dataset_row) && "Species" %in% names(dataset_row)) {
      .normalize_species(dataset_row$Species[1])
    } else if (!is.null(provenance) && !is.null(provenance$species)) {
      .normalize_species(provenance$species)[1]
    } else {
      .infer_species_from_entrez(one$ENTREZID)
    }
    technology <- if (nrow(dataset_row) && "Technology" %in% names(dataset_row)) {
      tolower(dataset_row$Technology[1])
    } else if (!is.null(provenance) && !is.null(provenance$source) &&
               grepl("RNA-seq|counts", provenance$source, ignore.case = TRUE)) {
      "rnaseq"
    } else {
      if (any(grepl("CountRaw|RNA", basename(path), ignore.case = TRUE))) "rnaseq" else "array"
    }
    if (is.na(species)) {
      stop("Unable to determine species for ", dataset, ".", call. = FALSE)
    }
    aggregate_function <- if (technology %in% c("rnaseq", "rna-seq", "rna seq")) {
      function(x) sum(x, na.rm = TRUE)
    } else {
      function(x) mean(x, na.rm = TRUE)
    }
    identifier_columns <- intersect(c("ENTREZID", "SYMBOL", "GENENAME"), names(one))
    identifiers <- unique(one[, ..identifier_columns], by = "ENTREZID")
    values <- one[, lapply(.SD, aggregate_function), by = ENTREZID,
                  .SDcols = sample_columns]
    one <- merge(identifiers, values, by = "ENTREZID", all.y = TRUE)
    one[, `:=`(
      SPECIES = species,
      FEATURE_ID = paste(species, ENTREZID, sep = ":")
    )]
    data.table::setcolorder(one, c("FEATURE_ID", "SPECIES", "ENTREZID", "SYMBOL", sample_columns))

    if (is.null(existing)) {
      existing <- one
    } else {
      if (!"FEATURE_ID" %in% names(existing)) {
        if (!any(c("ENTREZID_Human", "ENTREZID_Mouse") %in% names(existing))) {
          stop(
            "Existing compiled data lack FEATURE_ID and legacy Entrez columns.",
            call. = FALSE
          )
        }
        human_id <- if ("ENTREZID_Human" %in% names(existing)) existing$ENTREZID_Human else NA
        mouse_id <- if ("ENTREZID_Mouse" %in% names(existing)) existing$ENTREZID_Mouse else NA
        existing[, `:=`(
          SPECIES = ifelse(!is.na(human_id), "Human", "Mouse"),
          ENTREZID = as.character(ifelse(!is.na(human_id), human_id, mouse_id)),
          FEATURE_ID = paste(SPECIES, ENTREZID, sep = ":")
        )]
      }
      duplicate_samples <- intersect(sample_columns, names(existing))
      for (sample in duplicate_samples) {
        comparison <- merge(
          one[, c("FEATURE_ID", sample), with = FALSE],
          existing[, c("FEATURE_ID", sample), with = FALSE],
          by = "FEATURE_ID", suffixes = c(".new", ".existing")
        )
        new_features <- one$FEATURE_ID[!is.na(one[[sample]])]
        existing_features <- existing$FEATURE_ID[!is.na(existing[[sample]])]
        same_features <- setequal(new_features, existing_features)
        same <- same_features && isTRUE(all.equal(
          comparison[[paste0(sample, ".new")]],
          comparison[[paste0(sample, ".existing")]],
          tolerance = 0, check.attributes = FALSE
        ))
        if (same) {
          one[, (sample) := NULL]
        } else {
          replacement <- make.unique(c(names(existing), sample))[length(names(existing)) + 1L]
          data.table::setnames(one, sample, replacement)
        }
      }
      annotation <- intersect(
        c("FEATURE_ID", "SPECIES", "ENTREZID", "SYMBOL", "GENENAME"),
        unique(c(names(existing), names(one)))
      )
      if (ncol(one) > length(annotation)) {
        annotation_existing <- existing[, intersect(annotation, names(existing)), with = FALSE]
        annotation_new <- one[, intersect(annotation, names(one)), with = FALSE]
        annotation_rows <- unique(
          data.table::rbindlist(
            list(annotation_existing, annotation_new), use.names = TRUE, fill = TRUE
          ),
          by = "FEATURE_ID"
        )
        existing_values <- existing[, setdiff(names(existing), setdiff(annotation, "FEATURE_ID")), with = FALSE]
        one_values <- one[, setdiff(names(one), setdiff(annotation, "FEATURE_ID")), with = FALSE]
        existing <- merge(existing_values, one_values, by = "FEATURE_ID", all = TRUE)
        existing <- merge(
          existing[, setdiff(names(existing), setdiff(annotation, "FEATURE_ID")), with = FALSE],
          annotation_rows, by = "FEATURE_ID", all.x = TRUE
        )
        data.table::setcolorder(
          existing,
          c(annotation, setdiff(names(existing), annotation))
        )
      }
    }
    message("Incorporated ", index, " of ", length(files), ": ", basename(path))
  }
  .write_table_once(existing, outPath)
  invisible(existing)
}

.dataset_species <- function(overview, dataset) {
  if (!all(c("ID", "Species") %in% names(overview))) return(NA_character_)
  value <- overview$Species[match(dataset, overview$ID)]
  .normalize_species(value[1])
}

.generate_expression_collection <- function(df, ArrayDT, overview) {
  ArrayDT <- data.table::as.data.table(ArrayDT)
  annotation_columns <- intersect(
    c("FEATURE_ID", "SPECIES", "ENTREZID", "ENTREZID_Human", "ENTREZID_Mouse",
      "SYMBOL", "GENENAME", "Length_Human", "Length_Mouse"),
    names(ArrayDT)
  )
  sample_names <- setdiff(names(ArrayDT), annotation_columns)
  .validate_expression_metadata(df, sample_names)
  if (!identical(as.character(df$rawcolumnnames), rownames(df))) {
    stop(
      "`df$rawcolumnnames` must match its sample row names in the same order.",
      call. = FALSE
    )
  }
  df <- df[match(sample_names, rownames(df)), , drop = FALSE]
  if (!identical(rownames(df), sample_names)) {
    stop("Expression columns and metadata rows could not be aligned.", call. = FALSE)
  }
  if (!"FEATURE_ID" %in% names(ArrayDT)) {
    species <- ifelse(!is.na(ArrayDT$ENTREZID_Human), "Human", "Mouse")
    entrez <- ifelse(!is.na(ArrayDT$ENTREZID_Human),
                     ArrayDT$ENTREZID_Human, ArrayDT$ENTREZID_Mouse)
    ArrayDT[, `:=`(
      FEATURE_ID = make.unique(paste(species, entrez, sep = ":")),
      SPECIES = species,
      ENTREZID = as.character(entrez)
    )]
    annotation_columns <- unique(c("FEATURE_ID", "SPECIES", "ENTREZID", annotation_columns))
  }
  if (anyDuplicated(ArrayDT$FEATURE_ID)) {
    stop("Compiled expression data contain duplicated FEATURE_ID values.", call. = FALSE)
  }
  row_annotation <- as.data.frame(
    ArrayDT[, intersect(
      c("FEATURE_ID", "SPECIES", "ENTREZID", "ENTREZID_Human", "ENTREZID_Mouse",
        "SYMBOL", "GENENAME"), names(ArrayDT)
    ), with = FALSE]
  )
  rownames(row_annotation) <- row_annotation$FEATURE_ID
  matrix_all <- as.matrix(ArrayDT[, ..sample_names])
  storage.mode(matrix_all) <- "numeric"
  rownames(matrix_all) <- ArrayDT$FEATURE_ID

  rna_experiments <- list()
  array_experiments <- list()
  exploration_matrices <- list()
  datasets <- unique(as.character(df$dataset))
  for (dataset in datasets) {
    columns <- rownames(df)[df$dataset == dataset]
    columns <- intersect(columns, colnames(matrix_all))
    if (!length(columns)) next
    metadata <- df[columns, , drop = FALSE]
    assay_type <- unique(as.character(metadata$assay.type))
    assay_type <- assay_type[!is.na(assay_type)]
    if (length(assay_type) != 1L) {
      stop("Dataset ", dataset, " must have exactly one assay type.", call. = FALSE)
    }
    values <- matrix_all[, columns, drop = FALSE]
    dataset_species <- .dataset_species(overview, dataset)
    species_rows <- if (!is.na(dataset_species) && "SPECIES" %in% names(row_annotation)) {
      row_annotation$SPECIES == dataset_species
    } else {
      rep(TRUE, nrow(row_annotation))
    }
    species_rows[is.na(species_rows)] <- FALSE
    values <- values[species_rows, , drop = FALSE]
    assay_is_rnaseq <- tolower(assay_type) %in% c("rna-seq", "rnaseq", "rna seq")
    keep <- if (assay_is_rnaseq) {
      rowSums(!is.na(values)) == ncol(values)
    } else {
      rowSums(!is.na(values)) > 0L
    }
    values <- values[keep, , drop = FALSE]
    feature_data <- row_annotation[rownames(values), , drop = FALSE]
    if (assay_is_rnaseq) {
      if (!"condition" %in% names(metadata) ||
          anyNA(metadata$condition) || any(!nzchar(trimws(as.character(metadata$condition))))) {
        stop(
          "RNA-seq condition metadata are missing for dataset ", dataset,
          ". Conditions must be supplied explicitly; they are never inferred.",
          call. = FALSE
        )
      }
      if (any(values < 0, na.rm = TRUE) ||
          any(abs(values - round(values)) > sqrt(.Machine$double.eps), na.rm = TRUE)) {
        stop("RNA-seq data for ", dataset, " must be non-negative integer counts.", call. = FALSE)
      }
      counts_matrix <- round(values)
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts_matrix,
        colData = metadata,
        design = ~ 1
      )
      dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
      normalized <- DESeq2::counts(dds, normalized = TRUE)
      transformation <- if (nrow(counts_matrix) < 20L) "log2 normalized counts" else "DESeq2 VST"
      transformed <- if (nrow(counts_matrix) < 20L) {
        log2(normalized + 1)
      } else {
        tryCatch(
          SummarizedExperiment::assay(
            DESeq2::varianceStabilizingTransformation(dds, blind = TRUE)
          ),
          error = function(e) {
            transformation <<- "log2 normalized counts (VST fallback)"
            log2(normalized + 1)
          }
        )
      }
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(
          counts = counts_matrix,
          normalized_counts = normalized,
          transformed = transformed
        ),
        rowData = feature_data,
        colData = metadata,
        metadata = list(
          dataset = dataset,
          species = dataset_species,
          transformation = transformation,
          units = c(counts = "integer counts", normalized_counts = "DESeq2 size-factor normalized counts", transformed = transformation)
        )
      )
      rna_experiments[[dataset]] <- se
      exploration_matrices[[dataset]] <- .row_zscores(transformed)
    } else if (tolower(assay_type) == "array") {
      se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(expression = values),
        rowData = feature_data,
        colData = metadata,
        metadata = list(
          dataset = dataset,
          species = dataset_species,
          units = "source-normalized microarray expression"
        )
      )
      array_experiments[[dataset]] <- se
      exploration_matrices[[dataset]] <- .row_zscores(values)
    } else {
      stop("Unsupported assay type `", assay_type, "` for ", dataset, ".", call. = FALSE)
    }
  }
  exploration <- .bind_feature_matrices(exploration_matrices)
  if (is.null(exploration)) {
    stop("No RNA-seq or array datasets were available.", call. = FALSE)
  }
  exploration_coldata <- df[match(colnames(exploration), rownames(df)), , drop = FALSE]
  exploration_rowdata <- row_annotation[rownames(exploration), , drop = FALSE]
  exploration_se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(Expression = exploration),
    rowData = exploration_rowdata,
    colData = exploration_coldata,
    metadata = list(
      purpose = "Cross-dataset qualitative exploration only",
      units = "within-dataset feature z-score",
      warning = "Do not use this assay for differential-expression inference or absolute cross-platform abundance comparisons."
    )
  )
  collection <- list(RNAseq = rna_experiments, Array = array_experiments)
  list(
    RawSE = collection,
    RNAseqSE = rna_experiments,
    ArraySE = array_experiments,
    ExplorationSE = exploration_se,
    ExpressionSE = exploration_se,
    RPKMSE = exploration_se
  )
}

.standardize_de_result <- function(x, source_name) {
  provenance <- attr(x, "provenance", exact = TRUE)
  x <- data.table::as.data.table(x)
  names(x) <- sub("^.+_(logFC|Pvalue|AdjPValue)$", "\\1", names(x))
  required_statistics <- c("logFC", "Pvalue", "AdjPValue")
  missing_statistics <- setdiff(required_statistics, names(x))
  if (length(missing_statistics)) {
    stop(
      source_name, " is missing differential-expression columns: ",
      paste(missing_statistics, collapse = ", "),
      call. = FALSE
    )
  }
  stable_column <- if ("ENTREZID" %in% names(x)) {
    "ENTREZID"
  } else if ("ID" %in% names(x)) {
    "ID"
  } else {
    stop(source_name, " must contain `ENTREZID` or `ID`.", call. = FALSE)
  }
  stable_id <- as.character(x[[stable_column]])
  species <- if (!is.null(provenance) && !is.null(provenance$species)) {
    .normalize_species(provenance$species)
  } else if (any(grepl("^ENSMUSG", stable_id))) {
    "Mouse"
  } else if (any(grepl("^ENSRNOG", stable_id))) {
    "Rat"
  } else if (any(grepl("^ENSG", stable_id))) {
    "Human"
  } else {
    .infer_species_from_entrez(stable_id)
  }
  if (length(species) != 1L || is.na(species)) {
    stop("Unable to determine species for ", source_name, ".", call. = FALSE)
  }
  x[, FEATURE_ID := paste(species, stable_id, sep = ":")]
  x[, SPECIES := species]
  if (!"SYMBOL" %in% names(x)) x[, SYMBOL := NA_character_]
  if (!"GENENAME" %in% names(x)) x[, GENENAME := NA_character_]
  if (!"ID" %in% names(x)) x[, ID := stable_id]
  if (!"ENTREZID" %in% names(x)) x[, ENTREZID := NA_character_]
  x <- x[!is.na(FEATURE_ID) & nzchar(FEATURE_ID)]
  x <- x[order(Pvalue, na.last = TRUE)]
  unique(x, by = "FEATURE_ID")
}

.compile_de_experiments <- function(DEGDatapath, SEPath) {
  files <- list.files(DEGDatapath, pattern = "[.]rds$", full.names = TRUE)
  files <- files[!grepl("_Raw|_Spotfire|_RPKMRaw|_CountRaw", basename(files))]
  if (!length(files)) {
    stop("No differential-expression RDS files were found in `DEGDatapath`.", call. = FALSE)
  }
  results <- lapply(files, function(path) {
    .standardize_de_result(readRDS(path), basename(path))
  })
  names(results) <- tools::file_path_sans_ext(basename(files))
  row_annotations <- data.table::rbindlist(
    lapply(results, function(x) {
      x[, .(FEATURE_ID, SPECIES, ID, ENTREZID, SYMBOL, GENENAME)]
    }),
    use.names = TRUE, fill = TRUE
  )
  row_annotations <- unique(row_annotations, by = "FEATURE_ID")
  feature_ids <- sort(unique(row_annotations$FEATURE_ID))
  row_annotations <- row_annotations[match(feature_ids, FEATURE_ID)]
  row_data <- as.data.frame(row_annotations)
  rownames(row_data) <- row_data$FEATURE_ID
  assays <- lapply(results, function(x) {
    values <- matrix(
      NA_real_, nrow = length(feature_ids), ncol = 3L,
      dimnames = list(feature_ids, c("logFC", "Pvalue", "AdjPValue"))
    )
    matched <- match(x$FEATURE_ID, feature_ids)
    values[matched, ] <- as.matrix(x[, .(logFC, Pvalue, AdjPValue)])
    values
  })
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = S4Vectors::SimpleList(assays),
    rowData = row_data,
    metadata = list(
      feature_key = "SPECIES:stable_identifier",
      sources = normalizePath(files, winslash = "/", mustWork = FALSE)
    )
  )
  directory <- dirname(SEPath)
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  saveRDS(se, SEPath)
  se
}

#' Validate a dataset manifest
#'
#' Checks the dataset-level configuration used by the batch-processing workflow.
#' The function validates identifiers, technologies, species, comparison codes,
#' and logical processing flags without downloading any data.
#'
#' @param manifest A data frame or data table containing one row per dataset or
#'   dataset comparison.
#'
#' @return The validated manifest as a data frame, with normalized `Technology`
#'   and `Species` values.
#' @examples
#' manifest <- data.frame(
#'   ID = "GSE00000",
#'   Technology = "RNA-seq",
#'   Species = "Homo sapiens",
#'   GEO2R = FALSE
#' )
#' ValidateDatasetManifest(manifest)
#' @family setup functions
#' @export
ValidateDatasetManifest <- function(manifest) {
  manifest <- as.data.frame(manifest, stringsAsFactors = FALSE)
  required <- c("ID", "Technology", "Species")
  missing <- setdiff(required, names(manifest))
  if (length(missing)) {
    stop(
      "Manifest is missing required columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (!nrow(manifest)) {
    stop("Manifest contains no datasets.", call. = FALSE)
  }
  if (anyNA(manifest$ID) || any(!nzchar(trimws(manifest$ID)))) {
    stop("Every manifest row must have a non-empty `ID`.", call. = FALSE)
  }
  technology <- tolower(trimws(manifest$Technology))
  technology[technology %in% c("rnaseq", "rna-seq", "rna seq")] <- "RNAseq"
  technology[technology %in% c("array", "microarray")] <- "Array"
  technology[technology == "mass spectrometry"] <- "Mass Spectrometry"
  technology[technology == "bisulfite sequencing"] <- "Bisulfite Sequencing"
  invalid_technology <- is.na(technology) | !technology %in% c(
    "RNAseq", "Array", "Mass Spectrometry", "Bisulfite Sequencing"
  )
  if (any(invalid_technology)) {
    stop(
      "Unsupported Technology values: ",
      paste(unique(manifest$Technology[invalid_technology]), collapse = ", "),
      call. = FALSE
    )
  }
  species <- .normalize_species(manifest$Species)
  if (anyNA(species)) {
    stop(
      "Unsupported Species values: ",
      paste(unique(manifest$Species[is.na(species)]), collapse = ", "),
      call. = FALSE
    )
  }
  manifest$Technology <- technology
  manifest$Species <- species
  logical_columns <- intersect(
    c("DownloadMetaData", "DownloadRawData", "GEO2R"), names(manifest)
  )
  for (column in logical_columns) {
    values <- manifest[[column]]
    if (any(!is.na(values) & !values %in% c(TRUE, FALSE, 0, 1))) {
      stop("`", column, "` must contain only TRUE, FALSE, or NA.", call. = FALSE)
    }
    manifest[[column]] <- as.logical(values)
  }
  if ("GEO2R" %in% names(manifest) && "ComparisonVector" %in% names(manifest)) {
    selected <- which(!is.na(manifest$GEO2R) & manifest$GEO2R)
    missing_comparison <- selected[
      is.na(manifest$ComparisonVector[selected]) |
        !nzchar(trimws(manifest$ComparisonVector[selected]))
    ]
    if (length(missing_comparison)) {
      stop(
        "GEO2R rows require a ComparisonVector. Missing for: ",
        paste(manifest$ID[missing_comparison], collapse = ", "),
        call. = FALSE
      )
    }
    invalid_comparison <- selected[
      !grepl("^[[:space:]01Xx,;]+$", manifest$ComparisonVector[selected])
    ]
    if (length(invalid_comparison)) {
      stop(
        "ComparisonVector may contain only 0, 1, X, whitespace, commas, and semicolons. Invalid for: ",
        paste(manifest$ID[invalid_comparison], collapse = ", "),
        call. = FALSE
      )
    }
  }
  manifest
}
