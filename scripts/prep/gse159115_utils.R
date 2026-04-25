options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("hdf5r", quietly = TRUE)) {
    stop(
      "Package 'hdf5r' is required. Install it with install.packages('hdf5r').",
      call. = FALSE
    )
  }
  library(hdf5r)
  library(Matrix)
})

DATA_SET_GSE159115 <- "GSE159115"
DEFAULT_RUN_SLUG_GSE159115 <- "ccrcc_tumor_vs_ptb_ptc_normal"
DEFAULT_TUMOR_ANNOS_GSE159115 <- c("Tumor")
DEFAULT_NORMAL_ANNOS_GSE159115 <- c("PT-B", "PT-C")

resolve_repo_root_gse159115 <- function() {
  env_root <- Sys.getenv("CSCN_REPO_ROOT", unset = "")
  if (nzchar(env_root)) {
    return(normalizePath(env_root, winslash = "/", mustWork = TRUE))
  }

  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    script_path <- normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = TRUE)
    return(normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE))
  }

  normalizePath(".", winslash = "/", mustWork = TRUE)
}

default_data_dir_gse159115 <- function() {
  repo_root <- resolve_repo_root_gse159115()
  candidates <- c(
    file.path(repo_root, "data", "GSE159115"),
    file.path(repo_root, "data", "gse159115")
  )
  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(candidates[[1L]], winslash = "/", mustWork = FALSE)
}

log_gse159115 <- function(message) {
  cat(sprintf("[%s] %s\n", DATA_SET_GSE159115, message))
}

log_stage_gse159115 <- function(title) {
  cat("\n")
  cat(sprintf("=== %s ===\n", title))
}

validate_required_file_gse159115 <- function(path, description) {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s: %s", description, path), call. = FALSE)
  }
  log_gse159115(sprintf("%s: %s", description, path))
}

parse_tar_list_name_to_sample_gse159115 <- function(member_name) {
  sub("^GSM[0-9]+_(SI_[^_]+)_filtered_gene_bc_matrices_h5\\.h5$", "\\1", basename(member_name))
}

extract_h5_tar_gse159115 <- function(raw_tar_path, out_dir = NULL) {
  if (is.null(out_dir)) {
    out_dir <- file.path(dirname(raw_tar_path), ".gse159115_h5_cache")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  h5_files <- list.files(out_dir, pattern = "\\.h5$", full.names = TRUE)
  if (length(h5_files) == 0L) {
    utils::untar(raw_tar_path, exdir = out_dir)
    h5_files <- list.files(out_dir, pattern = "\\.h5$", full.names = TRUE)
  }
  if (length(h5_files) == 0L) {
    stop(sprintf("No .h5 files found after extracting %s", raw_tar_path), call. = FALSE)
  }
  sample_to_h5 <- stats::setNames(h5_files, vapply(h5_files, function(path) {
    parse_tar_list_name_to_sample_gse159115(basename(path))
  }, character(1)))
  list(out_dir = out_dir, sample_to_h5 = sample_to_h5)
}

read_annotation_table_gse159115 <- function(path, disease, allowed_annos) {
  df <- read.csv(path, check.names = FALSE)
  required <- c("cell", "sample", "anno", "patient", "doublet")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0L) {
    stop(
      sprintf("Annotation table %s is missing columns: %s", path, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  df <- df[df$anno %in% allowed_annos & df$doublet == "FALSE", , drop = FALSE]
  df$disease <- disease
  df
}

prepare_selected_annotations_gse159115 <- function(
    data_dir,
    tumor_annos = DEFAULT_TUMOR_ANNOS_GSE159115,
    normal_annos = DEFAULT_NORMAL_ANNOS_GSE159115,
    paired_only = TRUE,
    min_cells_per_sample = 20L) {
  ccrcc_path <- file.path(data_dir, "GSE159115_ccRCC_anno.csv.gz")
  normal_path <- file.path(data_dir, "GSE159115_normal_anno.csv.gz")

  validate_required_file_gse159115(ccrcc_path, "ccRCC annotation")
  validate_required_file_gse159115(normal_path, "normal annotation")

  tumor_df <- read_annotation_table_gse159115(ccrcc_path, "tumor", tumor_annos)
  normal_df <- read_annotation_table_gse159115(normal_path, "normal", normal_annos)
  tumor_df$sample <- as.character(tumor_df$sample)
  tumor_df$patient <- as.character(tumor_df$patient)
  tumor_df$anno <- as.character(tumor_df$anno)
  normal_df$sample <- as.character(normal_df$sample)
  normal_df$patient <- as.character(normal_df$patient)
  normal_df$anno <- as.character(normal_df$anno)

  summary <- list(
    tumor = list(
      total_cells = nrow(tumor_df),
      sample_counts = sort(table(tumor_df$sample), decreasing = TRUE),
      patient_counts = sort(table(tumor_df$patient), decreasing = TRUE),
      anno_counts = sort(table(tumor_df$anno), decreasing = TRUE)
    ),
    normal = list(
      total_cells = nrow(normal_df),
      sample_counts = sort(table(normal_df$sample), decreasing = TRUE),
      patient_counts = sort(table(normal_df$patient), decreasing = TRUE),
      anno_counts = sort(table(normal_df$anno), decreasing = TRUE)
    )
  )

  selected_df <- rbind(tumor_df, normal_df)
  selected_df$group <- selected_df$disease

  selected_df$sample <- as.character(selected_df$sample)
  selected_df$patient <- as.character(selected_df$patient)
  selected_df$anno <- as.character(selected_df$anno)

  kept_patients <- sort(unique(selected_df$patient))
  if (isTRUE(paired_only)) {
    kept_patients <- sort(intersect(unique(tumor_df$patient), unique(normal_df$patient)))
    selected_df <- selected_df[selected_df$patient %in% kept_patients, , drop = FALSE]
  }

  sample_counts_df <- aggregate(
    cell ~ sample + patient + disease,
    data = selected_df,
    FUN = length
  )
  colnames(sample_counts_df)[ncol(sample_counts_df)] <- "n_cells"
  sample_counts_df$keep_after_min_cells <- sample_counts_df$n_cells >= min_cells_per_sample

  kept_samples_df <- sample_counts_df[sample_counts_df$keep_after_min_cells, , drop = FALSE]
  if (isTRUE(paired_only)) {
    keep_patient_table <- table(kept_samples_df$patient, kept_samples_df$disease)
    complete_patients <- rownames(keep_patient_table)[
      keep_patient_table[, "normal", drop = TRUE] > 0 &
      keep_patient_table[, "tumor", drop = TRUE] > 0
    ]
    kept_patients <- sort(complete_patients)
    kept_samples_df <- kept_samples_df[kept_samples_df$patient %in% kept_patients, , drop = FALSE]
  } else {
    kept_patients <- sort(unique(kept_samples_df$patient))
  }

  kept_samples <- sort(unique(as.character(kept_samples_df$sample)))
  filtered_df <- selected_df[selected_df$sample %in% kept_samples, , drop = FALSE]

  excluded_samples <- sample_counts_df[!sample_counts_df$sample %in% kept_samples, , drop = FALSE]
  if (nrow(excluded_samples) > 0L) {
    excluded_samples$reason <- ifelse(
      excluded_samples$n_cells < min_cells_per_sample,
      sprintf("n_cells<%d", min_cells_per_sample),
      "unpaired_after_filter"
    )
  } else {
    excluded_samples$reason <- character(0)
  }

  list(
    selected_df = filtered_df,
    tumor_df = tumor_df,
    normal_df = normal_df,
    summary = summary,
    sample_counts = sample_counts_df,
    kept_patients = kept_patients,
    kept_samples = kept_samples,
    excluded_samples = excluded_samples
  )
}

read_10x_h5_old_gse159115 <- function(h5_path) {
  h5 <- H5File$new(h5_path, mode = "r")
  on.exit(h5$close_all(), add = TRUE)
  root_name <- names(h5)[[1L]]
  grp <- h5[[root_name]]

  genes <- as.character(grp[["gene_names"]][])
  barcodes <- as.character(grp[["barcodes"]][])
  indices <- as.integer(grp[["indices"]][])
  indptr <- as.integer(grp[["indptr"]][])
  shape <- as.integer(grp[["shape"]][])
  data <- as.numeric(grp[["data"]][])

  mat <- sparseMatrix(
    i = indices + 1L,
    p = indptr,
    x = data,
    dims = shape,
    dimnames = list(genes, barcodes),
    index1 = TRUE
  )

  list(matrix = mat, genes = genes, barcodes = barcodes)
}

aggregate_gene_vector_gse159115 <- function(values, gene_names) {
  aggregated <- rowsum(matrix(as.numeric(values), ncol = 1), group = gene_names, reorder = FALSE)
  as.numeric(aggregated[, 1])
}

aggregate_sample_counts_gse159115 <- function(h5_path, sample_id, selected_cells_df) {
  sample_obj <- read_10x_h5_old_gse159115(h5_path)
  barcodes <- sub(paste0("^", sample_id, "_"), "", selected_cells_df$cell)
  matches <- match(barcodes, colnames(sample_obj$matrix))
  if (anyNA(matches)) {
    stop(
      sprintf(
        "Some selected cells were not found in H5 for sample %s: %s",
        sample_id,
        paste(head(selected_cells_df$cell[is.na(matches)], 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  gene_sums <- Matrix::rowSums(sample_obj$matrix[, matches, drop = FALSE])
  aggregated <- rowsum(matrix(as.numeric(gene_sums), ncol = 1), group = rownames(sample_obj$matrix), reorder = FALSE)
  setNames(as.numeric(aggregated[, 1]), rownames(aggregated))
}

extract_topgene_counts_for_sample_gse159115 <- function(h5_path, sample_id, selected_cells_df, top_genes) {
  sample_obj <- read_10x_h5_old_gse159115(h5_path)
  barcodes <- sub(paste0("^", sample_id, "_"), "", selected_cells_df$cell)
  matches <- match(barcodes, colnames(sample_obj$matrix))
  if (anyNA(matches)) {
    stop(
      sprintf(
        "Some selected cells were not found in H5 for sample %s: %s",
        sample_id,
        paste(head(selected_cells_df$cell[is.na(matches)], 5), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  totals <- as.numeric(Matrix::colSums(sample_obj$matrix[, matches, drop = FALSE]))

  row_idx <- which(rownames(sample_obj$matrix) %in% top_genes)
  small <- sample_obj$matrix[row_idx, matches, drop = FALSE]
  small_dense <- as.matrix(small)
  aggregated <- rowsum(small_dense, group = rownames(small), reorder = FALSE)

  out <- matrix(
    0,
    nrow = nrow(selected_cells_df),
    ncol = length(top_genes),
    dimnames = list(selected_cells_df$cell, top_genes)
  )
  present <- intersect(top_genes, rownames(aggregated))
  if (length(present) > 0L) {
    out[, present] <- t(aggregated[present, , drop = FALSE])
  }

  list(counts = out, totals = totals)
}

normalize_log1p_gse159115 <- function(count_matrix, totals, target_sum = 1e6) {
  totals[totals == 0] <- 1
  normalized <- sweep(count_matrix, 1, totals, "/") * target_sum
  log1p(normalized)
}
