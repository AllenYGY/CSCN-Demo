options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop("Package 'Matrix' is required.", call. = FALSE)
  }
  library(Matrix)
})

DATA_SET_SCP259 <- "SCP259"
RUN_SLUG_SCP259 <- "inflamed_vs_healthy_crypt_prolif_epi"
CRYPT_PROLIF_CLUSTERS_SCP259 <- c(
  "Stem",
  "Cycling TA",
  "TA 1",
  "TA 2",
  "Enterocyte Progenitors",
  "Secretory TA"
)
DIFF_EPITHELIAL_CLUSTERS_SCP259 <- c(
  "Best4+ Enterocytes",
  "Enterocytes",
  "Enteroendocrine",
  "Goblet",
  "Immature Enterocytes 1",
  "Immature Enterocytes 2",
  "Immature Goblet",
  "M cells",
  "Tuft"
)
ALL_EPI_CLUSTERS_SCP259 <- c(CRYPT_PROLIF_CLUSTERS_SCP259, DIFF_EPITHELIAL_CLUSTERS_SCP259)

resolve_repo_root_scp259 <- function() {
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

default_data_dir_scp259 <- function() {
  repo_root <- resolve_repo_root_scp259()
  candidates <- c(file.path(repo_root, "data", "scp259"), file.path(repo_root, "data", "SCP259"))
  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  normalizePath(candidates[[1L]], winslash = "/", mustWork = FALSE)
}

log_scp259 <- function(message) {
  cat(sprintf("[%s] %s\n", DATA_SET_SCP259, message))
}

log_stage_scp259 <- function(title) {
  cat("\n")
  cat(sprintf("=== %s ===\n", title))
}

validate_required_file_scp259 <- function(path, description) {
  if (!file.exists(path)) {
    stop(sprintf("Missing %s: %s", description, path), call. = FALSE)
  }
  log_scp259(sprintf("%s: %s", description, path))
}

read_scp259_metadata <- function(metadata_path) {
  df <- read.delim(metadata_path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("NAME", "Cluster", "Subject", "Health", "Location", "Sample")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0L) {
    stop(
      sprintf("Metadata file %s is missing columns: %s", metadata_path, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  df
}

filter_scp259_epithelial <- function(metadata_df, health_values, clusters = ALL_EPI_CLUSTERS_SCP259) {
  metadata_df <- metadata_df[metadata_df$Location == "Epi", , drop = FALSE]
  metadata_df <- metadata_df[metadata_df$Health %in% health_values, , drop = FALSE]
  metadata_df <- metadata_df[metadata_df$Cluster %in% clusters, , drop = FALSE]
  metadata_df
}

read_scp259_epi_matrix <- function(data_dir) {
  expr_dir <- file.path(data_dir, "expression", "5cdc540d328cee7a2efc2348")
  matrix_path <- file.path(expr_dir, "gene_sorted-Epi.matrix.mtx")
  genes_path <- file.path(expr_dir, "Epi.genes.tsv")
  barcodes_path <- file.path(expr_dir, "Epi.barcodes2.tsv")

  validate_required_file_scp259(matrix_path, "Epi matrix")
  validate_required_file_scp259(genes_path, "Epi genes")
  validate_required_file_scp259(barcodes_path, "Epi barcodes")

  genes <- readLines(genes_path, warn = FALSE)
  barcodes <- readLines(barcodes_path, warn = FALSE)
  mat <- readMM(matrix_path)

  if (nrow(mat) != length(genes) || ncol(mat) != length(barcodes)) {
    stop(
      sprintf(
        "Matrix dimensions %s do not match genes (%d) / barcodes (%d)",
        paste(dim(mat), collapse = " x "),
        length(genes),
        length(barcodes)
      ),
      call. = FALSE
    )
  }

  rownames(mat) <- genes
  colnames(mat) <- barcodes
  list(matrix = mat, genes = genes, barcodes = barcodes)
}

sample_cells_by_health <- function(metadata_df, sample_size, random_seed) {
  set.seed(random_seed)
  sampled <- list()
  for (health in c("Healthy", "Inflamed")) {
    subset_df <- metadata_df[metadata_df$Health == health, , drop = FALSE]
    if (nrow(subset_df) < sample_size) {
      stop(
        sprintf("Health group %s has only %d cells, fewer than sample-size %d", health, nrow(subset_df), sample_size),
        call. = FALSE
      )
    }
    sampled_idx <- sort(sample(seq_len(nrow(subset_df)), size = sample_size, replace = FALSE))
    sampled[[health]] <- subset_df[sampled_idx, , drop = FALSE]
  }
  sampled
}

subject_pseudobulk_scp259 <- function(mat, selected_df) {
  selected_df$Subject <- as.character(selected_df$Subject)
  selected_df$NAME <- as.character(selected_df$NAME)

  selected_col_idx <- match(selected_df$NAME, colnames(mat))
  if (anyNA(selected_col_idx)) {
    missing <- head(selected_df$NAME[is.na(selected_col_idx)], 5)
    stop(
      sprintf("Missing selected cells in matrix: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }

  subject_levels <- unique(selected_df$Subject)
  subj_map <- match(selected_df$Subject, subject_levels)
  counts <- vector("list", length(subject_levels))
  names(counts) <- subject_levels
  for (subject in subject_levels) {
    cols <- selected_col_idx[selected_df$Subject == subject]
    counts[[subject]] <- Matrix::rowSums(mat[, cols, drop = FALSE])
  }

  gene_names <- rownames(mat)
  count_mat <- do.call(cbind, counts)
  rownames(count_mat) <- gene_names
  list(count_mat = count_mat, subject_levels = subject_levels)
}

normalize_log1p_matrix_scp259 <- function(count_matrix, target_sum = 1e6) {
  totals <- rowSums(count_matrix)
  totals[totals == 0] <- 1
  normalized <- sweep(count_matrix, 1, totals, "/") * target_sum
  log1p(normalized)
}
