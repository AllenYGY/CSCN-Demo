args_all_prepare <- commandArgs(trailingOnly = FALSE)
file_arg_prepare <- grep("^--file=", args_all_prepare, value = TRUE)
repo_root_prepare <- Sys.getenv("CSCN_REPO_ROOT", unset = "")
if (nzchar(repo_root_prepare)) {
  source(file.path(repo_root_prepare, "scripts", "prep", "gse159115_utils.R"), chdir = FALSE)
} else {
  script_path_prepare <- normalizePath(sub("^--file=", "", file_arg_prepare[[1L]]), winslash = "/", mustWork = TRUE)
  source(file.path(dirname(script_path_prepare), "gse159115_utils.R"), chdir = FALSE)
}

repo_root <- resolve_repo_root_gse159115()
run_slug <- DEFAULT_RUN_SLUG_GSE159115

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir_gse159115(),
    min_cells_per_sample = 20L
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected positional argument: %s", arg), call. = FALSE)
    }
    key <- substring(arg, 3L)
    if (key %in% c("help", "h")) {
      parsed$help <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for --%s", key), call. = FALSE)
    }
    value <- args[[i + 1L]]
    if (key == "data-dir") {
      parsed$data_dir <- value
    } else if (key == "min-cells-per-sample") {
      parsed$min_cells_per_sample <- as.integer(value)
    } else {
      stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    }
    i <- i + 2L
  }
  parsed
}

print_help <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript prepare_GSE159115.R [options]",
      "",
      "Options:",
      sprintf("  --data-dir PATH              Dataset directory. Default: %s", default_data_dir_gse159115()),
      "  --min-cells-per-sample INT   Minimum cells per sample to keep. Default: 20",
      "  --help                       Show this message",
      sep = "\n"
    )
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  print_help()
  quit(save = "no", status = 0L)
}

data_dir <- normalizePath(args$data_dir, winslash = "/", mustWork = TRUE)
output_dir <- file.path(data_dir, "output_deseq")
raw_tar_path <- file.path(data_dir, "GSE159115_RAW.tar")
covariates_path <- file.path(data_dir, sprintf("%s_%s_covariates.csv.gz", DATA_SET_GSE159115, run_slug))
metadata_path <- file.path(output_dir, sprintf("metadata_%s_paired_by_sample.csv", run_slug))
count_matrix_path <- file.path(output_dir, sprintf("count_matrix_%s_paired_by_sample.csv", run_slug))
excluded_samples_path <- file.path(output_dir, sprintf("excluded_samples_%s_paired_by_sample.csv", run_slug))

log_stage_gse159115("Configuration")
log_gse159115(sprintf("dataset dir: %s", data_dir))
log_gse159115(sprintf("output dir: %s", output_dir))
log_gse159115(sprintf("run slug: %s", run_slug))
log_gse159115(sprintf("normal annos: %s", paste(DEFAULT_NORMAL_ANNOS_GSE159115, collapse = ", ")))
log_gse159115(sprintf("tumor annos: %s", paste(DEFAULT_TUMOR_ANNOS_GSE159115, collapse = ", ")))
log_gse159115(sprintf("min cells per sample: %d", args$min_cells_per_sample))
validate_required_file_gse159115(raw_tar_path, "RAW tar")

log_stage_gse159115("Select Cells")
selection <- prepare_selected_annotations_gse159115(
  data_dir = data_dir,
  tumor_annos = DEFAULT_TUMOR_ANNOS_GSE159115,
  normal_annos = DEFAULT_NORMAL_ANNOS_GSE159115,
  paired_only = TRUE,
  min_cells_per_sample = args$min_cells_per_sample
)

for (group_name in c("tumor", "normal")) {
  group_summary <- selection$summary[[group_name]]
  log_gse159115(sprintf("%s selected cells before pairing/filter: %d", group_name, group_summary$total_cells))
  log_gse159115(sprintf("%s sample counts: %s", group_name, paste(names(group_summary$sample_counts), as.integer(group_summary$sample_counts), sep = "=", collapse = ", ")))
  log_gse159115(sprintf("%s anno counts: %s", group_name, paste(names(group_summary$anno_counts), as.integer(group_summary$anno_counts), sep = "=", collapse = ", ")))
}

log_gse159115(sprintf("kept paired patients: %s", paste(selection$kept_patients, collapse = ", ")))
log_gse159115(sprintf("kept samples: %s", paste(selection$kept_samples, collapse = ", ")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(selection$excluded_samples, excluded_samples_path, row.names = FALSE)

covariates_df <- selection$selected_df[, c("cell", "sample", "patient", "anno", "disease", "group")]
gz_con <- gzfile(covariates_path, open = "wt")
write.csv(covariates_df, gz_con, row.names = FALSE)
close(gz_con)
log_gse159115(sprintf("saved covariates: %s", covariates_path))

metadata_df <- aggregate(cell ~ sample + patient + disease, data = selection$selected_df, FUN = length)
colnames(metadata_df)[ncol(metadata_df)] <- "n_cells"
anno_by_sample <- aggregate(anno ~ sample + patient + disease, data = selection$selected_df, FUN = function(x) paste(unique(x), collapse = "+"))
metadata_df <- merge(metadata_df, anno_by_sample, by = c("sample", "patient", "disease"), all.x = TRUE, sort = FALSE)
metadata_df <- metadata_df[order(metadata_df$patient, metadata_df$disease), ]
write.csv(metadata_df, metadata_path, row.names = FALSE)
log_gse159115(sprintf("saved metadata: %s", metadata_path))

log_stage_gse159115("Extract H5 Files")
tmp_info <- extract_h5_tar_gse159115(raw_tar_path)
sample_to_h5 <- tmp_info$sample_to_h5

log_stage_gse159115("Build Sample-Level Pseudobulk")
sample_counts_list <- list()
log_gse159115(sprintf("available H5 samples: %s", paste(sort(names(sample_to_h5)), collapse = ", ")))
for (sample_id in selection$kept_samples) {
  sample_cells_df <- selection$selected_df[selection$selected_df$sample == sample_id, , drop = FALSE]
  h5_index <- match(as.character(sample_id), names(sample_to_h5))
  if (is.na(h5_index)) {
    stop(sprintf("No H5 file found for sample %s", sample_id), call. = FALSE)
  }
  h5_path <- unname(sample_to_h5[h5_index])
  if (!file.exists(h5_path)) {
    stop(sprintf("Resolved H5 path does not exist for sample %s: %s", sample_id, h5_path), call. = FALSE)
  }
  sample_counts_list[[sample_id]] <- aggregate_sample_counts_gse159115(
    h5_path = h5_path,
    sample_id = sample_id,
    selected_cells_df = sample_cells_df
  )
  log_gse159115(sprintf("aggregated pseudobulk counts for %s (%d cells)", sample_id, nrow(sample_cells_df)))
}

all_genes <- unique(unlist(lapply(sample_counts_list, names), use.names = FALSE))
count_mat <- matrix(0, nrow = length(all_genes), ncol = length(sample_counts_list), dimnames = list(all_genes, names(sample_counts_list)))
for (sample_id in names(sample_counts_list)) {
  vec <- sample_counts_list[[sample_id]]
  count_mat[names(vec), sample_id] <- vec
}

count_df <- data.frame(gene = rownames(count_mat), count_mat, check.names = FALSE)
write.csv(count_df, count_matrix_path, row.names = FALSE)
log_gse159115(sprintf("saved count matrix: %s", count_matrix_path))
log_gse159115(sprintf("genes written: %d", nrow(count_df)))
