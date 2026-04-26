args_all_prepare <- commandArgs(trailingOnly = FALSE)
file_arg_prepare <- grep("^--file=", args_all_prepare, value = TRUE)
repo_root_prepare <- Sys.getenv("CSCN_REPO_ROOT", unset = "")
if (nzchar(repo_root_prepare)) {
  source(file.path(repo_root_prepare, "scripts", "prep", "scp259_utils.R"), chdir = FALSE)
} else {
  script_path_prepare <- normalizePath(sub("^--file=", "", file_arg_prepare[[1L]]), winslash = "/", mustWork = TRUE)
  source(file.path(dirname(script_path_prepare), "scp259_utils.R"), chdir = FALSE)
}

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir_scp259()
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
      "Rscript prepare_SCP259.R [options]",
      "",
      "Options:",
      sprintf("  --data-dir PATH   Dataset directory. Default: %s", default_data_dir_scp259()),
      "  --help            Show this message",
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
metadata_path <- file.path(data_dir, "metadata", "all.meta2.txt")
count_matrix_path <- file.path(output_dir, sprintf("count_matrix_%s_by_subject.csv", RUN_SLUG_SCP259))
metadata_out_path <- file.path(output_dir, sprintf("metadata_%s_by_subject.csv", RUN_SLUG_SCP259))
covariates_path <- file.path(output_dir, sprintf("covariates_%s_selected_cells.csv", RUN_SLUG_SCP259))

log_stage_scp259("Configuration")
log_scp259(sprintf("dataset dir: %s", data_dir))
log_scp259(sprintf("output dir: %s", output_dir))
validate_required_file_scp259(metadata_path, "metadata table")

log_stage_scp259("Select Cells")
all_meta <- read_scp259_metadata(metadata_path)
selected_meta <- filter_scp259_epithelial(
  metadata_df = all_meta,
  health_values = c("Healthy", "Inflamed"),
  clusters = CRYPT_PROLIF_CLUSTERS_SCP259
)
selected_meta$NAME <- as.character(selected_meta$NAME)
selected_meta$Subject <- as.character(selected_meta$Subject)
selected_meta$Health <- as.character(selected_meta$Health)
selected_meta$Sample <- as.character(selected_meta$Sample)
selected_meta$Cluster <- as.character(selected_meta$Cluster)
log_scp259(sprintf("selected healthy cells: %d", sum(selected_meta$Health == "Healthy")))
log_scp259(sprintf("selected inflamed cells: %d", sum(selected_meta$Health == "Inflamed")))

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(
  selected_meta[, c("NAME", "Subject", "Health", "Location", "Sample", "Cluster")],
  covariates_path,
  row.names = FALSE
)
log_scp259(sprintf("saved selected covariates: %s", covariates_path))

log_stage_scp259("Read Matrix")
epi <- read_scp259_epi_matrix(data_dir)

log_stage_scp259("Build Subject Pseudobulk")
pseudo <- subject_pseudobulk_scp259(epi$matrix, selected_meta)
count_df <- data.frame(gene = rownames(pseudo$count_mat), pseudo$count_mat, check.names = FALSE)
write.csv(count_df, count_matrix_path, row.names = FALSE)

meta_by_subject <- aggregate(
  NAME ~ Subject + Health,
  data = selected_meta,
  FUN = length
)
colnames(meta_by_subject)[ncol(meta_by_subject)] <- "n_cells"
write.csv(meta_by_subject, metadata_out_path, row.names = FALSE)
log_scp259(sprintf("saved count matrix: %s", count_matrix_path))
log_scp259(sprintf("saved metadata: %s", metadata_out_path))
log_scp259(sprintf("subjects kept: %d", nrow(meta_by_subject)))
log_scp259(sprintf("genes written: %d", nrow(count_df)))
