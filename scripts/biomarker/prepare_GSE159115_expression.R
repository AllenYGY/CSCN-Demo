options(stringsAsFactors = FALSE)

args_all_expr <- commandArgs(trailingOnly = FALSE)
file_arg_expr <- grep("^--file=", args_all_expr, value = TRUE)
repo_root_expr <- Sys.getenv("CSCN_REPO_ROOT", unset = "")
if (nzchar(repo_root_expr)) {
  source(file.path(repo_root_expr, "scripts", "prep", "gse159115_utils.R"), chdir = FALSE)
} else {
  script_path_expr <- normalizePath(sub("^--file=", "", file_arg_expr[[1L]]), winslash = "/", mustWork = TRUE)
  source(file.path(dirname(script_path_expr), "..", "prep", "gse159115_utils.R"), chdir = FALSE)
}

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir_gse159115(),
    gene_list_path = NULL,
    sample_size = 100L,
    random_seed = 42L,
    run_name = sprintf("%s_%s", DATA_SET_GSE159115, DEFAULT_RUN_SLUG_GSE159115)
  )
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected positional argument: %s", arg), call. = FALSE)
    }
    key <- substring(arg, 3L)
    if (i == length(args)) {
      stop(sprintf("Missing value for --%s", key), call. = FALSE)
    }
    value <- args[[i + 1L]]
    if (key == "data-dir") {
      parsed$data_dir <- value
    } else if (key == "gene-list-path") {
      parsed$gene_list_path <- value
    } else if (key == "sample-size") {
      parsed$sample_size <- as.integer(value)
    } else if (key == "random-seed") {
      parsed$random_seed <- as.integer(value)
    } else if (key == "run-name") {
      parsed$run_name <- value
    } else {
      stop(sprintf("Unknown argument: --%s", key), call. = FALSE)
    }
    i <- i + 2L
  }
  parsed
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

if (is.null(args$gene_list_path) || identical(args$gene_list_path, "")) {
  stop("--gene-list-path is required", call. = FALSE)
}

data_dir <- normalizePath(args$data_dir, winslash = "/", mustWork = TRUE)
gene_list_path <- normalizePath(args$gene_list_path, winslash = "/", mustWork = TRUE)
output_dir <- file.path(data_dir, "output_deseq")
raw_tar_path <- file.path(data_dir, "GSE159115_RAW.tar")

validate_required_file_gse159115(raw_tar_path, "RAW tar")
validate_required_file_gse159115(gene_list_path, "top-gene list")

top_gene_df <- read.csv(gene_list_path, stringsAsFactors = FALSE)
if (!"gene" %in% colnames(top_gene_df)) {
  stop(sprintf("Missing 'gene' column in %s", gene_list_path), call. = FALSE)
}
top_genes <- unique(top_gene_df$gene[!is.na(top_gene_df$gene) & top_gene_df$gene != ""])
if (length(top_genes) == 0L) {
  stop("Top-gene list is empty", call. = FALSE)
}

selection <- prepare_selected_annotations_gse159115(
  data_dir = data_dir,
  tumor_annos = DEFAULT_TUMOR_ANNOS_GSE159115,
  normal_annos = DEFAULT_NORMAL_ANNOS_GSE159115,
  paired_only = TRUE,
  min_cells_per_sample = 1L
)

eligible_df <- selection$selected_df
group_tables <- split(eligible_df, eligible_df$disease)

if (!all(c("normal", "tumor") %in% names(group_tables))) {
  stop("Expected both normal and tumor eligible cells", call. = FALSE)
}

set.seed(args$random_seed)
sampled_cells <- list()
for (group_name in c("normal", "tumor")) {
  group_df <- group_tables[[group_name]]
  if (nrow(group_df) < args$sample_size) {
    stop(
      sprintf("Group %s has only %d eligible cells, fewer than requested sample size %d", group_name, nrow(group_df), args$sample_size),
      call. = FALSE
    )
  }
  sampled_rows <- sort(sample(seq_len(nrow(group_df)), size = args$sample_size, replace = FALSE))
  sampled_cells[[group_name]] <- group_df[sampled_rows, , drop = FALSE]
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

tmp_info <- extract_h5_tar_gse159115(raw_tar_path)
sample_to_h5 <- tmp_info$sample_to_h5

normalized_matrices <- list()
for (group_name in c("normal", "tumor")) {
  group_df <- sampled_cells[[group_name]]
  split_by_sample <- split(group_df, group_df$sample)
  group_mat <- NULL
  for (sample_id in names(split_by_sample)) {
    h5_path <- sample_to_h5[[sample_id]]
    if (is.null(h5_path) || !file.exists(h5_path)) {
      stop(sprintf("No H5 file found for sample %s", sample_id), call. = FALSE)
    }
    extracted <- extract_topgene_counts_for_sample_gse159115(
      h5_path = h5_path,
      sample_id = sample_id,
      selected_cells_df = split_by_sample[[sample_id]],
      top_genes = top_genes
    )
    norm_mat <- normalize_log1p_gse159115(extracted$counts, extracted$totals)
    if (is.null(group_mat)) {
      group_mat <- norm_mat
    } else {
      group_mat <- rbind(group_mat, norm_mat)
    }
  }
  normalized_matrices[[group_name]] <- group_mat
}

combined <- rbind(normalized_matrices$normal, normalized_matrices$tumor)
keep_genes <- colSums(combined) > 0
used_genes <- colnames(combined)[keep_genes]

for (group_name in c("normal", "tumor")) {
  mat <- normalized_matrices[[group_name]][, keep_genes, drop = FALSE]
  expr_path <- file.path(output_dir, sprintf("%s_%s_expression.csv.gz", args$run_name, group_name))
  expr_df <- data.frame(cell_id = rownames(mat), mat, check.names = FALSE)
  gz_con <- gzfile(expr_path, open = "wt")
  write.csv(expr_df, gz_con, row.names = FALSE)
  close(gz_con)

  sampled_cells_path <- file.path(output_dir, sprintf("%s_%s_sampled_cells.csv", args$run_name, group_name))
  write.csv(data.frame(cell_id = rownames(mat), stringsAsFactors = FALSE), sampled_cells_path, row.names = FALSE)
  log_gse159115(sprintf("saved %s expression matrix: %s", group_name, expr_path))
}

used_genes_path <- file.path(output_dir, sprintf("%s_top%d_genes_used.csv", args$run_name, length(used_genes)))
write.csv(data.frame(gene = used_genes, stringsAsFactors = FALSE), used_genes_path, row.names = FALSE)
log_gse159115(sprintf("saved used-gene list: %s", used_genes_path))
