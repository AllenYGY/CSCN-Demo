args_all_expr <- commandArgs(trailingOnly = FALSE)
file_arg_expr <- grep("^--file=", args_all_expr, value = TRUE)
repo_root_expr <- Sys.getenv("CSCN_REPO_ROOT", unset = "")
if (nzchar(repo_root_expr)) {
  source(file.path(repo_root_expr, "scripts", "prep", "scp259_utils.R"), chdir = FALSE)
} else {
  script_path_expr <- normalizePath(sub("^--file=", "", file_arg_expr[[1L]]), winslash = "/", mustWork = TRUE)
  source(file.path(dirname(script_path_expr), "..", "prep", "scp259_utils.R"), chdir = FALSE)
}

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir_scp259(),
    gene_list_path = NULL,
    sample_size = 4000L,
    random_seed = 42L,
    run_name = sprintf("SCP259_%s", RUN_SLUG_SCP259)
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

print_help <- function() {
  cat(
    paste(
      "Usage:",
      "Rscript prepare_SCP259_expression.R [options]",
      "",
      "Options:",
      sprintf("  --data-dir PATH        Dataset directory. Default: %s", default_data_dir_scp259()),
      "  --gene-list-path PATH  Top-gene CSV from DESeq2 (required).",
      "  --sample-size INT      Number of cells per group. Default: 4000",
      "  --random-seed INT      Random seed. Default: 42",
      "  --run-name NAME        Output prefix. Default: SCP259_inflamed_vs_healthy_crypt_prolif_epi",
      "  --help                 Show this message",
      sep = "\n"
    )
  )
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  print_help()
  quit(save = "no", status = 0L)
}

if (is.null(args$gene_list_path) || identical(args$gene_list_path, "")) {
  stop("--gene-list-path is required", call. = FALSE)
}

data_dir <- normalizePath(args$data_dir, winslash = "/", mustWork = TRUE)
gene_list_path <- normalizePath(args$gene_list_path, winslash = "/", mustWork = TRUE)
output_dir <- file.path(data_dir, "output_deseq")
metadata_path <- file.path(data_dir, "metadata", "all.meta2.txt")
epi_matrix_path <- file.path(data_dir, "expression", "5cdc540d328cee7a2efc2348", "gene_sorted-Epi.matrix.mtx")
epi_genes_path <- file.path(data_dir, "expression", "5cdc540d328cee7a2efc2348", "Epi.genes.tsv")
epi_barcodes_path <- file.path(data_dir, "expression", "5cdc540d328cee7a2efc2348", "Epi.barcodes2.tsv")

validate_required_file_scp259(metadata_path, "metadata table")
validate_required_file_scp259(gene_list_path, "top-gene list")
validate_required_file_scp259(epi_matrix_path, "Epi matrix")
validate_required_file_scp259(epi_genes_path, "Epi genes")
validate_required_file_scp259(epi_barcodes_path, "Epi barcodes")

top_gene_df <- read.csv(gene_list_path, stringsAsFactors = FALSE)
if (!"gene" %in% colnames(top_gene_df)) {
  stop(sprintf("Missing 'gene' column in %s", gene_list_path), call. = FALSE)
}
top_genes <- unique(top_gene_df$gene[!is.na(top_gene_df$gene) & top_gene_df$gene != ""])
if (length(top_genes) == 0L) {
  stop("Top-gene list is empty", call. = FALSE)
}

set.seed(args$random_seed)

meta <- read_scp259_metadata(metadata_path)
selected_meta <- filter_scp259_epithelial(
  metadata_df = meta,
  health_values = c("Healthy", "Inflamed"),
  clusters = CRYPT_PROLIF_CLUSTERS_SCP259
)
selected_meta$NAME <- as.character(selected_meta$NAME)
selected_meta$Health <- as.character(selected_meta$Health)

sampled <- sample_cells_by_health(selected_meta, sample_size = args$sample_size, random_seed = args$random_seed)

epi <- read_scp259_epi_matrix(data_dir)
gene_idx <- match(top_genes, rownames(epi$matrix))
gene_idx <- gene_idx[!is.na(gene_idx)]
used_genes <- rownames(epi$matrix)[gene_idx]
used_genes_path <- file.path(output_dir, sprintf("%s_top%d_genes_used.csv", args$run_name, length(used_genes)))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (group_name in names(sampled)) {
  group_df <- sampled[[group_name]]
  col_idx <- match(group_df$NAME, colnames(epi$matrix))
  if (anyNA(col_idx)) {
    missing <- head(group_df$NAME[is.na(col_idx)], 5)
    stop(sprintf("Missing sampled cells in matrix for %s: %s", group_name, paste(missing, collapse = ", ")), call. = FALSE)
  }

  small <- epi$matrix[gene_idx, col_idx, drop = FALSE]
  small <- as.matrix(small)
  totals <- colSums(small)
  totals[totals == 0] <- 1
  normalized <- sweep(small, 2, totals, "/") * 1e6
  normalized <- log1p(t(normalized))

  out_df <- data.frame(cell_id = group_df$NAME, normalized, check.names = FALSE)
  expr_path <- file.path(output_dir, sprintf("%s_%s_expression.csv.gz", args$run_name, tolower(group_name)))
  gz_con <- gzfile(expr_path, open = "wt")
  write.csv(out_df, gz_con, row.names = FALSE)
  close(gz_con)
  log_scp259(sprintf("saved %s expression matrix: %s", group_name, expr_path))

  sampled_path <- file.path(output_dir, sprintf("%s_%s_sampled_cells.csv", args$run_name, tolower(group_name)))
  write.csv(data.frame(cell_id = group_df$NAME, stringsAsFactors = FALSE), sampled_path, row.names = FALSE)
  log_scp259(sprintf("saved %s sampled cells: %s", group_name, sampled_path))
}

write.csv(data.frame(gene = used_genes, stringsAsFactors = FALSE), used_genes_path, row.names = FALSE)
log_scp259(sprintf("saved used-gene list: %s", used_genes_path))
