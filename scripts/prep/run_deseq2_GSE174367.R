options(stringsAsFactors = FALSE)

normalize_slug <- function(value) {
  slug <- tolower(gsub("[^A-Za-z0-9]+", "_", value))
  slug <- gsub("^_+|_+$", "", slug)
  slug
}

build_run_slug <- function(cell_type, case_group, control_group) {
  sprintf(
    "%s_%s_vs_%s",
    normalize_slug(cell_type),
    normalize_slug(case_group),
    normalize_slug(control_group)
  )
}

resolve_repo_root <- function() {
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

REPO_ROOT <- resolve_repo_root()

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop(
      "Package 'DESeq2' is required. Install it with BiocManager::install('DESeq2').",
      call. = FALSE
    )
  }
  library(DESeq2)
})

default_data_dir <- normalizePath(file.path(REPO_ROOT, "data", "GSE174367"), winslash = "/", mustWork = FALSE)

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir,
    cell_type = "ODC",
    case_group = "AD",
    control_group = "Control",
    top_n = 150L,
    min_count = 10L,
    min_samples = 3L,
    out_prefix = NULL
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
    } else if (key == "cell-type") {
      parsed$cell_type <- value
    } else if (key == "case-group") {
      parsed$case_group <- value
    } else if (key == "control-group") {
      parsed$control_group <- value
    } else if (key == "top-n") {
      parsed$top_n <- as.integer(value)
    } else if (key == "min-count") {
      parsed$min_count <- as.integer(value)
    } else if (key == "min-samples") {
      parsed$min_samples <- as.integer(value)
    } else if (key == "out-prefix") {
      parsed$out_prefix <- value
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
      "Rscript run_deseq2_GSE174367.R [options]",
      "",
      "Options:",
      sprintf("  --data-dir PATH         Dataset directory. Default: %s", default_data_dir),
      "  --cell-type NAME        Cell.Type to analyze. Default: ODC",
      "  --case-group NAME       Case diagnosis label. Default: AD",
      "  --control-group NAME    Control diagnosis label. Default: Control",
      "  --top-n INT             Number of genes to save. Default: 150",
      "  --min-count INT         Minimum count threshold for filtering. Default: 10",
      "  --min-samples INT       Minimum samples passing min-count. Default: 3",
      "  --out-prefix PREFIX     Optional output prefix. Default: deseq2_<celltype>_<case>_vs_<control>",
      "  --help                  Show this message",
      sep = "\n"
    )
  )
}

ensure_required_columns <- function(df, columns, label) {
  missing <- setdiff(columns, colnames(df))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", label, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
}

resolve_design_terms <- function(metadata) {
  terms <- character(0)

  if ("batch" %in% colnames(metadata) && length(unique(metadata$batch)) > 1L) {
    metadata$batch <- factor(metadata$batch)
    terms <- c(terms, "batch")
  }

  if ("sex" %in% colnames(metadata)) {
    non_missing <- !is.na(metadata$sex) & metadata$sex != ""
    if (all(non_missing) && length(unique(metadata$sex)) > 1L) {
      metadata$sex <- factor(metadata$sex)
      terms <- c(terms, "sex")
    }
  }

  if ("age" %in% colnames(metadata)) {
    metadata$age <- as.numeric(metadata$age)
    if (all(!is.na(metadata$age)) && length(unique(metadata$age)) > 1L) {
      terms <- c(terms, "age")
    }
  }

  metadata$diagnosis <- factor(metadata$diagnosis)
  terms <- c(terms, "diagnosis")
  list(metadata = metadata, terms = terms)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  print_help()
  quit(save = "no", status = 0L)
}

if (args$case_group == args$control_group) {
  stop("--case-group and --control-group must be different", call. = FALSE)
}

if (is.na(args$top_n) || args$top_n <= 0L) {
  stop("--top-n must be a positive integer", call. = FALSE)
}

if (is.na(args$min_count) || args$min_count <= 0L) {
  stop("--min-count must be a positive integer", call. = FALSE)
}

if (is.na(args$min_samples) || args$min_samples <= 0L) {
  stop("--min-samples must be a positive integer", call. = FALSE)
}

run_slug <- build_run_slug(args$cell_type, args$case_group, args$control_group)
data_dir <- normalizePath(args$data_dir, winslash = "/", mustWork = TRUE)
output_dir <- file.path(data_dir, "output_deseq")
count_path <- file.path(output_dir, sprintf("count_matrix_%s_by_sample.csv", run_slug))
metadata_path <- file.path(output_dir, sprintf("metadata_%s_by_sample.csv", run_slug))
out_prefix <- args$out_prefix
if (is.null(out_prefix) || identical(out_prefix, "")) {
  out_prefix <- sprintf("deseq2_%s", run_slug)
}
top_path <- file.path(output_dir, sprintf("%s_top%d_genes.csv", out_prefix, args$top_n))
all_path <- file.path(output_dir, sprintf("%s_all_genes.csv", out_prefix))

if (!file.exists(count_path)) {
  stop(sprintf("Missing count matrix: %s", count_path), call. = FALSE)
}
if (!file.exists(metadata_path)) {
  stop(sprintf("Missing metadata table: %s", metadata_path), call. = FALSE)
}

count_df <- read.csv(count_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)

ensure_required_columns(count_df, "gene", "count matrix")
ensure_required_columns(metadata, c("sample", "diagnosis", "n_cells"), "metadata")

metadata <- metadata[metadata$diagnosis %in% c(args$control_group, args$case_group), , drop = FALSE]
if (nrow(metadata) == 0L) {
  stop("No metadata rows remain after filtering to the requested diagnosis groups", call. = FALSE)
}

sample_ids <- unique(metadata$sample)
missing_samples <- setdiff(sample_ids, colnames(count_df))
if (length(missing_samples) > 0L) {
  stop(
    sprintf("Count matrix is missing metadata samples: %s", paste(missing_samples, collapse = ", ")),
    call. = FALSE
  )
}

count_df <- count_df[, c("gene", sample_ids), drop = FALSE]
count_mat <- as.matrix(count_df[, -1, drop = FALSE])
rownames(count_mat) <- count_df$gene
storage.mode(count_mat) <- "integer"

metadata <- metadata[match(sample_ids, metadata$sample), , drop = FALSE]
rownames(metadata) <- metadata$sample
metadata$diagnosis <- factor(metadata$diagnosis, levels = c(args$control_group, args$case_group))
design_resolution <- resolve_design_terms(metadata)
metadata <- design_resolution$metadata
design_formula <- as.formula(sprintf("~ %s", paste(design_resolution$terms, collapse = " + ")))

cat("=== Configuration ===\n")
cat("data_dir:", data_dir, "\n")
cat("cell_type:", args$cell_type, "\n")
cat("count matrix:", count_path, "\n")
cat("metadata:", metadata_path, "\n")
cat("case group:", args$case_group, "\n")
cat("control group:", args$control_group, "\n")
cat("top_n:", args$top_n, "\n")
cat("min_count:", args$min_count, "\n")
cat("min_samples:", args$min_samples, "\n")
cat("design:", deparse(design_formula), "\n")
cat("top output:", top_path, "\n")
cat("full output:", all_path, "\n\n")

cat("=== Samples ===\n")
print(metadata)
cat("\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = design_formula
)

keep <- rowSums(counts(dds) >= args$min_count) >= args$min_samples
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n\n")

dds <- DESeq(dds, sfType = "poscounts")
res <- results(dds, contrast = c("diagnosis", args$case_group, args$control_group))

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[!is.na(res_df$padj), , drop = FALSE]
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), , drop = FALSE]
top_df <- head(res_df, args$top_n)

write.csv(res_df, all_path, row.names = FALSE)
write.csv(top_df, top_path, row.names = FALSE)

cat("Wrote full DESeq2 results to:", all_path, "\n")
cat("Wrote top gene list to:", top_path, "\n")
cat("Top gene rows:", nrow(top_df), "\n\n")
print(head(top_df, 10L))
