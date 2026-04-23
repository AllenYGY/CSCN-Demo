options(stringsAsFactors = FALSE)

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
RUN_SLUG <- "nlung_epithelial_vs_tlung_malignant"

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop(
      "Package 'DESeq2' is required. Install it with BiocManager::install('DESeq2').",
      call. = FALSE
    )
  }
  library(DESeq2)
})

default_data_dir <- normalizePath(file.path(REPO_ROOT, "data", "GSE131907"), winslash = "/", mustWork = FALSE)

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir,
    case_group = "cancer",
    control_group = "normal",
    top_n = 150L,
    min_count = 10L,
    min_samples = 2L,
    out_prefix = "deseq2_nlung_epithelial_vs_tlung_malignant"
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
      "Rscript run_deseq2_GSE131907.R [options]",
      "",
      "Options:",
      sprintf("  --data-dir PATH         Dataset directory. Default: %s/data/GSE131907", REPO_ROOT),
      "  --case-group NAME       Case disease label. Default: cancer",
      "  --control-group NAME    Control disease label. Default: normal",
      "  --top-n INT             Number of genes to save. Default: 150",
      "  --min-count INT         Minimum count threshold for filtering. Default: 10",
      "  --min-samples INT       Minimum samples passing min-count. Default: 2",
      "  --out-prefix PREFIX     Output prefix. Default: deseq2_nlung_epithelial_vs_tlung_malignant",
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

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  print_help()
  quit(save = "no", status = 0L)
}

if (args$case_group == args$control_group) {
  stop("--case-group and --control-group must be different", call. = FALSE)
}

data_dir <- normalizePath(args$data_dir, winslash = "/", mustWork = TRUE)
output_dir <- file.path(data_dir, "output_deseq")
count_path <- file.path(output_dir, sprintf("count_matrix_%s_by_sample.csv", RUN_SLUG))
metadata_path <- file.path(output_dir, sprintf("metadata_%s_by_sample.csv", RUN_SLUG))
top_path <- file.path(output_dir, sprintf("%s_top%d_genes.csv", args$out_prefix, args$top_n))
all_path <- file.path(output_dir, sprintf("%s_all_genes.csv", args$out_prefix))

if (!file.exists(count_path)) {
  stop(sprintf("Missing count matrix: %s", count_path), call. = FALSE)
}

if (!file.exists(metadata_path)) {
  stop(sprintf("Missing metadata table: %s", metadata_path), call. = FALSE)
}

cat("=== Configuration ===\n")
cat("data_dir:", data_dir, "\n")
cat("count matrix:", count_path, "\n")
cat("metadata:", metadata_path, "\n")
cat("case group:", args$case_group, "\n")
cat("control group:", args$control_group, "\n")
cat("top_n:", args$top_n, "\n")
cat("min_count:", args$min_count, "\n")
cat("min_samples:", args$min_samples, "\n")
cat("top output:", top_path, "\n")
cat("full output:", all_path, "\n\n")

count_df <- read.csv(count_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)

ensure_required_columns(count_df, "gene", "count matrix")
ensure_required_columns(metadata, c("sample", "disease"), "metadata")

metadata <- metadata[metadata$disease %in% c(args$control_group, args$case_group), , drop = FALSE]
if (nrow(metadata) == 0L) {
  stop("No metadata rows remain after filtering to the requested disease groups", call. = FALSE)
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
metadata$disease <- factor(metadata$disease, levels = c(args$control_group, args$case_group))

cat("=== Samples ===\n")
print(metadata)
cat("\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = ~ disease
)

keep <- rowSums(counts(dds) >= args$min_count) >= args$min_samples
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Design: ~ disease\n\n")

dds <- DESeq(dds, sfType = "poscounts")
res <- results(dds, contrast = c("disease", args$case_group, args$control_group))

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
