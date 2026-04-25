options(stringsAsFactors = FALSE)

args_all_deseq <- commandArgs(trailingOnly = FALSE)
file_arg_deseq <- grep("^--file=", args_all_deseq, value = TRUE)
repo_root_deseq <- Sys.getenv("CSCN_REPO_ROOT", unset = "")
if (nzchar(repo_root_deseq)) {
  source(file.path(repo_root_deseq, "scripts", "prep", "gse159115_utils.R"), chdir = FALSE)
} else {
  script_path_deseq <- normalizePath(sub("^--file=", "", file_arg_deseq[[1L]]), winslash = "/", mustWork = TRUE)
  source(file.path(dirname(script_path_deseq), "gse159115_utils.R"), chdir = FALSE)
}

suppressPackageStartupMessages({
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop(
      "Package 'DESeq2' is required. Install it with BiocManager::install('DESeq2').",
      call. = FALSE
    )
  }
  library(DESeq2)
})

run_slug <- DEFAULT_RUN_SLUG_GSE159115

parse_args <- function(args) {
  parsed <- list(
    data_dir = default_data_dir_gse159115(),
    case_group = "tumor",
    control_group = "normal",
    top_n = 150L,
    min_count = 10L,
    min_samples = 2L,
    out_prefix = "deseq2_ccrcc_tumor_vs_ptb_ptc_normal"
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
      "Rscript run_deseq2_GSE159115.R [options]",
      "",
      "Options:",
      sprintf("  --data-dir PATH         Dataset directory. Default: %s", default_data_dir_gse159115()),
      "  --case-group NAME       Case disease label. Default: tumor",
      "  --control-group NAME    Control disease label. Default: normal",
      "  --top-n INT             Number of genes to save. Default: 150",
      "  --min-count INT         Minimum count threshold for filtering. Default: 10",
      "  --min-samples INT       Minimum samples passing min-count. Default: 2",
      "  --out-prefix PREFIX     Output prefix. Default: deseq2_ccrcc_tumor_vs_ptb_ptc_normal",
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

data_dir <- normalizePath(args$data_dir, winslash = "/", mustWork = TRUE)
output_dir <- file.path(data_dir, "output_deseq")
count_path <- file.path(output_dir, sprintf("count_matrix_%s_paired_by_sample.csv", run_slug))
metadata_path <- file.path(output_dir, sprintf("metadata_%s_paired_by_sample.csv", run_slug))
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
cat("design: ~ patient + disease\n")
cat("top output:", top_path, "\n")
cat("full output:", all_path, "\n\n")

count_df <- read.csv(count_path, check.names = FALSE)
metadata <- read.csv(metadata_path, check.names = FALSE)

ensure_required_columns(count_df, "gene", "count matrix")
ensure_required_columns(metadata, c("sample", "patient", "disease"), "metadata")

sample_ids <- metadata$sample
count_df <- count_df[, c("gene", sample_ids), drop = FALSE]
count_mat <- as.matrix(count_df[, -1, drop = FALSE])
rownames(count_mat) <- count_df$gene
storage.mode(count_mat) <- "integer"

rownames(metadata) <- metadata$sample
metadata$disease <- factor(metadata$disease, levels = c(args$control_group, args$case_group))
metadata$patient <- factor(metadata$patient)

cat("=== Samples ===\n")
print(metadata)
cat("\n")

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = ~ patient + disease
)

keep <- rowSums(counts(dds) >= args$min_count) >= args$min_samples
dds <- dds[keep, ]

cat("Genes after filtering:", nrow(dds), "\n")
cat("Samples:", ncol(dds), "\n")
cat("Patients:", length(unique(metadata$patient)), "\n\n")

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
