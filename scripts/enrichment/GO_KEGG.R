# ======================================
#     人类基因 GO / KEGG 富集分析脚本 (修复版)
#     数据来源: CSV 文件中的 gene 列
# ======================================

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

parse_args <- function(args) {
  parsed <- list(
    dataset = "GSE121893",
    data_dir = NULL,
    biomarker_file = NULL,
    used_gene_file = NULL,
    run_name = NULL,
    help = FALSE
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
    if (key == "dataset") {
      parsed$dataset <- value
    } else if (key == "data-dir") {
      parsed$data_dir <- value
    } else if (key == "biomarker-file") {
      parsed$biomarker_file <- value
    } else if (key == "used-gene-file") {
      parsed$used_gene_file <- value
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
      "Rscript GO_KEGG.R [options]",
      "",
      "Options:",
      "  --dataset NAME         Dataset key. Default: GSE121893",
      "  --data-dir PATH        Dataset directory. Default: <repo>/data/<dataset>",
      "  --biomarker-file PATH  Optional explicit biomarker CSV path",
      "  --used-gene-file PATH  Optional explicit top-genes-used CSV path",
      "  --run-name NAME        Output prefix. Default: dataset-specific preset",
      "  --help                 Show this message",
      sep = "\n"
    )
  )
}

resolve_defaults <- function(dataset, data_dir) {
  if (dataset == "GSE121893") {
    return(list(
      biomarker_file = file.path(data_dir, "Biomarkers_GSE121893_dhf_vs_n_all_all.csv"),
      used_gene_file = file.path(data_dir, "output_deseq", "GSE121893_dhf_vs_n_all_all_top150_genes_used.csv"),
      run_name = "GSE121893_dhf_vs_n_all_all"
    ))
  }

  if (dataset == "GSE138852") {
    return(list(
      biomarker_file = file.path(data_dir, "Biomarkers.csv"),
      used_gene_file = file.path(data_dir, "output_deseq", "GSE138852_top150_genes_used.csv"),
      run_name = "GSE138852"
    ))
  }

  return(list(
    biomarker_file = file.path(data_dir, "Biomarkers.csv"),
    used_gene_file = file.path(data_dir, "output_deseq", sprintf("%s_top150_genes_used.csv", dataset)),
    run_name = dataset
  ))
}

# ========= 1. 参数 =========
args <- parse_args(commandArgs(trailingOnly = TRUE))
if (isTRUE(args$help)) {
  print_help()
  quit(save = "no", status = 0L)
}

dataset_name <- args$dataset
dataset_dir <- if (is.null(args$data_dir) || identical(args$data_dir, "")) {
  file.path(REPO_ROOT, "data", dataset_name)
} else {
  normalizePath(args$data_dir, winslash = "/", mustWork = FALSE)
}

defaults <- resolve_defaults(dataset_name, dataset_dir)
run_name <- if (is.null(args$run_name) || identical(args$run_name, "")) defaults$run_name else args$run_name
biomarker_file_path <- if (is.null(args$biomarker_file) || identical(args$biomarker_file, "")) {
  defaults$biomarker_file
} else {
  args$biomarker_file
}
used_gene_file_path <- if (is.null(args$used_gene_file) || identical(args$used_gene_file, "")) {
  defaults$used_gene_file
} else {
  args$used_gene_file
}

# ========= 2. 安装 & 加载 R 包 =========
required_bioc_packages <- c(
  "org.Hs.eg.db", "clusterProfiler", "enrichplot",
  "DOSE", "ComplexHeatmap", "GOplot"
)

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in required_bioc_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    message(paste("正在安装", pkg, "..."))
    BiocManager::install(pkg, ask = FALSE)
    library(pkg, character.only = TRUE)
  }
}

cran_packages <- c("ggplot2", "stringr", "circlize", "dplyr", "tidyr", "ggnewscale")
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# ========= 3. 输入路径 & 基因集准备 =========
current_script_dir <- REPO_ROOT
input_base_dir <- file.path(dataset_dir, "enrichment_inputs")
output_base_dir <- file.path(dataset_dir, "enrichment_results")

if (!dir.exists(input_base_dir)) {
  dir.create(input_base_dir, recursive = TRUE)
  message(paste0("已创建输入目录: ", input_base_dir))
}

if (!dir.exists(output_base_dir)) {
  dir.create(output_base_dir, recursive = TRUE)
  message(paste0("已创建结果目录: ", output_base_dir))
}

GO_database <- "org.Hs.eg.db" # 人类 GO 数据库
KEGG_database <- "hsa" # 人类 KEGG 代码

unique_preserve_order <- function(values) {
  values[!duplicated(values)]
}

read_gene_column <- function(csv_path) {
  if (!file.exists(csv_path)) {
    stop(paste("❌ 文件未找到:", csv_path))
  }

  gene_df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (!"gene" %in% colnames(gene_df)) {
    stop(paste("❌ CSV 文件中未找到 'gene' 列:", csv_path))
  }

  genes <- toupper(na.omit(gene_df$gene))
  genes <- genes[genes != ""]
  unique_preserve_order(genes)
}

biomarker_genes <- read_gene_column(biomarker_file_path)
used_genes <- read_gene_column(used_gene_file_path)
deseq_only_genes <- used_genes[!used_genes %in% biomarker_genes]

message(paste("✅ Biomarker 基因数:", length(biomarker_genes)))
message(paste("✅ CSCN 实际使用 top 基因数:", length(used_genes)))
message(paste("✅ DESeq2_only 基因数:", length(deseq_only_genes)))
message(paste("✅ dataset:", dataset_name))
message(paste("✅ run name:", run_name))

if (length(biomarker_genes) == 0) {
  stop("❌ Biomarker 基因列表为空。")
}

if (length(deseq_only_genes) == 0) {
  stop("❌ DESeq2_only 基因列表为空。")
}

biomarker_export_path <- file.path(input_base_dir, "Biomarkers.csv")
deseq_only_export_path <- file.path(input_base_dir, "DESeq2_only_genes.csv")

write.csv(data.frame(gene = biomarker_genes), biomarker_export_path, row.names = FALSE)
write.csv(data.frame(gene = deseq_only_genes), deseq_only_export_path, row.names = FALSE)

message(paste("✅ 已导出 Biomarker 输入文件:", biomarker_export_path))
message(paste("✅ 已导出 DESeq2_only 输入文件:", deseq_only_export_path))

gene_sets <- list(
  Biomarkers = list(
    genes = biomarker_genes,
    output_dir = file.path(output_base_dir, "Biomarkers")
  ),
  DESeq2_only_genes = list(
    genes = deseq_only_genes,
    output_dir = file.path(output_base_dir, "DESeq2_only_genes")
  )
)

# ========= 4. 测试基因 ID 转换 =========
message("🔍 测试基因 ID 转换...")
all_input_genes <- unique_preserve_order(c(biomarker_genes, deseq_only_genes))
test_conversion <- tryCatch(
  {
    bitr(all_input_genes[1:min(5, length(all_input_genes))],
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = GO_database
    )
  },
  error = function(e) {
    stop(paste("基因 ID 转换失败:", e$message))
  }
)

if (nrow(test_conversion) > 0) {
  message(paste(
    "✅ 基因 ID 转换成功！测试的",
    nrow(test_conversion),
    "个基因转换通过。"
  ))
} else {
  stop("❌ 基因 ID 转换测试失败，没有基因成功转换。")
}

# ========= 5. 修复后的富集分析函数 =========
run_enrichment_analysis <- function(gene_list_vector, method_name, type_value, output_dir_for_run) {
  message(paste0("\n--- 🔄 正在处理: ", method_name, "_Gene_", type_value, " ---"))

  if (!dir.exists(output_dir_for_run)) {
    dir.create(output_dir_for_run, recursive = TRUE)
  }

  gene_list_vector <- na.omit(gene_list_vector)
  gene_list_vector <- gene_list_vector[gene_list_vector != ""]
  if (length(gene_list_vector) == 0) {
    return(NULL)
  }

  gene_entrez <- tryCatch(
    {
      bitr(gene_list_vector, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = GO_database)
    },
    error = function(e) {
      warning(paste0("❌ 基因 ID 转换失败: ", e$message))
      return(NULL)
    }
  )

  if (is.null(gene_entrez) || nrow(gene_entrez) == 0) {
    warning("⚠️ 没有基因成功转换为 ENTREZID，跳过。")
    return(NULL)
  } else {
    message(paste0("✅ 成功转换 ", nrow(gene_entrez), " 个基因 ID 为 ENTREZID。"))
  }

  file_prefix <- paste0(method_name, "_Gene_", type_value)

  # --- GO 富集 (修复版) ---
  message("▶ 进行 GO 富集分析...")
  GO_result <- enrichGO(gene_entrez$ENTREZID,
    OrgDb = GO_database,
    keyType = "ENTREZID",
    ont = "ALL",
    readable = TRUE,
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1
  )

  if (!is.null(GO_result) && nrow(as.data.frame(GO_result)) > 0) {
    write.csv(as.data.frame(GO_result),
      file = file.path(output_dir_for_run, paste0(file_prefix, "_GO_results.csv")),
      row.names = FALSE
    )

    show_num <- min(15, nrow(as.data.frame(GO_result)))

    # 处理长标签
    go_data <- as.data.frame(GO_result)
    max_desc_length <- max(nchar(go_data$Description))
    message(paste("最长的GO term描述字符数:", max_desc_length))

    # 截断过长的描述
    if (max_desc_length > 60) {
      go_data$Description <- ifelse(nchar(go_data$Description) > 60,
        paste0(substr(go_data$Description, 1, 57), "..."),
        go_data$Description
      )
      GO_result@result$Description <- go_data$Description
    }

    # 优化的 barplot
    tryCatch(
      {
        p_bar <- barplot(GO_result,
          x = "GeneRatio", split = "ONTOLOGY",
          showCategory = show_num,
          title = "GO Enrichment",
          label_format = 100
        ) +
          facet_grid(ONTOLOGY ~ ., space = "free_y", scales = "free_y") +
          theme_minimal() +
          theme(
            plot.margin = margin(20, 40, 20, 20, "pt"), # 增加边距
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            strip.text = element_text(size = 9),
            legend.text = element_text(size = 8)
          )

        ggsave(file.path(output_dir_for_run, paste0(file_prefix, "_GO_barplot.png")),
          plot = p_bar, width = 14, height = 12, dpi = 300
        )
        message("✅ GO barplot 已保存")
      },
      error = function(e) {
        message(paste("⚠️ GO barplot 生成失败:", e$message))
      }
    )

    # 优化的 dotplot
    tryCatch(
      {
        p_dot <- dotplot(GO_result,
                         x = "GeneRatio", split = "ONTOLOGY",
                         showCategory = show_num,
                         title = "GO Enrichment",
                         label_format = 100
        ) +
          facet_grid(ONTOLOGY ~ ., space = "free_y", scales = "free_y") +
          theme_minimal() +
          theme(
            plot.margin = margin(20, 40, 20, 20, "pt"), # 增加边距，特别是右边距
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.y = element_text(
              size = 10, 
              face = "bold",
              angle = 15,        # 旋转角度（度数）
              hjust = 1          # 水平对齐方式
            ),
            axis.text.x = element_text(size = 8),
            strip.text = element_text(size = 9),
            legend.position = "right",
            legend.text = element_text(size = 8),
            legend.title = element_text(size = 9),
            legend.key.size = grid::unit(0.8, "cm")
          )
        
        ggsave(file.path(output_dir_for_run, paste0(file_prefix, "_GO_dotplot.png")),
          plot = p_dot, width = 16, height = 12, dpi = 300
        ) # 增加宽度
        message("✅ GO dotplot 已保存")
      },
      error = function(e) {
        message(paste("⚠️ GO dotplot 生成失败:", e$message))
      }
    )
  } else {
    message("⚠️ GO 富集无显著通路。")
  }

  # --- KEGG 富集 (修复版) ---
  message("▶ 进行 KEGG 富集分析...")
  KEGG_result <- tryCatch(
    {
      enrichKEGG(gene_entrez$ENTREZID,
        organism = KEGG_database,
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.1
      )
    },
    error = function(e) {
      message(paste("⚠️ KEGG 富集分析失败:", e$message))
      return(NULL)
    }
  )

  if (!is.null(KEGG_result) && nrow(as.data.frame(KEGG_result)) > 0) {
    write.csv(as.data.frame(KEGG_result),
      file = file.path(output_dir_for_run, paste0(file_prefix, "_KEGG_results.csv")),
      row.names = FALSE
    )

    show_num <- min(15, nrow(as.data.frame(KEGG_result)))

    # 处理KEGG长标签
    kegg_data <- as.data.frame(KEGG_result)
    if (max(nchar(kegg_data$Description)) > 60) {
      kegg_data$Description <- ifelse(nchar(kegg_data$Description) > 60,
        paste0(substr(kegg_data$Description, 1, 57), "..."),
        kegg_data$Description
      )
      KEGG_result@result$Description <- kegg_data$Description
    }

    # 优化的 KEGG barplot
    tryCatch(
      {
        p_kegg_bar <- barplot(KEGG_result,
          showCategory = show_num,
          title = "KEGG Enrichment",
          label_format = 100
        ) +
          theme_minimal() +
          theme(
            plot.margin = margin(20, 40, 20, 20, "pt"),
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.y = element_text(size = 8),
            axis.text.x = element_text(size = 8)
          )

        ggsave(file.path(output_dir_for_run, paste0(file_prefix, "_KEGG_barplot.png")),
          plot = p_kegg_bar, width = 12, height = 8, dpi = 300
        )
        message("✅ KEGG barplot 已保存")
      },
      error = function(e) {
        message(paste("⚠️ KEGG barplot 生成失败:", e$message))
      }
    )

    # 优化的 KEGG dotplot
    tryCatch(
      {
        p_kegg_dot <- dotplot(KEGG_result,
                              showCategory = show_num,
                              title = "KEGG Enrichment",
                              label_format = 100
        ) +
          theme_minimal() +
          theme(
            plot.margin = margin(20, 40, 20, 20, "pt"),
            plot.title = element_text(size = 12, hjust = 0.5),
            axis.text.y = element_text(
              size = 15,         # 放大字体
              face = "bold",     # 加粗
              angle = 15,        # 旋转角度
              hjust = 1          # 水平对齐
            ),
            axis.text.x = element_text(size = 8),
            legend.position = "right",
            legend.text = element_text(size = 8)
          )
        

        ggsave(file.path(output_dir_for_run, paste0(file_prefix, "_KEGG_dotplot.png")),
          plot = p_kegg_dot, width = 14, height = 8, dpi = 300
        )
        message("✅ KEGG dotplot 已保存")
      },
      error = function(e) {
        message(paste("⚠️ KEGG dotplot 生成失败:", e$message))
      }
    )
  } else {
    message("⚠️ KEGG 富集无显著通路。")
  }

  message(paste0("--- ✅ 完成: ", file_prefix, " ---"))
  return(list(GO_result = GO_result, KEGG_result = KEGG_result, file_prefix = file_prefix))
}

# ========= 6. 执行分析 =========
collected_results <- list()
method_name_for_output <- run_name

for (gene_set_name in names(gene_sets)) {
  gene_set <- gene_sets[[gene_set_name]]
  result <- run_enrichment_analysis(
    gene_list_vector = gene_set$genes,
    method_name = method_name_for_output,
    type_value = gene_set_name,
    output_dir_for_run = gene_set$output_dir
  )

  if (!is.null(result)) {
    collected_results[[result$file_prefix]] <- result
  }
}

# ========= 7. 汇总报告 =========
message("\n========= 分析汇总 =========")
message(
  paste(
    "分析的基因数量:",
    paste(
      names(gene_sets),
      vapply(gene_sets, function(item) length(item$genes), integer(1)),
      sep = "=",
      collapse = ", "
    )
  )
)
message("输出目录:", output_base_dir)

if (length(collected_results) > 0) {
  for (result_name in names(collected_results)) {
    result <- collected_results[[result_name]]
    go_count <- ifelse(is.null(result$GO_result) || nrow(as.data.frame(result$GO_result)) == 0, 0,
      nrow(as.data.frame(result$GO_result))
    )
    kegg_count <- ifelse(is.null(result$KEGG_result) || nrow(as.data.frame(result$KEGG_result)) == 0, 0,
      nrow(as.data.frame(result$KEGG_result))
    )
    message(paste(result_name, "- GO通路:", go_count, "个, KEGG通路:", kegg_count, "个"))
  }
} else {
  message("⚠️ 未生成任何富集分析结果。")
}

message("\n✅ 所有人类基因富集分析已完成！")
message(paste("📂 结果目录:", output_base_dir))
