# Time + celltype consensus GO & KEGG enrichment (minimal deps)

required_packages <- list(
  cran_packages = c("ggplot2", "stringr", "dplyr"),
  bioc_packages = c("clusterProfiler", "org.Hs.eg.db", "enrichplot")
)

check_and_install_minimal <- function(packages_list) {
  message("=== Checking minimal dependencies ===\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing BiocManager...")
    install.packages("BiocManager")
  }

  missing_cran <- c()
  for (pkg in packages_list$cran_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_cran <- c(missing_cran, pkg)
    }
  }
  if (length(missing_cran) > 0) {
    message(paste0("Installing CRAN packages: ", paste(missing_cran, collapse = ", ")))
    install.packages(missing_cran)
  }

  missing_bioc <- c()
  for (pkg in packages_list$bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_bioc <- c(missing_bioc, pkg)
    }
  }
  if (length(missing_bioc) > 0) {
    message(paste0("Installing Bioconductor packages: ", paste(missing_bioc, collapse = ", ")))
    BiocManager::install(missing_bioc, ask = FALSE)
  }

  critical_packages <- c("clusterProfiler", "org.Hs.eg.db")
  missing_critical <- c()
  for (pkg in critical_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_critical <- c(missing_critical, pkg)
    }
  }
  if (length(missing_critical) > 0) {
    message("\n❌ Missing critical packages:")
    message(paste0("Missing: ", paste(missing_critical, collapse = ", ")))
    message("\nPlease install manually:")
    message("BiocManager::install(c('clusterProfiler', 'org.Hs.eg.db'))")
    return(FALSE)
  }

  message("✅ Dependencies ready\n")
  return(TRUE)
}

packages_ok <- check_and_install_minimal(required_packages)
if (!packages_ok) stop("Missing required packages")

library(stringr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

ggplot_available <- requireNamespace("ggplot2", quietly = TRUE)
enrichplot_available <- requireNamespace("enrichplot", quietly = TRUE)
if (ggplot_available) library(ggplot2)
if (enrichplot_available) library(enrichplot)

clean_gene_id <- function(gene_id) {
  cleaned <- str_trim(gene_id)
  cleaned <- str_split(cleaned, "\t", simplify = TRUE)[,1]
  cleaned <- cleaned[cleaned != "" & !is.na(cleaned)]
  return(cleaned)
}

convert_gene_ids <- function(genes) {
  # Try SYMBOL first, then ENSEMBL
  gene_entrez <- tryCatch({
    bitr(genes, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)

  if (!is.null(gene_entrez) && nrow(gene_entrez) > 0) {
    return(gene_entrez)
  }

  gene_entrez <- tryCatch({
    bitr(genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  }, error = function(e) NULL)

  if (!is.null(gene_entrez) && nrow(gene_entrez) > 0) {
    return(gene_entrez)
  }

  return(NULL)
}

run_enrichment <- function(gene_list, output_dir, prefix) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  genes_cleaned <- clean_gene_id(gene_list)
  if (length(genes_cleaned) == 0) {
    message("⚠️  Empty gene list, skip")
    return(NULL)
  }

  gene_entrez <- convert_gene_ids(genes_cleaned)
  if (is.null(gene_entrez) || nrow(gene_entrez) == 0) {
    message("⚠️  Gene ID conversion failed, skip")
    return(NULL)
  }

  message(paste0("Converted ", nrow(gene_entrez), " genes"))

  GO_result <- tryCatch({
    enrichGO(gene_entrez$ENTREZID,
             OrgDb = org.Hs.eg.db,
             keyType = "ENTREZID",
             ont = "ALL",
             readable = TRUE,
             pvalueCutoff = 0.05,
             qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(GO_result) && nrow(as.data.frame(GO_result)) > 0) {
    go_file <- file.path(output_dir, paste0(prefix, "_GO_results.csv"))
    write.csv(as.data.frame(GO_result), file = go_file, row.names = FALSE)
    if (ggplot_available && enrichplot_available) {
      go_plot_file <- file.path(output_dir, paste0(prefix, "_GO_dotplot.png"))
      tryCatch({
        p <- dotplot(GO_result, x = "GeneRatio", split = "ONTOLOGY", showCategory = 8) +
          facet_grid(ONTOLOGY ~ ., scales = "free_y") +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 14, face = "bold"))
        ggsave(go_plot_file, plot = p, width = 12, height = 14, dpi = 300)
      }, error = function(e) NULL)
    }
  } else {
    message("GO: no significant results")
  }

  KEGG_result <- tryCatch({
    enrichKEGG(gene_entrez$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05, qvalueCutoff = 0.2)
  }, error = function(e) NULL)

  if (!is.null(KEGG_result) && nrow(as.data.frame(KEGG_result)) > 0) {
    kegg_file <- file.path(output_dir, paste0(prefix, "_KEGG_results.csv"))
    write.csv(as.data.frame(KEGG_result), file = kegg_file, row.names = FALSE)
    if (ggplot_available && enrichplot_available) {
      kegg_plot_file <- file.path(output_dir, paste0(prefix, "_KEGG_dotplot.png"))
      tryCatch({
        p <- dotplot(KEGG_result, showCategory = 15) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 12, face = "bold"))
        ggsave(kegg_plot_file, plot = p, width = 10, height = 8, dpi = 300)
      }, error = function(e) NULL)
    }
  } else {
    message("KEGG: no significant results")
  }
}

# ===== Main =====
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  data_set <- args[1]
} else {
  data_set <- "E-GEOD-93593"
}

base_dir <- getwd()
input_root <- file.path(base_dir, "data", data_set, "visualizations")
if (!dir.exists(input_root)) stop(paste0("Missing: ", input_root))

method_dirs <- list.dirs(input_root, recursive = FALSE, full.names = TRUE)

for (method_dir in method_dirs) {
  method_name <- basename(method_dir)
  input_dir <- file.path(method_dir, "enrichment_inputs_consensus")
  if (!dir.exists(input_dir)) next

  files <- list.files(input_dir, pattern = "__consensus\\.csv$", full.names = TRUE)
  if (length(files) == 0) next

  summary_rows <- data.frame(
    method = character(),
    group = character(),
    genes = integer(),
    stringsAsFactors = FALSE
  )

  for (file in files) {
    base <- basename(file)
    parts <- str_split(base, "__", simplify = TRUE)
    group <- parts[1]

    genes <- read.csv(file, stringsAsFactors = FALSE)$gene
    genes <- genes[!is.na(genes) & genes != ""]
    if (length(genes) < 3) {
      message(paste0("Skip ", base, " (", length(genes), " genes)"))
      next
    }

    output_dir <- file.path(method_dir, "enrichment_results_consensus", group)
    prefix <- paste0(method_name, "_", group, "_consensus")
    message(paste0("\n=== ", prefix, " ==="))
    run_enrichment(genes, output_dir, prefix)

    summary_rows <- rbind(summary_rows, data.frame(
      method = method_name,
      group = group,
      genes = length(genes),
      stringsAsFactors = FALSE
    ))
  }

  summary_dir <- file.path(method_dir, "enrichment_results_consensus")
  if (!dir.exists(summary_dir)) dir.create(summary_dir, recursive = TRUE)
  summary_path <- file.path(summary_dir, "enrichment_results_consensus_summary.csv")
  write.csv(summary_rows, file = summary_path, row.names = FALSE)
  message(paste0("Summary saved to: ", summary_path))
}
