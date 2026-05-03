options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(clusterProfiler)
  library(enrichplot)
})

parse_ratio <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NA_real_)
  num <- suppressWarnings(as.numeric(parts[1]))
  den <- suppressWarnings(as.numeric(parts[2]))
  if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
  num / den
}

wrap_label <- function(x, width = 30) {
  if (is.na(x) || x == "") return(x)
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

make_go_dotplot <- function(csv_path, out_path, title_text, top_n = 8) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  if (!"Description" %in% names(df)) stop("Missing Description column")
  if (!"p.adjust" %in% names(df) && !"pvalue" %in% names(df)) stop("Missing p.adjust/pvalue column")

  p_col <- if ("p.adjust" %in% names(df)) "p.adjust" else "pvalue"
  df$pval <- suppressWarnings(as.numeric(df[[p_col]]))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  df$GeneRatio_num <- vapply(df$GeneRatio, parse_ratio, numeric(1))
  df <- df[!is.na(df$pval) & !is.na(df$Count) & !is.na(df$GeneRatio_num) & df$Description != "", , drop = FALSE]
  if (nrow(df) == 0) stop("No GO BP rows available")

  # reconstruct a lightweight enrichResult-like object by using the CSV as the plotting backend
  # and reapply the original plotting grammar from GO_KEGG.R
  if ("ONTOLOGY" %in% names(df)) {
    split_df <- split(df, df$ONTOLOGY)
    split_df <- lapply(split_df, function(chunk) {
      chunk <- chunk[order(chunk$pval), , drop = FALSE]
      head(chunk, top_n)
    })
    plot_df <- do.call(rbind, split_df)
  } else {
    plot_df <- df[order(df$pval), , drop = FALSE]
    plot_df <- head(plot_df, top_n)
  }

  plot_df$Description <- vapply(plot_df$Description, wrap_label, character(1), width = 30)
  plot_df$Description <- factor(plot_df$Description, levels = rev(unique(plot_df$Description)))
  plot_df$logp <- -log10(pmax(plot_df$pval, 1e-300))

  p <- ggplot(plot_df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = logp)) +
    scale_size(range = c(4, 10)) +
    scale_color_gradient(low = "blue", high = "red") +
    theme_minimal() +
    theme(
      plot.margin = margin(20, 40, 70, 24, "pt"),
      plot.title = element_text(size = 30, hjust = 0.5),
      axis.text.y = element_text(size = 20, face = "bold", angle = 15, hjust = 1),
      axis.text.x = element_text(size = 8),
      strip.text = element_text(size = 9),
      legend.position = "right",
      legend.text = element_text(size = 8),
      legend.title = element_text(size = 9),
      legend.key.size = grid::unit(0.8, "cm")
    ) +
    labs(
      title = title_text,
      x = "GeneRatio",
      y = NULL,
      color = "-log10(p.adjust)",
      size = "Count"
    ) +
    coord_cartesian(clip = "off")

  if ("ONTOLOGY" %in% names(plot_df)) {
    p <- p + facet_grid(ONTOLOGY ~ ., space = "free_y", scales = "free_y")
  }

  plot_height <- max(12, 0.9 * nrow(plot_df) + 4.5)
  ggsave(out_path, plot = p, width = 16, height = plot_height, dpi = 300, limitsize = FALSE)
}

jobs <- list(
  list(
    csv = "/Users/allenygy/Research/CSCN/data/SCP259/enrichment_results/Biomarkers/SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_Biomarkers_GO_results.csv",
    out = "/Users/allenygy/Research/CSCN/data/SCP259/enrichment_results/Biomarkers/SCP259_inflamed_vs_healthy_crypt_prolif_epi_Gene_Biomarkers_GO_dotplot.png",
    title = "SCP259 biomarker GO enrichment"
  ),
  list(
    csv = "/Users/allenygy/Research/CSCN/data/GSE159115/enrichment_results/Biomarkers/GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_Biomarkers_GO_results.csv",
    out = "/Users/allenygy/Research/CSCN/data/GSE159115/enrichment_results/Biomarkers/GSE159115_ccrcc_tumor_vs_ptb_ptc_normal_Gene_Biomarkers_GO_dotplot.png",
    title = "GSE159115 biomarker GO enrichment"
  ),
  list(
    csv = "/Users/allenygy/Research/CSCN/data/GSE138852/enrichment_results/Biomarkers/GSE138852_Gene_Biomarkers_GO_results.csv",
    out = "/Users/allenygy/Research/CSCN/data/GSE138852/enrichment_results/Biomarkers/GSE138852_Gene_Biomarkers_GO_dotplot.png",
    title = "GSE138852 biomarker GO enrichment"
  )
)

for (job in jobs) {
  message("Regenerating: ", job$out)
  make_go_dotplot(job$csv, job$out, job$title)
}

message("Done.")
