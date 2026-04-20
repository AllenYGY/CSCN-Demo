# Regenerate GO/KEGG dotplots from existing CSVs (no enrichment rerun)

required_packages <- c("ggplot2")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}
library(ggplot2)

args <- commandArgs(trailingOnly = TRUE)
data_set <- if (length(args) > 0) args[1] else "E-GEOD-93593"

base_dir <- getwd()
input_root <- file.path(base_dir, "data", data_set, "visualizations")
if (!dir.exists(input_root)) {
  stop(paste0("Missing: ", input_root))
}

parse_ratio <- function(x) {
  if (is.na(x) || x == "") return(NA_real_)
  parts <- strsplit(as.character(x), "/", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(NA_real_)
  num <- suppressWarnings(as.numeric(parts[1]))
  den <- suppressWarnings(as.numeric(parts[2]))
  if (is.na(num) || is.na(den) || den == 0) return(NA_real_)
  return(num / den)
}

wrap_label <- function(x, width = 24) {
  if (is.na(x) || x == "") return(x)
  paste(strwrap(as.character(x), width = width), collapse = "\n")
}

make_dotplot <- function(df, out_path, title, kind) {
  if (nrow(df) == 0) return()

  p_col <- if ("p.adjust" %in% names(df)) "p.adjust" else if ("pvalue" %in% names(df)) "pvalue" else NA
  if (is.na(p_col)) return()

  df$pval <- suppressWarnings(as.numeric(df[[p_col]]))
  df$Description <- as.character(df$Description)
  df <- df[!is.na(df$pval) & df$Description != "", ]
  if (nrow(df) == 0) return()

  df$GeneRatio_num <- vapply(df$GeneRatio, parse_ratio, numeric(1))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  df$logp <- -log10(df$pval)

  # select top terms by p-value
  max_terms <- if (kind == "GO") 8 else 15
  df <- df[order(df$pval), ]
  df <- head(df, max_terms)

  df$Description <- vapply(df$Description, wrap_label, character(1), width = 24)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))

  p <- ggplot(df, aes(x = GeneRatio_num, y = Description)) +
    geom_point(aes(size = Count, color = logp)) +
    scale_color_gradient(low = "#A0CBE2", high = "#B22222") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 12),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.title.position = "panel",
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      plot.margin = margin(12, 30, 12, 24)
    ) +
    labs(
      title = title,
      x = "GeneRatio",
      y = NULL,
      color = "-log10(p)"
    )

  if ("ONTOLOGY" %in% names(df)) {
    p <- p + facet_grid(ONTOLOGY ~ ., scales = "free_y")
  }

  height <- max(7, 0.45 * nrow(df) + 5)
  ggsave(out_path, plot = p, width = 12, height = height, dpi = 300)
}

csv_files <- list.files(
  input_root,
  pattern = "(_GO_results|_KEGG_results)\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

for (file in csv_files) {
  df <- tryCatch(read.csv(file, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df)) next

  kind <- if (grepl("_GO_results\\.csv$", file)) "GO" else "KEGG"
  out_path <- sub("_GO_results\\.csv$", "_GO_dotplot.png", file)
  out_path <- sub("_KEGG_results\\.csv$", "_KEGG_dotplot.png", out_path)

  title <- basename(out_path)
  title <- sub("_GO_dotplot\\.png$", "", title)
  title <- sub("_KEGG_dotplot\\.png$", "", title)
  title <- sub("^(kTotal|kWithin|Module Correlation)_", "", title)
  title <- gsub("_", " ", title)
  has_sender <- grepl("\\bsender\\b", title, ignore.case = TRUE)
  has_receiver <- grepl("\\breceiver\\b", title, ignore.case = TRUE)
  title <- gsub("\\bconsensus\\b", "", title, ignore.case = TRUE)
  title <- gsub("\\bsender\\b|\\breceiver\\b", "", title, ignore.case = TRUE)
  title <- gsub("\\s+", " ", title)
  title <- trimws(title)
  role <- ""
  if (has_sender) role <- "sender"
  if (has_receiver) role <- "receiver"
  if (role != "") {
    title <- paste0(title, " (", role, ", ", kind, ")")
  } else {
    title <- paste0(title, " (", kind, ")")
  }

  make_dotplot(df, out_path, title, kind)
}

message("Dotplots regenerated from CSVs.")
