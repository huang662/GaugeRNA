###########################################################################
# Figure 1: GIS Distribution Across Methods
# Reproduces Figure 1 from the GaugeRNA paper
# Shows the distribution of Gauge Instability Scores across different
# differential expression methods for multiple datasets
###########################################################################

library(GaugeRNA)

# If ggplot2 is available, use it; otherwise use base R
use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

# Reproduce the same pre-computed results
paper_data <- GaugeRNA::gis_paper_results

cat("========================================\n")
cat("  Figure 1: GIS Distribution Across Methods\n")
cat("========================================\n\n")

# Print the data
print(paper_data)

cat("\nKey findings:\n")
cat("- Gamma-projected method consistently shows lowest GIS\n")
cat("- Raw logFC shows highest GIS across all datasets\n")
cat("- scRNA-seq data shows higher overall GIS than bulk RNA-seq\n")

# Create plot
if (use_ggplot) {
  library(ggplot2)
  p <- ggplot(paper_data, aes(x = Dataset, y = GIS, fill = Method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
    geom_errorbar(aes(ymin = GIS - sd_stability, ymax = GIS + sd_stability),
                  position = position_dodge(width = 0.8), width = 0.25) +
    labs(title = "Figure 1: GIS Distribution Across Methods",
         x = "Dataset", y = "Gauge Instability Score (GIS)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave("figure1_gis_distribution.pdf", plot = p, width = 8, height = 5)
  cat("\nFigure saved to figure1_gis_distribution.pdf\n")
} else {
  # Base R plot
  pdf("figure1_gis_distribution.pdf", width = 8, height = 5)
  datasets <- unique(paper_data$Dataset)
  methods <- unique(paper_data$Method)
  colors <- c("steelblue", "darkgreen", "darkred")
  
  plot_data <- matrix(paper_data$GIS, nrow = length(methods), ncol = length(datasets),
                      byrow = FALSE)
  barplot(plot_data, beside = TRUE, names.arg = datasets,
          col = colors, main = "Figure 1: GIS Distribution Across Methods",
          ylab = "Gauge Instability Score (GIS)", ylim = c(0, 0.6))
  legend("topright", legend = methods, fill = colors)
  dev.off()
  cat("\nFigure saved to figure1_gis_distribution.pdf\n")
}

cat("\nDone.\n")