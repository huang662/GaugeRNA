###########################################################################
# Figure 4: Gauge Sensitivity vs Expression Level
# Reproduces Figure 4 from the GaugeRNA paper
# Shows the relationship between gene expression level and GIS
###########################################################################

library(GaugeRNA)

cat("========================================\n")
cat("  Figure 4: Gauge Sensitivity vs Expression Level\n")
cat("========================================\n\n")

# Simulate data with varying expression levels
set.seed(42)
n_genes <- 1000
n_samples <- 10

# Generate genes with varying mean expression
mean_expr <- 10^seq(0, 4, length.out = n_genes)
expr_mat <- matrix(0, nrow = n_genes, ncol = n_samples)
for (g in seq_len(n_genes)) {
  expr_mat[g, ] <- rnbinom(n_samples, mu = mean_expr[g], size = 10)
}

condition <- c(rep(0, 5), rep(1, 5))

# Generate a full audit report
report <- audit_report(expr_mat, condition, n_perturbations = 15)

# Extract per-gene GIS and expression
per_gene_gis <- report$classification$GIS
mean_expr_log <- log10(rowMeans(expr_mat) + 1)

cat("Correlation between expression level and GIS:\n")
cor_test <- cor.test(mean_expr_log, per_gene_gis, method = "spearman")
cat("  Spearman rho =", round(cor_test$estimate, 4), "\n")
cat("  p-value =", format(cor_test$p.value, digits = 4), "\n")

# Create figure
use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

if (use_ggplot) {
  library(ggplot2)
  df <- data.frame(
    Expression = mean_expr_log,
    GIS = per_gene_gis,
    Classification = report$classification$classification
  )
  
  p <- ggplot(df, aes(x = Expression, y = GIS, color = Classification)) +
    geom_point(alpha = 0.5, size = 0.8) +
    geom_smooth(method = "loess", se = TRUE, color = "black") +
    labs(title = "Figure 4: Gauge Sensitivity vs Expression Level",
         x = "log10(Mean Expression + 1)", y = "Per-Gene GIS") +
    theme_minimal()
  
  ggsave("figure4_gis_vs_expression.pdf", plot = p, width = 8, height = 5)
  cat("\nFigure saved to figure4_gis_vs_expression.pdf\n")
} else {
  pdf("figure4_gis_vs_expression.pdf", width = 8, height = 5)
  plot(mean_expr_log, per_gene_gis,
       col = ifelse(report$classification$classification == "gauge-sensitive",
                    rgb(0.8, 0, 0, 0.3), rgb(0, 0.5, 0, 0.3)),
       pch = 16, cex = 0.5,
       main = "Figure 4: Gauge Sensitivity vs Expression Level",
       xlab = "log10(Mean Expression + 1)", ylab = "Per-Gene GIS")
  legend("topright", legend = c("Gauge-sensitive", "Gauge-stable"),
         col = c(rgb(0.8, 0, 0, 0.8), rgb(0, 0.5, 0, 0.8)), pch = 16)
  # Add loess smooth
  lo <- loess(per_gene_gis ~ mean_expr_log)
  idx <- order(mean_expr_log)
  lines(mean_expr_log[idx], predict(lo)[idx], col = "black", lwd = 2)
  dev.off()
  cat("\nFigure saved to figure4_gis_vs_expression.pdf\n")
}

cat("\nDone.\n")