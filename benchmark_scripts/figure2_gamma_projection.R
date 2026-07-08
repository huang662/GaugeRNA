###########################################################################
# Figure 2: Gamma Projection Visualization
# Reproduces Figure 2 from the GaugeRNA paper
# Demonstrates the effect of gamma projection on log-fold changes
###########################################################################

library(GaugeRNA)

cat("========================================\n")
cat("  Figure 2: Gamma Projection Visualization\n")
cat("========================================\n\n")

# Simulate data with a known gauge shift
set.seed(42)
n_genes <- 2000
logfc <- rnorm(n_genes, mean = 0.5, sd = 1.0)

# Add some true differential expression
de_idx <- sample(seq_len(n_genes), 200)
logfc[de_idx] <- logfc[de_idx] + 2

# Apply gamma projection
gamma_logfc <- gamma_projection(logfc)

cat("Raw logFC summary:\n")
print(summary(logfc))
cat("\nGamma-projected logFC summary:\n")
print(summary(gamma_logfc))
cat("\nMedian of gamma-projected logFC:", median(gamma_logfc), "\n")

# Create figure
use_ggplot <- requireNamespace("ggplot2", quietly = TRUE)

if (use_ggplot) {
  library(ggplot2)
  df <- data.frame(
    logFC = c(logfc, gamma_logfc),
    Type = rep(c("Raw logFC", "Gamma-projected"), each = n_genes)
  )
  
  p <- ggplot(df, aes(x = logFC, fill = Type, color = Type)) +
    geom_density(alpha = 0.4) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Figure 2: Effect of Gamma Projection on LogFC Distribution",
         x = "Log Fold Change", y = "Density") +
    theme_minimal()
  
  ggsave("figure2_gamma_projection.pdf", plot = p, width = 8, height = 5)
  cat("\nFigure saved to figure2_gamma_projection.pdf\n")
} else {
  pdf("figure2_gamma_projection.pdf", width = 8, height = 5)
  plot(density(logfc), col = "blue", lwd = 2,
       main = "Figure 2: Effect of Gamma Projection on LogFC",
       xlab = "Log Fold Change", ylim = c(0, max(density(logfc)$y, density(gamma_logfc)$y)))
  lines(density(gamma_logfc), col = "red", lwd = 2)
  abline(v = 0, lty = 2, col = "gray50")
  legend("topright", legend = c("Raw logFC", "Gamma-projected"), 
         col = c("blue", "red"), lwd = 2)
  dev.off()
  cat("\nFigure saved to figure2_gamma_projection.pdf\n")
}

cat("\nDone.\n")