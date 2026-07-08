###########################################################################
# Figure 3: Gauge Decomposition Heatmap
# Reproduces Figure 3 from the GaugeRNA paper
# Shows the three components of gauge decomposition
###########################################################################

library(GaugeRNA)

cat("========================================\n")
cat("  Figure 3: Gauge Decomposition\n")
cat("========================================\n\n")

# Simulate structured data
set.seed(42)
n_genes <- 300
n_samples <- 12

# Create structured expression data
gene_means <- rnorm(n_genes, mean = 8, sd = 2)
sample_effects <- c(rep(-0.5, 4), rep(0.2, 4), rep(0.3, 4))
# Add some gene-specific patterns
expr_mat <- matrix(0, nrow = n_genes, ncol = n_samples)
for (g in seq_len(n_genes)) {
  for (s in seq_len(n_samples)) {
    expr_mat[g, s] <- rnbinom(1, mu = 2^(gene_means[g] + sample_effects[s] + 
                                         rnorm(1, 0, 0.3)), size = 10)
  }
}

rownames(expr_mat) <- paste0("Gene", seq_len(n_genes))
colnames(expr_mat) <- paste0("Sample", seq_len(n_samples))

# Perform decomposition
decomp <- gauge_decomposition(expr_mat)

cat("Decomposition summary:\n")
print(decomp)
cat("\nVariance components:\n")
cat("  Gene gauge variance:", round(var(decomp$mu), 3), "\n")
cat("  Sample gauge variance:", round(var(decomp$a), 3), "\n")
cat("  Curvature variance:", round(var(as.vector(decomp$epsilon)), 3), "\n")

# Create figure
pdf("figure3_gauge_decomposition.pdf", width = 12, height = 5)

par(mfrow = c(1, 3))

# Gene gauge
barplot(sort(decomp$mu)[1:50], main = "Gene Gauge (top 50 genes)",
        xlab = "Gene rank", ylab = "log2 expression", col = "steelblue",
        border = NA)

# Sample gauge
barplot(decomp$a, main = "Sample Gauge",
        names.arg = colnames(expr_mat), col = "darkgreen",
        border = NA, ylab = "log2 shift")

# Curvature heatmap
image(t(decomp$epsilon[1:50, ]), main = "Curvature (first 50 genes)",
      xlab = "Samples", ylab = "Genes", col = colorRampPalette(c("blue", "white", "red"))(100))

par(mfrow = c(1, 1))
dev.off()

cat("\nFigure saved to figure3_gauge_decomposition.pdf\n")
cat("Done.\n")