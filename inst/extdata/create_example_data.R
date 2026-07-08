# Create example dataset: simulated_counts_matrix
# This is a small simulated RNA-seq count matrix for package examples

set.seed(1234)
n_genes <- 100
n_samples <- 6

# Base expression
simulated_counts_matrix <- matrix(
  rnbinom(n_genes * n_samples, mu = 100, size = 10),
  nrow = n_genes, ncol = n_samples
)

# Add some differential expression
simulated_counts_matrix[1:10, 4:6] <- simulated_counts_matrix[1:10, 4:6] * 3
simulated_counts_matrix[11:20, 4:6] <- simulated_counts_matrix[11:20, 4:6] * 2

# Name rows and columns
rownames(simulated_counts_matrix) <- paste0("Gene", seq_len(n_genes))
colnames(simulated_counts_matrix) <- paste0("Sample", seq_len(n_samples))

# Save as R data file
save(simulated_counts_matrix, 
     file = "inst/extdata/simulated_counts_matrix.rda")

cat("Example dataset created: simulated_counts_matrix.rda\n")
cat("  Dimensions:", n_genes, "x", n_samples, "\n")
cat("  Genes 1-10: 3x upregulation in samples 4-6\n")
cat("  Genes 11-20: 2x upregulation in samples 4-6\n")