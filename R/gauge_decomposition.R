#' Gauge Decomposition of Log-Count Matrix
#'
#' Decomposes a log-transformed count matrix into three components:
#' gene gauge (mu), sample gauge (a), and curvature (epsilon).
#'
#' The gauge decomposition model is:
#' \deqn{Y_{gs} = \mu_g + a_s + \epsilon_{gs}}
#' where \eqn{Y_{gs} = \log(K_{gs} + 1)} is the log-count for gene g and
#' sample s. The gene gauge \eqn{\mu_g} captures the average expression
#' level of each gene. The sample gauge \eqn{a_s} captures the global
#' shift common to all genes in a given sample. The curvature
#' \eqn{\epsilon_{gs}} is the residual that captures gene-sample-specific
#' variation not explained by the gauge components.
#'
#' @param expr_mat A numeric matrix of expression counts (genes x samples).
#'   Rows are genes, columns are samples.
#' @param pseudo_count Numeric. Pseudocount added before log transformation.
#'   Default is 1.
#'
#' @return A list with components:
#'   \item{mu}{Numeric vector of gene gauges (length = n_genes).}
#'   \item{a}{Numeric vector of sample gauges (length = n_samples).}
#'   \item{epsilon}{Numeric matrix of curvature residuals (genes x samples).}
#'   \item{log_counts}{The log-transformed count matrix.}
#'   \item{gene_names}{Character vector of gene names.}
#'   \item{sample_names}{Character vector of sample names.}
#'
#' @examples
#' set.seed(42)
#' expr_mat <- matrix(rnbinom(100 * 6, mu = 100, size = 10), nrow = 100, ncol = 6)
#' decomp <- gauge_decomposition(expr_mat)
#' str(decomp)
#' # The curvature matrix should have mean near 0
#' mean(decomp$epsilon)
#'
#' @export
gauge_decomposition <- function(expr_mat, pseudo_count = 1) {
  if (!is.matrix(expr_mat)) {
    stop("expr_mat must be a numeric matrix (genes x samples)")
  }
  if (pseudo_count < 0) {
    stop("pseudo_count must be non-negative")
  }

  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)

  # Log transform
  log_counts <- log2(expr_mat + pseudo_count)

  # Gene gauge: mean expression per gene
  mu <- rowMeans(log_counts)

  # Sample gauge: mean expression per sample (centered)
  sample_means <- colMeans(log_counts)
  a <- sample_means - mean(sample_means)

  # Curvature: residual after removing gene and sample gauges
  mu_mat <- matrix(mu, nrow = n_genes, ncol = n_samples)
  a_mat <- matrix(a, nrow = n_genes, ncol = n_samples, byrow = TRUE)
  epsilon <- log_counts - mu_mat - a_mat

  # Get names
  gene_names <- rownames(expr_mat)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene", seq_len(n_genes))
  }

  sample_names <- colnames(expr_mat)
  if (is.null(sample_names)) {
    sample_names <- paste0("Sample", seq_len(n_samples))
  }

  # Set names
  names(mu) <- gene_names
  names(a) <- sample_names
  rownames(epsilon) <- gene_names
  colnames(epsilon) <- sample_names

  result <- list(
    mu = mu,
    a = a,
    epsilon = epsilon,
    log_counts = log_counts,
    gene_names = gene_names,
    sample_names = sample_names
  )

  class(result) <- "gaugerna_decomposition"
  return(result)
}

#' @export
print.gaugerna_decomposition <- function(x, ...) {
  cat("Gauge Decomposition\n")
  cat("  Genes: ", length(x$mu), "\n")
  cat("  Samples: ", length(x$a), "\n")
  cat("  Gene gauge (mu) range: [",
      sprintf("%.3f", min(x$mu)), ", ",
      sprintf("%.3f", max(x$mu)), "]\n", sep = "")
  cat("  Sample gauge (a) range: [",
      sprintf("%.3f", min(x$a)), ", ",
      sprintf("%.3f", max(x$a)), "]\n", sep = "")
  cat("  Curvature (epsilon) range: [",
      sprintf("%.3f", min(x$epsilon)), ", ",
      sprintf("%.3f", max(x$epsilon)), "]\n", sep = "")
  invisible(x)
}