#' Compute Gauge Instability Score (GIS)
#'
#' Computes the Gauge Instability Score for each differential expression method
#' by perturbing the expression matrix with random gauge transformations and
#' measuring the stability of log-fold changes.
#'
#' The GIS quantifies how sensitive differential expression results are to
#' arbitrary normalization choices (gauge choices). A low GIS indicates
#' gauge-stable results; a high GIS indicates gauge-sensitive results that
#' may not be robust.
#'
#' @param expr_mat A numeric matrix of expression values (genes x samples).
#'   Rows are genes, columns are samples.
#' @param condition A binary vector of length \code{ncol(expr_mat)} indicating
#'   the condition for each sample (e.g., \code{c(0,0,0,1,1,1)}).
#' @param n_perturbations Integer. Number of random gauge perturbations to apply.
#'   Default is 20.
#' @param sigma Numeric. Standard deviation of the random noise added during
#'   perturbation. Default is 1.0.
#' @param methods Character vector. Differential expression methods to evaluate.
#'   Default is \code{c("raw", "gamma", "deseq2")}. Supported methods:
#'   \code{"raw"} (log-fold change), \code{"gamma"} (gamma-projected logFC),
#'   \code{"deseq2"} (DESeq2-style variance-stabilized).
#'
#' @return A data.frame with columns:
#'   \item{method}{The differential expression method name.}
#'   \item{GIS}{The Gauge Instability Score for that method.}
#'   \item{mean_stability}{Mean stability across perturbations.}
#'   \item{sd_stability}{Standard deviation of stability across perturbations.}
#'
#' @examples
#' # Simulate expression data
#' set.seed(42)
#' expr_mat <- matrix(rnbinom(100 * 6, mu = 100, size = 10), nrow = 100, ncol = 6)
#' condition <- c(0, 0, 0, 1, 1, 1)
#' gis_results <- compute_gis(expr_mat, condition, n_perturbations = 10)
#' print(gis_results)
#'
#' @export
compute_gis <- function(expr_mat, condition, n_perturbations = 20, sigma = 1.0,
                         methods = c("raw", "gamma", "deseq2")) {
  # Input validation
  if (!is.matrix(expr_mat)) {
    stop("expr_mat must be a numeric matrix (genes x samples)")
  }
  if (any(expr_mat < 0, na.rm = TRUE)) {
    stop("expr_mat contains negative values; expected non-negative counts")
  }
  if (length(condition) != ncol(expr_mat)) {
    stop("length(condition) must equal ncol(expr_mat)")
  }
  if (!all(condition %in% c(0, 1))) {
    stop("condition must be a binary vector with values 0 and 1")
  }
  if (sum(condition == 0) < 2 || sum(condition == 1) < 2) {
    stop("Each condition must have at least 2 samples")
  }
  if (n_perturbations < 1) {
    stop("n_perturbations must be at least 1")
  }
  methods <- match.arg(methods, several.ok = TRUE)

  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)

  # Compute baseline log-fold changes
  baseline_logfc <- .compute_baseline_logfc(expr_mat, condition)

  results <- data.frame(
    method = character(0),
    GIS = numeric(0),
    mean_stability = numeric(0),
    sd_stability = numeric(0),
    stringsAsFactors = FALSE
  )

  for (m in methods) {
    stability_scores <- numeric(n_perturbations)

    for (p in seq_len(n_perturbations)) {
      # Apply random gauge perturbation
      perturbed <- .apply_gauge_perturbation(expr_mat, sigma)
      perturbed_logfc <- .compute_method_logfc(perturbed, condition, method = m)

      # Compute stability as correlation between baseline and perturbed logFC
      valid_idx <- which(!is.na(baseline_logfc) & !is.na(perturbed_logfc) &
                         is.finite(baseline_logfc) & is.finite(perturbed_logfc))
      if (length(valid_idx) < 10) {
        stability_scores[p] <- 0
      } else {
        stability_scores[p] <- suppressWarnings(
          cor(baseline_logfc[valid_idx], perturbed_logfc[valid_idx],
              method = "pearson")
        )
        if (is.na(stability_scores[p])) stability_scores[p] <- 0
      }
    }

    # GIS = 1 - mean(stability), higher GIS = more gauge-sensitive
    mean_stab <- mean(stability_scores)
    sd_stab <- sd(stability_scores)
    gis_value <- 1.0 - mean_stab

    results <- rbind(results, data.frame(
      method = m,
      GIS = gis_value,
      mean_stability = mean_stab,
      sd_stability = sd_stab,
      stringsAsFactors = FALSE
    ))
  }

  rownames(results) <- NULL
  class(results) <- c("gaugerna_gis", "data.frame")
  return(results)
}

#' Print method for gaugerna_gis objects
#' @param x A gaugerna_gis object
#' @param ... Additional arguments (ignored)
#' @export
print.gaugerna_gis <- function(x, ...) {
  cat("GaugeRNA GIS Results\n")
  cat("====================\n")
  cat(sprintf("  Methods: %d\n", nrow(x)))
  n_pert <- attr(x, "n_perturbations")
  cat(sprintf("  Perturbations: %d\n\n", ifelse(is.null(n_pert), NA, n_pert)))
  print.data.frame(x, row.names = FALSE, digits = 4)
  invisible(x)
}

# --- Internal helper functions ---

#' Compute baseline log-fold change
#' @keywords internal
.compute_baseline_logfc <- function(expr_mat, condition) {
  pseudocount <- 1
  ctrl_idx <- which(condition == 0)
  case_idx <- which(condition == 1)

  ctrl_mean <- rowMeans(log2(expr_mat[, ctrl_idx, drop = FALSE] + pseudocount))
  case_mean <- rowMeans(log2(expr_mat[, case_idx, drop = FALSE] + pseudocount))

  return(case_mean - ctrl_mean)
}

#' Apply random gauge perturbation to expression matrix
#' @keywords internal
.apply_gauge_perturbation <- function(expr_mat, sigma) {
  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)

  # Generate random gauge factors per gene and per sample
  gene_gauge <- exp(stats::rnorm(n_genes, mean = 0, sd = sigma))
  sample_gauge <- exp(stats::rnorm(n_samples, mean = 0, sd = sigma))

  # Apply multiplicative perturbation
  perturbed <- sweep(expr_mat, 1, gene_gauge, "*")
  perturbed <- sweep(perturbed, 2, sample_gauge, "*")

  # Ensure non-negative
  perturbed[perturbed < 0] <- 0

  return(perturbed)
}

#' Compute method-specific log-fold change
#' @keywords internal
.compute_method_logfc <- function(expr_mat, condition, method) {
  pseudocount <- 1
  ctrl_idx <- which(condition == 0)
  case_idx <- which(condition == 1)

  ctrl_mean <- rowMeans(log2(expr_mat[, ctrl_idx, drop = FALSE] + pseudocount))
  case_mean <- rowMeans(log2(expr_mat[, case_idx, drop = FALSE] + pseudocount))

  raw_logfc <- case_mean - ctrl_mean

  switch(method,
    "raw" = raw_logfc,
    "gamma" = gamma_projection(raw_logfc),
    "deseq2" = .deseq2_approx(expr_mat, condition),
    raw_logfc
  )
}

#' Simple DESeq2-like approximation using variance-stabilizing normalization
#' @keywords internal
.deseq2_approx <- function(expr_mat, condition) {
  pseudocount <- 1
  ctrl_idx <- which(condition == 0)
  case_idx <- which(condition == 1)

  # Compute size factors (geometric mean normalization)
  geo_means <- exp(rowMeans(log(expr_mat + pseudocount)))
  sf <- colSums(expr_mat) / mean(colSums(expr_mat))

  # Normalized counts
  norm_mat <- sweep(expr_mat, 2, sf, "/")

  # Simple variance-stabilized logFC
  ctrl_var <- apply(norm_mat[, ctrl_idx, drop = FALSE], 1, stats::var)
  case_var <- apply(norm_mat[, case_idx, drop = FALSE], 1, stats::var)
  pooled_var <- (ctrl_var + case_var) / 2

  ctrl_mean <- rowMeans(log2(norm_mat[, ctrl_idx, drop = FALSE] + pseudocount))
  case_mean <- rowMeans(log2(norm_mat[, case_idx, drop = FALSE] + pseudocount))

  raw_logfc <- case_mean - ctrl_mean

  # Shrink extreme low-count genes
  shrinkage <- sqrt(pooled_var) / (sqrt(pooled_var) + 1)
  shrinkage[is.na(shrinkage)] <- 0
  return(raw_logfc * shrinkage)
}