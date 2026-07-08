#' Generate a Complete GaugeRNA Audit Report
#'
#' Generates a comprehensive audit report for RNA-seq differential expression
#' analysis. The report includes Gauge Instability Scores (GIS), gene
#' classification, gauge decomposition, and actionable recommendations.
#'
#' This is the main entry point for the GaugeRNA framework. It performs all
#' analyses in sequence and returns a structured report object that can be
#' printed, summarized, or used to generate tables and figures.
#'
#' @param expr_mat A numeric matrix of expression counts (genes x samples).
#'   Rows are genes, columns are samples.
#' @param condition A binary vector of length \code{ncol(expr_mat)} indicating
#'   the condition for each sample.
#' @param n_perturbations Integer. Number of random gauge perturbations.
#'   Default is 20.
#' @param sigma Numeric. Noise standard deviation for perturbations. Default is 1.0.
#' @param methods Character vector. Methods to evaluate. Default is
#'   \code{c("raw", "gamma", "deseq2")}.
#' @param gis_threshold Classification threshold. Passed to \code{classify_genes()}.
#'   Default is \code{"auto"}.
#'
#' @return A list of class \code{gaugerna_report} with components:
#'   \item{summary}{A data.frame with overall summary statistics.}
#'   \item{gis_results}{The GIS results from \code{compute_gis()}.}
#'   \item{classification}{The gene classification from \code{classify_genes()}.}
#'   \item{decomposition}{The gauge decomposition from \code{gauge_decomposition()}.}
#'   \item{recommendations}{Character vector of actionable recommendations.}
#'   \item{parameters}{List of parameters used.}
#'   \item{timestamp}{When the report was generated.}
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' expr_mat <- matrix(rnbinom(500 * 10, mu = 100, size = 10), nrow = 500, ncol = 10)
#' condition <- c(rep(0, 5), rep(1, 5))
#' report <- audit_report(expr_mat, condition, n_perturbations = 15)
#' print(report)
#' }
#'
#' @export
audit_report <- function(expr_mat, condition, n_perturbations = 20, sigma = 1.0,
                          methods = c("raw", "gamma", "deseq2"),
                          gis_threshold = "auto") {
  # Input validation
  if (!is.matrix(expr_mat)) {
    stop("expr_mat must be a numeric matrix (genes x samples)")
  }
  if (length(condition) != ncol(expr_mat)) {
    stop("length(condition) must equal ncol(expr_mat)")
  }
  if (!all(condition %in% c(0, 1))) {
    stop("condition must be a binary vector with values 0 and 1")
  }
  methods <- match.arg(methods, several.ok = TRUE)

  # Step 1: Compute GIS
  gis_results <- compute_gis(
    expr_mat = expr_mat,
    condition = condition,
    n_perturbations = n_perturbations,
    sigma = sigma,
    methods = methods
  )

  # Step 2: Compute per-gene GIS
  # For per-gene GIS, use the gamma method by default
  per_gene_gis <- .compute_per_gene_gis(expr_mat, condition, sigma, n_perturbations)

  # Step 3: Classify genes
  classification <- classify_genes(
    gis_vector = per_gene_gis,
    threshold = gis_threshold
  )

  # Step 4: Gauge decomposition
  decomposition <- gauge_decomposition(expr_mat)

  # Step 5: Generate recommendations
  recommendations <- .generate_recommendations(gis_results, classification, decomposition)

  # Step 6: Summary
  summary <- .generate_summary(gis_results, classification, decomposition)

  report <- list(
    summary = summary,
    gis_results = gis_results,
    classification = classification,
    decomposition = decomposition,
    recommendations = recommendations,
    parameters = list(
      n_genes = nrow(expr_mat),
      n_samples = ncol(expr_mat),
      n_condition_0 = sum(condition == 0),
      n_condition_1 = sum(condition == 1),
      n_perturbations = n_perturbations,
      sigma = sigma,
      methods = methods,
      gis_threshold = gis_threshold
    ),
    timestamp = Sys.time()
  )

  class(report) <- "gaugerna_report"
  return(report)
}

#' Compute per-gene GIS values
#' @keywords internal
.compute_per_gene_gis <- function(expr_mat, condition, sigma, n_perturbations) {
  n_genes <- nrow(expr_mat)
  n_samples <- ncol(expr_mat)

  baseline_logfc <- .compute_baseline_logfc(expr_mat, condition)

  # Per-gene stability across perturbations
  all_logfc <- matrix(NA, nrow = n_genes, ncol = n_perturbations)

  for (p in seq_len(n_perturbations)) {
    perturbed <- .apply_gauge_perturbation(expr_mat, sigma)
    perturbed_logfc <- .compute_baseline_logfc(perturbed, condition)
    all_logfc[, p] <- perturbed_logfc
  }

  # For each gene, compute stability as the correlation between baseline
  # and a reference (use mean of perturbations as reference)
  ref_logfc <- rowMeans(all_logfc, na.rm = TRUE)

  per_gene_gis <- numeric(n_genes)
  for (g in seq_len(n_genes)) {
    if (is.na(baseline_logfc[g])) {
      per_gene_gis[g] <- 1.0
      next
    }
    gene_vals <- all_logfc[g, ]
    valid_idx <- which(!is.na(gene_vals) & is.finite(gene_vals))
    if (length(valid_idx) < 3) {
      per_gene_gis[g] <- 1.0
      next
    }
    # GIS = variance of logFC across perturbations (normalized)
    per_gene_gis[g] <- stats::var(gene_vals[valid_idx]) /
                        (stats::var(gene_vals[valid_idx]) + 1)
  }

  gene_names <- rownames(expr_mat)
  if (is.null(gene_names)) {
    gene_names <- paste0("Gene", seq_len(n_genes))
  }
  names(per_gene_gis) <- gene_names

  return(per_gene_gis)
}

#' Generate recommendations based on audit results
#' @keywords internal
.generate_recommendations <- function(gis_results, classification, decomposition) {
  recs <- character()

  # Overall GIS assessment
  best_method <- gis_results$method[which.min(gis_results$GIS)]
  worst_method <- gis_results$method[which.max(gis_results$GIS)]
  min_gis <- min(gis_results$GIS)
  max_gis <- max(gis_results$GIS)

  recs <- c(recs, sprintf(
    "Overall GIS range: %.3f (%s) to %.3f (%s). Lower GIS indicates more gauge-stable results.",
    min_gis, best_method, max_gis, worst_method
  ))

  if (min_gis < 0.1) {
    recs <- c(recs, sprintf(
      "RECOMMENDATION: Method '%s' shows excellent gauge stability (GIS = %.3f). Use this method for primary analysis.",
      best_method, min_gis
    ))
  } else if (min_gis < 0.3) {
    recs <- c(recs, sprintf(
      "RECOMMENDATION: Method '%s' shows moderate gauge stability (GIS = %.3f). Results should be interpreted with some caution.",
      best_method, min_gis
    ))
  } else {
    recs <- c(recs, "WARNING: No method shows strong gauge stability (all GIS > 0.3). Consider alternative normalization strategies.")
  }

  # Classification summary
  n_sensitive <- sum(classification$classification == "gauge-sensitive")
  n_stable <- sum(classification$classification == "gauge-stable")
  n_total <- nrow(classification)

  recs <- c(recs, sprintf(
    "Gene classification: %d/%d (%.1f%%) genes are gauge-sensitive, %d/%d (%.1f%%) are gauge-stable.",
    n_sensitive, n_total, 100 * n_sensitive / n_total,
    n_stable, n_total, 100 * n_stable / n_total
  ))

  if (n_sensitive / n_total > 0.5) {
    recs <- c(recs, "WARNING: More than half of genes are gauge-sensitive. The dataset may have high technical variability.")
  }

  # Decomposition assessment
  curvature_var <- stats::var(as.vector(decomposition$epsilon))
  sample_gauge_var <- stats::var(decomposition$a)
  recs <- c(recs, sprintf(
    "Gauge decomposition: sample gauge variance = %.3f, curvature variance = %.3f.",
    sample_gauge_var, curvature_var
  ))

  if (sample_gauge_var > curvature_var) {
    recs <- c(recs, "NOTE: Sample gauge dominates curvature. Strong batch effects may be present.")
  }

  # General best practices
  recs <- c(recs, "BEST PRACTICE: Report GIS values alongside differential expression results.")
  recs <- c(recs, "BEST PRACTICE: Flag gauge-sensitive genes in supplementary tables.")
  recs <- c(recs, "BEST PRACTICE: Use gamma-projected log-fold changes for downstream analysis.")

  return(recs)
}

#' Generate summary statistics
#' @keywords internal
.generate_summary <- function(gis_results, classification, decomposition) {
  n_sensitive <- sum(classification$classification == "gauge-sensitive")
  n_stable <- sum(classification$classification == "gauge-stable")
  n_total <- nrow(classification)

  summary_df <- data.frame(
    metric = c(
      "n_genes",
      "n_samples",
      "best_method",
      "best_GIS",
      "worst_method",
      "worst_GIS",
      "n_gauge_sensitive",
      "n_gauge_stable",
      "pct_sensitive",
      "gene_gauge_mean",
      "sample_gauge_sd",
      "curvature_sd"
    ),
    value = c(
      nrow(decomposition$epsilon),
      ncol(decomposition$epsilon),
      gis_results$method[which.min(gis_results$GIS)],
      sprintf("%.4f", min(gis_results$GIS)),
      gis_results$method[which.max(gis_results$GIS)],
      sprintf("%.4f", max(gis_results$GIS)),
      n_sensitive,
      n_stable,
      sprintf("%.1f%%", 100 * n_sensitive / n_total),
      sprintf("%.3f", mean(decomposition$mu)),
      sprintf("%.3f", stats::sd(decomposition$a)),
      sprintf("%.3f", stats::sd(as.vector(decomposition$epsilon)))
    ),
    stringsAsFactors = FALSE
  )

  return(summary_df)
}

#' @export
print.gaugerna_report <- function(x, ...) {
  cat("========================================\n")
  cat("  GaugeRNA Audit Report\n")
  cat("========================================\n")
  cat("Generated: ", format(x$timestamp), "\n\n")

  cat("--- Parameters ---\n")
  cat("  Genes:", x$parameters$n_genes, "\n")
  cat("  Samples:", x$parameters$n_samples, "\n")
  cat("  Condition 0:", x$parameters$n_condition_0, "\n")
  cat("  Condition 1:", x$parameters$n_condition_1, "\n")
  cat("  Perturbations:", x$parameters$n_perturbations, "\n")
  cat("  Sigma:", x$parameters$sigma, "\n")
  cat("  Methods:", paste(x$parameters$methods, collapse = ", "), "\n\n")

  cat("--- GIS Results ---\n")
  print(x$gis_results, row.names = FALSE)
  cat("\n")

  cat("--- Classification Summary ---\n")
  cat("  Gauge-sensitive:", sum(x$classification$classification == "gauge-sensitive"), "\n")
  cat("  Gauge-stable:", sum(x$classification$classification == "gauge-stable"), "\n\n")

  cat("--- Recommendations ---\n")
  for (i in seq_along(x$recommendations)) {
    cat("  ", i, ". ", x$recommendations[i], "\n", sep = "")
  }
  cat("\n========================================\n")

  invisible(x)
}

#' @export
summary.gaugerna_report <- function(object, ...) {
  cat("--- Audit Summary ---\n")
  cat("  Best method:", as.character(object$summary$value[object$summary$metric == "best_method"]), "\n")
  cat("  Best GIS:", as.character(object$summary$value[object$summary$metric == "best_GIS"]), "\n")
  cat("  Worst method:", as.character(object$summary$value[object$summary$metric == "worst_method"]), "\n")
  cat("  Worst GIS:", as.character(object$summary$value[object$summary$metric == "worst_GIS"]), "\n")
  cat("  Gauge-sensitive genes:", as.character(object$summary$value[object$summary$metric == "n_gauge_sensitive"]), "\n")
  cat("  Gauge-stable genes:", as.character(object$summary$value[object$summary$metric == "n_gauge_stable"]), "\n")
  cat("\n")
  invisible(object$summary)
}