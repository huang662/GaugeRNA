#' Classify Genes into Gauge-Sensitive and Gauge-Stable
#'
#' Classifies genes based on their Gauge Instability Scores (GIS) into
#' gauge-sensitive and gauge-stable categories.
#'
#' A gene is classified as gauge-sensitive if its GIS exceeds the specified
#' threshold, and gauge-stable otherwise. This classification helps identify
#' genes whose differential expression results should be interpreted with
#' caution due to gauge dependency.
#'
#' @param gis_vector A numeric vector of per-gene Gauge Instability Scores.
#'   Must be named with gene identifiers.
#' @param threshold Numeric. The GIS threshold for classification. Default is
#'   the 75th percentile of the GIS values (\code{threshold = "auto"}).
#'   If \code{"auto"}, uses the 75th percentile. If \code{"median"}, uses
#'   the median. If a numeric value, uses that value directly.
#'
#' @return A data.frame with columns:
#'   \item{gene}{Gene identifier (from names of \code{gis_vector}).}
#'   \item{GIS}{The Gauge Instability Score for the gene.}
#'   \item{classification}{Character: \code{"gauge-sensitive"} or
#'     \code{"gauge-stable"}.}
#'   \item{threshold_used}{The threshold value used for classification.}
#'
#' @examples
#' gis <- setNames(runif(100, 0, 0.5), paste0("Gene", 1:100))
#' classification <- classify_genes(gis, threshold = "auto")
#' table(classification$classification)
#'
#' @export
classify_genes <- function(gis_vector, threshold = "auto") {
  if (!is.numeric(gis_vector)) {
    stop("gis_vector must be a numeric vector")
  }
  if (length(gis_vector) == 0) {
    stop("gis_vector must have length > 0")
  }

  genes <- names(gis_vector)
  if (is.null(genes)) {
    genes <- paste0("Gene", seq_along(gis_vector))
    names(gis_vector) <- genes
  }

  # Determine threshold
  if (is.character(threshold)) {
    threshold <- match.arg(threshold, c("auto", "median"))
    if (threshold == "auto") {
      thresh_val <- stats::quantile(gis_vector, probs = 0.75, na.rm = TRUE)
    } else if (threshold == "median") {
      thresh_val <- stats::median(gis_vector, na.rm = TRUE)
    }
  } else if (is.numeric(threshold) && length(threshold) == 1) {
    thresh_val <- threshold
  } else {
    stop("threshold must be 'auto', 'median', or a single numeric value")
  }

  classification <- ifelse(gis_vector > thresh_val,
                           "gauge-sensitive", "gauge-stable")

  result <- data.frame(
    gene = genes,
    GIS = gis_vector,
    classification = classification,
    threshold_used = thresh_val,
    stringsAsFactors = FALSE
  )

  class(result) <- c("gaugerna_classification", "data.frame")
  return(result)
}

#' Print method for gaugerna_classification objects
#' @param x A gaugerna_classification object
#' @param ... Additional arguments (ignored)
#' @export
print.gaugerna_classification <- function(x, ...) {
  cat("GaugeRNA Gene Classification\n")
  cat("============================\n")
  n_sensitive <- sum(x$classification == "gauge-sensitive")
  n_stable <- sum(x$classification == "gauge-stable")
  cat(sprintf("  Gauge-sensitive: %d genes\n", n_sensitive))
  cat(sprintf("  Gauge-stable: %d genes\n", n_stable))
  cat(sprintf("  Threshold: GIS > %.2f\n\n", x$threshold_used[1]))
  if (n_sensitive > 0) {
    cat("  Top gauge-sensitive genes:\n")
    sensitive <- x[x$classification == "gauge-sensitive", ]
    sensitive <- sensitive[order(-sensitive$GIS), ]
    n_show <- min(10, nrow(sensitive))
    print.data.frame(sensitive[1:n_show, c("gene", "GIS", "classification")],
          row.names = FALSE, digits = 4)
  }
  invisible(x)
}