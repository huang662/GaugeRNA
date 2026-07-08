#' Compute Gamma Projection
#'
#' Computes the gamma projection of log-fold changes, which removes the
#' gene-independent gauge component from raw log-fold changes.
#'
#' The gamma projection is defined as:
#' \deqn{\Gamma_g = \Delta_g^{raw} - \text{median}(\Delta^{raw})}
#' where \eqn{\Delta_g^{raw}} is the raw log-fold change for gene g.
#'
#' This centering operation removes the global shift that is common to all
#' genes and represents the arbitrary gauge choice of normalization.
#'
#' @param logfc A numeric vector of raw log-fold changes (one per gene).
#'
#' @return A numeric vector of gamma-projected log-fold changes.
#'
#' @examples
#' logfc <- rnorm(100, mean = 0.5, sd = 1.0)
#' gamma <- gamma_projection(logfc)
#' # The median of gamma should be approximately 0
#' median(gamma)
#'
#' @export
gamma_projection <- function(logfc) {
  if (!is.numeric(logfc)) {
    stop("logfc must be a numeric vector")
  }
  if (length(logfc) == 0) {
    stop("logfc must have length > 0")
  }

  med <- stats::median(logfc, na.rm = TRUE)
  gamma <- logfc - med

  # Set NA/infinite values to 0
  gamma[!is.finite(gamma)] <- 0

  return(gamma)
}