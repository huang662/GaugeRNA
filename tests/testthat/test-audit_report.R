# Tests for audit_report

test_that("audit_report returns correct structure", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 8
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(rep(0, 4), rep(1, 4))

  report <- audit_report(expr_mat, condition, n_perturbations = 5)

  expect_s3_class(report, "gaugerna_report")
  expect_true("summary" %in% names(report))
  expect_true("gis_results" %in% names(report))
  expect_true("classification" %in% names(report))
  expect_true("decomposition" %in% names(report))
  expect_true("recommendations" %in% names(report))
  expect_true("parameters" %in% names(report))
  expect_true("timestamp" %in% names(report))
})

test_that("audit_report summary has expected metrics", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 8
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(rep(0, 4), rep(1, 4))

  report <- audit_report(expr_mat, condition, n_perturbations = 5)

  expected_metrics <- c("n_genes", "n_samples", "best_method", "best_GIS",
                        "worst_method", "worst_GIS", "n_gauge_sensitive",
                        "n_gauge_stable", "gene_gauge_mean", "sample_gauge_sd",
                        "curvature_sd")
  expect_true(all(expected_metrics %in% report$summary$metric))
})

test_that("audit_report recommendations are non-empty", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 8
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(rep(0, 4), rep(1, 4))

  report <- audit_report(expr_mat, condition, n_perturbations = 5)

  expect_true(length(report$recommendations) > 0)
  expect_type(report$recommendations, "character")
})

test_that("audit_report validates input", {
  expr_mat <- matrix(rnbinom(100 * 6, mu = 100, size = 10), nrow = 100, ncol = 6)
  condition <- c(0, 0, 0, 1, 1, 1)

  expect_error(audit_report(expr_mat, condition[1:3]),
               "length\\(condition\\)")
  expect_error(audit_report(expr_mat, c(0, 0, 0, 1, 1, 2)),
               "condition must be a binary")
})

test_that("audit_report print method works", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 8
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(rep(0, 4), rep(1, 4))

  report <- audit_report(expr_mat, condition, n_perturbations = 5)

  expect_output(print(report), "GaugeRNA")
})

test_that("audit_report summary method works", {
  set.seed(42)
  n_genes <- 200
  n_samples <- 8
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(rep(0, 4), rep(1, 4))

  report <- audit_report(expr_mat, condition, n_perturbations = 5)

  s <- suppressMessages(summary(report))
  expect_s3_class(s, "data.frame")
})