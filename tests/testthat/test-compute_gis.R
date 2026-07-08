# Tests for compute_gis

test_that("compute_gis returns correct structure", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(0, 0, 0, 1, 1, 1)

  result <- compute_gis(expr_mat, condition, n_perturbations = 5)

  expect_s3_class(result, "data.frame")
  expect_true("method" %in% colnames(result))
  expect_true("GIS" %in% colnames(result))
  expect_true("mean_stability" %in% colnames(result))
  expect_true("sd_stability" %in% colnames(result))
})

test_that("compute_gis GIS values are between 0 and 1", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(0, 0, 0, 1, 1, 1)

  result <- compute_gis(expr_mat, condition, n_perturbations = 5,
                        methods = c("raw", "gamma"))

  expect_true(all(result$GIS >= 0 & result$GIS <= 1))
})

test_that("compute_gis errors on invalid input", {
  expr_mat <- matrix(rnbinom(100 * 6, mu = 100, size = 10), nrow = 100, ncol = 6)
  condition <- c(0, 0, 0, 1, 1, 1)

  expect_error(compute_gis(expr_mat, condition[1:3]),
               "length\\(condition\\) must equal")
  expect_error(compute_gis(expr_mat, c(0, 0, 0, 1, 1, 2)),
               "condition must be a binary")
  expect_error(compute_gis(expr_mat, condition, n_perturbations = 0),
               "n_perturbations must be at least 1")
})

test_that("compute_gis handles different method specifications", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  condition <- c(0, 0, 0, 1, 1, 1)

  # Single method
  result1 <- compute_gis(expr_mat, condition, n_perturbations = 3,
                         methods = "raw")
  expect_equal(nrow(result1), 1)

  # All methods
  result2 <- compute_gis(expr_mat, condition, n_perturbations = 3,
                         methods = c("raw", "gamma", "deseq2"))
  expect_equal(nrow(result2), 3)
})