# Tests for gauge_decomposition

test_that("gauge_decomposition returns correct structure", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)

  result <- gauge_decomposition(expr_mat)

  expect_type(result, "list")
  expect_true("mu" %in% names(result))
  expect_true("a" %in% names(result))
  expect_true("epsilon" %in% names(result))
  expect_true("log_counts" %in% names(result))
})

test_that("gauge_decomposition dimensions are correct", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)

  result <- gauge_decomposition(expr_mat)

  expect_equal(length(result$mu), n_genes)
  expect_equal(length(result$a), n_samples)
  expect_equal(dim(result$epsilon), c(n_genes, n_samples))
})

test_that("gauge_decomposition curvature has mean near zero", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)

  result <- gauge_decomposition(expr_mat)

  expect_equal(mean(result$epsilon), 0, tolerance = 1e-10)
})

test_that("gauge_decomposition reconstruction is correct", {
  set.seed(42)
  n_genes <- 50
  n_samples <- 4
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)

  result <- gauge_decomposition(expr_mat)

  mu_mat <- matrix(result$mu, nrow = n_genes, ncol = n_samples)
  a_mat <- matrix(result$a, nrow = n_genes, ncol = n_samples, byrow = TRUE)
  reconstructed <- mu_mat + a_mat + result$epsilon

  expect_equal(as.numeric(reconstructed), as.numeric(result$log_counts), tolerance = 1e-10)
})

test_that("gauge_decomposition sample gauge is centered", {
  set.seed(42)
  n_genes <- 100
  n_samples <- 6
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)

  result <- gauge_decomposition(expr_mat)

  expect_equal(mean(result$a), 0, tolerance = 1e-10)
})

test_that("gauge_decomposition handles named matrix", {
  set.seed(42)
  n_genes <- 50
  n_samples <- 4
  expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10),
                     nrow = n_genes, ncol = n_samples)
  rownames(expr_mat) <- paste0("Gene", 1:n_genes)
  colnames(expr_mat) <- paste0("Sample", 1:n_samples)

  result <- gauge_decomposition(expr_mat)

  expect_equal(names(result$mu), rownames(expr_mat))
  expect_equal(names(result$a), colnames(expr_mat))
})