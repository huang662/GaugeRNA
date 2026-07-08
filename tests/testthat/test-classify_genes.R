# Tests for classify_genes

test_that("classify_genes returns correct structure", {
  gis <- setNames(runif(100, 0, 0.5), paste0("Gene", 1:100))
  result <- classify_genes(gis, threshold = "auto")

  expect_s3_class(result, "data.frame")
  expect_true("gene" %in% colnames(result))
  expect_true("GIS" %in% colnames(result))
  expect_true("classification" %in% colnames(result))
  expect_true("threshold_used" %in% colnames(result))
})

test_that("classify_genes auto threshold uses 75th percentile", {
  gis <- setNames(runif(100, 0, 0.5), paste0("Gene", 1:100))
  result <- classify_genes(gis, threshold = "auto")

  expected_thresh <- as.numeric(quantile(gis, probs = 0.75))
  expect_equal(as.numeric(result$threshold_used[1]), expected_thresh)
})

test_that("classify_genes median threshold works", {
  gis <- setNames(runif(100, 0, 0.5), paste0("Gene", 1:100))
  result <- classify_genes(gis, threshold = "median")

  expected_thresh <- as.numeric(median(gis))
  expect_equal(as.numeric(result$threshold_used[1]), expected_thresh)
})

test_that("classify_genes numeric threshold works", {
  gis <- setNames(runif(100, 0, 0.5), paste0("Gene", 1:100))
  result <- classify_genes(gis, threshold = 0.3)

  expect_equal(result$threshold_used[1], 0.3)
})

test_that("classify_genes produces correct classification counts", {
  gis <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  names(gis) <- paste0("Gene", 1:5)
  result <- classify_genes(gis, threshold = 0.25)

  n_sensitive <- sum(result$classification == "gauge-sensitive")
  n_stable <- sum(result$classification == "gauge-stable")

  expect_equal(n_sensitive, 3)  # 0.3, 0.4, 0.5
  expect_equal(n_stable, 2)     # 0.1, 0.2
})

test_that("classify_genes generates gene names if none provided", {
  gis <- runif(10, 0, 0.5)
  result <- classify_genes(gis, threshold = "auto")
  expect_true(all(grepl("^Gene", result$gene)))
})

test_that("classify_genes errors on empty input", {
  expect_error(classify_genes(numeric(0)), "gis_vector must have length")
})