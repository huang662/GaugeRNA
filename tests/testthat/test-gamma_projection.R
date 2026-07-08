# Tests for gamma_projection

test_that("gamma_projection returns correct length", {
  logfc <- rnorm(100, mean = 0.5, sd = 1.0)
  gamma <- gamma_projection(logfc)
  expect_equal(length(gamma), length(logfc))
})

test_that("gamma_projection has median approximately zero", {
  logfc <- rnorm(100, mean = 0.5, sd = 1.0)
  gamma <- gamma_projection(logfc)
  expect_equal(median(gamma), 0, tolerance = 1e-10)
})

test_that("gamma_projection handles NA values", {
  logfc <- c(rnorm(50), NA, rnorm(49))
  gamma <- gamma_projection(logfc)
  expect_equal(length(gamma), length(logfc))
  expect_false(any(is.na(gamma)))
})

test_that("gamma_projection maintains relative ordering", {
  logfc <- rnorm(100)
  gamma <- gamma_projection(logfc)
  # Differences should be preserved
  expect_equal(diff(logfc), diff(gamma))
})

test_that("gamma_projection errors on non-numeric input", {
  expect_error(gamma_projection(c("a", "b", "c")), "logfc must be a numeric vector")
})

test_that("gamma_projection errors on empty input", {
  expect_error(gamma_projection(numeric(0)), "logfc must have length > 0")
})