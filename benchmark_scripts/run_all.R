###########################################################################
# Run All Benchmarks
# Executes all benchmark scripts to reproduce paper figures
###########################################################################

cat("========================================\n")
cat("  GaugeRNA Benchmark Suite\n")
cat("  Reproducing all paper figures\n")
cat("========================================\n\n")

benchmark_dir <- "benchmark_scripts"

scripts <- c(
  "figure1_gis_distribution.R",
  "figure2_gamma_projection.R",
  "figure3_gauge_decomposition.R",
  "figure4_gis_vs_expression.R"
)

start_time <- Sys.time()

for (script in scripts) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n")
  cat("Running:", script, "\n")
  cat(paste(rep("=", 60), collapse = ""), "\n\n")
  
  script_path <- file.path(benchmark_dir, script)
  if (file.exists(script_path)) {
    source(script_path)
  } else {
    cat("WARNING: Script not found:", script_path, "\n")
  }
}

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("All benchmarks completed.\n")
cat("Total elapsed time:", round(elapsed, 2), "minutes\n")
cat(paste(rep("=", 60), collapse = ""), "\n")