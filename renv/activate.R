# renv activation script
# Run this to activate renv for the GaugeRNA project
# source("renv/activate.R")

local({
  # Check if renv is available
  if (requireNamespace("renv", quietly = TRUE)) {
    renv::activate()
    cat("renv activated for GaugeRNA project.\n")
  } else {
    cat("renv is not installed. To install: install.packages('renv')\n")
    cat("Then run: renv::activate()\n")
  }
})