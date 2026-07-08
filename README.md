# GaugeRNA: A Gauge-Invariant Audit Framework for RNA-seq Differential Expression

[![R >= 3.5.0](https://img.shields.io/badge/R-%3E%3D%203.5.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

GaugeRNA is an R package for auditing RNA-seq differential expression analysis through the lens of gauge invariance. The framework quantifies the sensitivity of differential expression results to arbitrary normalization choices (gauge choices) via the **Gauge Instability Score (GIS)**.

### Key Features

- **Gauge Instability Score (GIS)**: Quantify how much differential expression results change under random gauge perturbations
- **Gamma Projection**: Remove the gene-independent gauge component from log-fold changes
- **Gauge Decomposition**: Factorize the log-count matrix into gene gauge, sample gauge, and curvature components
- **Gene Classification**: Identify gauge-sensitive vs. gauge-stable genes
- **Audit Reports**: Generate comprehensive, reproducible audit reports

## Installation

### From Source

```r
# Install dependencies
install.packages(c("testthat", "knitr", "rmarkdown"))

# Install GaugeRNA
install.packages("path/to/GaugeRNA", repos = NULL, type = "source")

# Or using devtools
# devtools::install("path/to/GaugeRNA")
```

### Using Docker

```bash
docker build -t gaugerna .
docker run -it gaugerna R
```

## Quick Start

```r
library(GaugeRNA)

# Simulate expression data
set.seed(42)
n_genes <- 500
n_samples <- 10
expr_mat <- matrix(rnbinom(n_genes * n_samples, mu = 100, size = 10), 
                   nrow = n_genes, ncol = n_samples)
condition <- c(rep(0, 5), rep(1, 5))

# Generate a complete audit report
report <- audit_report(expr_mat, condition, n_perturbations = 20)

# Print the report
print(report)

# View summary
summary(report)

# Examine gene classification
head(report$classification)
```

## Core Functions

| Function | Description |
|----------|-------------|
| `compute_gis()` | Compute Gauge Instability Scores for multiple methods |
| `gamma_projection()` | Apply gamma projection to log-fold changes |
| `classify_genes()` | Classify genes as gauge-sensitive or gauge-stable |
| `gauge_decomposition()` | Decompose log-count matrix into gauge components |
| `audit_report()` | Generate a complete audit report |

## Documentation

After installation, access the vignette:

```r
browseVignettes("GaugeRNA")
```

Or view help for individual functions:

```r
?compute_gis
?gamma_projection
?classify_genes
?gauge_decomposition
?audit_report
```

## Package Structure

```
GaugeRNA/
  DESCRIPTION           # Package metadata
  NAMESPACE             # Exports and imports
  R/                    # Core R functions
  man/                  # Documentation (.Rd files)
  vignettes/            # Package vignettes
  tests/                # testthat unit tests
  inst/extdata/         # Example datasets
  data/                 # Pre-computed paper results
  benchmark_scripts/    # Scripts to reproduce paper figures
  Dockerfile            # Docker container definition
  renv/                 # Dependency management
```

## Reproducibility

The package includes:
- `inst/extdata/simulated_counts_matrix.rda`: Example simulated RNA-seq dataset
- `data/gis_paper_results.rda`: Pre-computed GIS results from the paper
- `benchmark_scripts/`: R scripts to reproduce all figures in the manuscript
- `Dockerfile`: Containerized execution environment
- `renv.lock`: Locked dependency versions

## Citation

If you use GaugeRNA in your research, please cite:

```
GaugeRNA Research Group. GaugeRNA: A Gauge-Invariant Audit Framework 
for RNA-seq Differential Expression. R package version 0.1.0.
```

## License

MIT License. See LICENSE file for details.