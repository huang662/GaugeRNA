# Dockerfile for GaugeRNA
# Reproducible execution environment for the GaugeRNA R package

FROM rocker/r-ver:4.3.3

LABEL maintainer="GaugeRNA Research Group"
LABEL description="GaugeRNA: A Gauge-Invariant Audit Framework for RNA-seq Differential Expression"
LABEL version="0.1.0"

# Install system dependencies
RUN apt-get update && apt-get install -y \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    && rm -rf /var/lib/apt/lists/*

# Install R package dependencies
RUN install2.r --error --skipinstalled \
    testthat \
    knitr \
    rmarkdown \
    ggplot2 \
    DESeq2

# Copy package source
COPY . /GaugeRNA

# Install GaugeRNA
RUN R -e 'install.packages("/GaugeRNA", repos = NULL, type = "source")'

# Set working directory
WORKDIR /workspace

# Default command
CMD ["R"]