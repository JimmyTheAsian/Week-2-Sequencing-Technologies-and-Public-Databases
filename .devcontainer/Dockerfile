# Use a comprehensive R development image
FROM rocker/rstudio:4.3.0

# Install system dependencies, including BLAST and general utilities
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    wget \
    curl \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages via Bioconductor
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('sangerseqR', 'Biostrings'))"

# Set a working directory
WORKDIR /workspace
