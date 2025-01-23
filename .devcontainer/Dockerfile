# Use a comprehensive R development image
FROM rocker/rstudio:4.3.0

# Install additional system dependencies
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    zlib1g-dev \
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Pre-install R packages using BiocManager to speed up container initialization
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('sangerseqR', 'Biostrings'))"

# Set the working directory
WORKDIR /workspace
