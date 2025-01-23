# Use a comprehensive R image that includes RStudio
FROM rocker/rstudio:4.3.0

# Install additional system dependencies
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    zlib1g-dev \       # Required for R package dependencies
    wget \
    git \
    && rm -rf /var/lib/apt/lists/*

# Install R libraries using BiocManager
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('sangerseqR', 'Biostrings'))"

# Set working directory
WORKDIR /workspace
