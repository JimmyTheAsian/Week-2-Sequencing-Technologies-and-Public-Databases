# Use a comprehensive R development image
FROM rocker/rstudio:4.3.0

# Install additional system dependencies, including zlib development library
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    wget \
    git \
    zlib1g-dev \  # Add this to install zlib development library
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN Rscript -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager'); BiocManager::install(c('sangerseqR', 'Biostrings'))"

# Set working directory
WORKDIR /workspace
