# Use a base image with R pre-installed
FROM rocker/r-ver:4.3.0

# Install system dependencies for R and BLAST
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \
    libxml2-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install required R packages
RUN Rscript -e "install.packages('BiocManager'); BiocManager::install(c('sangerseqR', 'Biostrings'))"

# Set working directory
WORKDIR /workspace
