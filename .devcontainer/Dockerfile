# Use a development-focused R image with common tools
FROM rocker/dev:4.3.0

# Install additional dependencies
RUN apt-get update && apt-get install -y \
    ncbi-blast+ \   # For BLAST tools
    wget \           # For downloading files
    && rm -rf /var/lib/apt/lists/*

# Install R libraries using BiocManager
RUN Rscript -e "install.packages('BiocManager'); BiocManager::install(c('sangerseqR', 'Biostrings'))"

# Set working directory inside the container
WORKDIR /workspace
