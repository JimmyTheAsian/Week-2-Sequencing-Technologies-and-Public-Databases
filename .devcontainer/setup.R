if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required R packages
BiocManager::install(c("sangerseqR", "Biostrings"))

# Add any other R packages you need
install.packages(c("tidyverse", "ggplot2"))
