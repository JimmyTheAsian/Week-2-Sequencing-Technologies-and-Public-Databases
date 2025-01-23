if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install required R packages
BiocManager::install(c("sangerseqR", "Biostrings"))

# Add additional R packages as needed
install.packages(c("tidyverse", "ggplot2"))
