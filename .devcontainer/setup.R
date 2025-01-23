if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("sangerseqR", "Biostrings"))

# Install other packages as needed
install.packages(c("tidyverse"))
