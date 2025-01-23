# Load required libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("sangerseqR", "Biostrings"))

library(sangerseqR)
library(Biostrings)

# Define constants
ab1_file <- "your_file.ab1"         # Replace with your actual AB1 file path
vector_fasta <- "pGEM-T.fasta"      # Replace with the FASTA file of the pGEM-T vector

# Step 1: Read the AB1 File
cat("Reading AB1 file...\n")
sanger_data <- read.abif(ab1_file)
sanger_seq <- sangerseq(sanger_data)

# Extract the primary sequence (raw base calls from the chromatogram)
raw_sequence <- primarySeq(sanger_seq)

# Save the raw sequence as a FASTA file for BLAST
cat("Saving raw sequence to FASTA...\n")
writeXStringSet(DNAStringSet(raw_sequence), "raw_sequence.fasta")

# Step 2: BLAST Against the Vector
cat("Performing BLAST alignment against vector...\n")

# Create a BLAST database for the vector
system(paste("makeblastdb -in", vector_fasta, "-dbtype nucl -out vector_db"), intern = TRUE)

# Run BLAST using the raw sequence
system(paste(
  "blastn -query raw_sequence.fasta -db vector_db -out blast_results.txt -outfmt '6 qstart qend sstart send'",
  sep = " "
), intern = TRUE)

# Step 3: Parse BLAST Results
cat("Parsing BLAST results...\n")
if (file.exists("blast_results.txt") && file.info("blast_results.txt")$size > 0) {
  blast_hits <- read.table("blast_results.txt", header = FALSE)
  colnames(blast_hits) <- c("qstart", "qend", "sstart", "send")

  # Ensure BLAST results are valid
  if (nrow(blast_hits) > 0) {
    vector_start <- min(blast_hits$qstart)
    vector_end <- max(blast_hits$qend)

    # Step 4: Validate Trimming Coordinates
    seq_length <- nchar(raw_sequence)

    # Ensure 'vector_start' and 'vector_end' are within valid bounds
    if (vector_start < 1 || vector_start > seq_length) {
      cat("Warning: 'vector_start' is out of bounds. Setting to 1.\n")
      vector_start <- 1
    }
    if (vector_end < 1 || vector_end > seq_length) {
      cat("Warning: 'vector_end' is out of bounds. Setting to sequence length.\n")
      vector_end <- seq_length
    }

    # Ensure 'vector_end' does not overlap 'vector_start'
    if (vector_end < vector_start) {
      stop("Error: BLAST trimming coordinates are invalid ('vector_end' < 'vector_start').")
    }

    # Step 5: Trim Vector Portions
    cat("Trimming vector portions...\n")
    trimmed_sequence <- subseq(raw_sequence, start = vector_end + 1, end = vector_start - 1)
    cat("Trimmed sequence:\n", as.character(trimmed_sequence), "\n")

    # Save the trimmed sequence as FASTA
    writeXStringSet(DNAStringSet(trimmed_sequence), "final_trimmed_sequence.fasta")
    cat("Final trimmed sequence saved as 'final_trimmed_sequence.fasta'\n")
  } else {
    stop("Error: No valid BLAST hits found.")
  }
} else {
  cat("No BLAST hits found. Saving raw sequence as 'final_trimmed_sequence.fasta'...\n")
  writeXStringSet(DNAStringSet(raw_sequence), "final_trimmed_sequence.fasta")
}

# Step 6: Plot Chromatogram
cat("Plotting chromatogram...\n")
chromatogram(sanger_seq)

# End of Script
cat("Processing complete.\n")
