#!/usr/bin/env Rscript

# Load required libraries
library(minfi)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylation450kmanifest)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript process_methylation.R <input_dir> <output_dir>")
}

input_dir <- args[1]
output_dir <- args[2]

# Process methylation data
process_methylation <- function(input_dir, output_dir) {
  tryCatch({
    # Load the data
    rgSet <- read.metharray.exp(base = input_dir)
    
    # Preprocess and normalize
    mSet <- preprocessFunnorm(rgSet)
    
    # Get beta values
    beta_values <- getBeta(mSet)
    
    # Quality control metrics
    qc <- getQC(mSet)
    
    # Save processed data
    saveRDS(beta_values, file = file.path(output_dir, "beta_values.rds"))
    saveRDS(qc, file = file.path(output_dir, "qc_metrics.rds"))
    
    # Generate QC plots
    pdf(file.path(output_dir, "qc_plots.pdf"))
    plotQC(qc)
    dev.off()
    
    message("Successfully processed methylation data")
  }, error = function(e) {
    stop("Failed to process methylation data: ", e$message)
  })
}

# Execute processing
process_methylation(input_dir, output_dir) 