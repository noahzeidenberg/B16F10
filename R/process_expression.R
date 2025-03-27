#!/usr/bin/env Rscript

# Load required libraries
library(affy)
library(mogene10sttranscriptcluster.db)  # Mouse Gene 1.0 ST annotation package
library(oligo)                           # For Affymetrix arrays

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: Rscript process_expression.R <input_dir> <output_dir>")
}

input_dir <- args[1]
output_dir <- args[2]

# Process expression data
process_expression <- function(input_dir, output_dir) {
  tryCatch({
    # Load the CEL files
    message("Reading CEL files...")
    cel_files <- list.files(input_dir, pattern = "\\.cel$|\\.CEL$|\\.cel.gz$|\\.CEL.gz$", full.names = TRUE)
    if(length(cel_files) == 0) {
      stop("No CEL files found in input directory")
    }
    
    # Read the data
    raw_data <- read.celfiles(cel_files)
    
    # Normalize data (RMA normalization)
    message("Normalizing data...")
    eset <- rma(raw_data)
    
    # Get expression values
    expression_matrix <- exprs(eset)
    
    # Save processed data
    message("Saving results...")
    saveRDS(expression_matrix, file = file.path(output_dir, "expression_matrix.rds"))
    write.csv(expression_matrix, file = file.path(output_dir, "expression_matrix.csv"))
    
    # Basic QC plots
    message("Generating QC plots...")
    pdf(file.path(output_dir, "qc_plots.pdf"))
    boxplot(expression_matrix, main="Expression Distribution", las=2)
    hist(expression_matrix, main="Expression Histogram")
    dev.off()
    
    message("Successfully processed expression data")
  }, error = function(e) {
    stop("Failed to process expression data: ", e$message)
  })
}

# Create output directory if it doesn't exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Execute processing
process_expression(input_dir, output_dir) 