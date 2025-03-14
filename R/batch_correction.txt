#!/usr/bin/env Rscript

# Load required libraries
library(sva)
library(limma)
library(DESeq2)

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
    stop("Usage: Rscript batch_correction.R data_files.txt batch_info.csv data_type output_dir")
}

data_files_path <- args[1]
batch_info_path <- args[2]
data_type <- args[3]
output_dir <- args[4]

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Read file paths
data_files <- readLines(data_files_path)
message("Processing files:")
for (file in data_files) {
    message(file)
}

# Read batch information
batch_info <- read.csv(batch_info_path)
batch <- batch_info$batch

# Function to read and combine data
combine_data <- function(files, data_type) {
    if (data_type == "rnaseq") {
        # Read and combine RNA-seq count data
        data_list <- lapply(files, function(f) {
            message("Reading ", f)
            read.csv(f, row.names = 1)
        })
        combined_data <- do.call(cbind, data_list)
    } else {
        # Read and combine methylation beta values
        data_list <- lapply(files, function(f) {
            message("Reading ", f)
            readRDS(f)
        })
        combined_data <- do.call(cbind, data_list)
    }
    return(combined_data)
}

# Combine data
message("Combining data matrices...")
combined_data <- combine_data(data_files, data_type)

# Perform batch correction
message("Performing batch correction...")
if (data_type == "rnaseq") {
    # RNA-seq specific batch correction
    # Convert to log2 counts
    log_counts <- log2(combined_data + 1)
    
    # Combat-seq batch correction
    corrected_data <- ComBat_seq(
        counts = combined_data,
        batch = batch,
        group = NULL
    )
    
    # Save results
    write.csv(
        corrected_data,
        file.path(output_dir, "batch_corrected_counts.csv")
    )
    
} else {
    # Methylation array batch correction
    # ComBat batch correction for beta values
    corrected_data <- ComBat(
        dat = combined_data,
        batch = batch,
        mod = NULL,
        par.prior = TRUE,
        prior.plots = FALSE
    )
    
    # Save results
    write.csv(
        corrected_data,
        file.path(output_dir, "batch_corrected_beta_values.csv")
    )
}

message("Batch correction completed successfully") 