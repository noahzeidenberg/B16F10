# RNA-seq Batch Correction Pipeline
# This script performs batch correction across all datasets

# Load required packages
cat("=== Starting RNA-seq Batch Correction Pipeline ===\n")
cat("Loading required packages...\n")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "sva",
  "edgeR",
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  cat(sprintf("Loading package: %s\n", pkg))
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Function to collect all feature counts files
collect_feature_counts <- function(base_dir) {
  cat("Collecting all feature counts files...\n")
  
  # Find all GSE directories
  gse_dirs <- list.dirs(base_dir, recursive = FALSE)
  gse_dirs <- gse_dirs[grep("^GSE", basename(gse_dirs))]
  
  cat(sprintf("Found %d GSE directories\n", length(gse_dirs)))
  
  # Collect all feature counts files
  all_counts <- list()
  all_samples <- c()
  all_gse_ids <- c()
  
  for (gse_dir in gse_dirs) {
    gse_id <- basename(gse_dir)
    counts_file <- file.path(gse_dir, "results", "counts", "feature_counts.rds")
    
    if (file.exists(counts_file)) {
      cat(sprintf("Loading counts from %s\n", gse_id))
      counts_obj <- readRDS(counts_file)
      
      # Extract counts and sample names
      counts <- counts_obj$counts
      samples <- colnames(counts)
      
      # Add GSE ID to sample names to make them unique
      new_samples <- paste0(gse_id, "_", samples)
      colnames(counts) <- new_samples
      
      # Add to lists
      all_counts[[gse_id]] <- counts
      all_samples <- c(all_samples, new_samples)
      all_gse_ids <- c(all_gse_ids, rep(gse_id, length(samples)))
    } else {
      cat(sprintf("Warning: No feature counts file found for %s\n", gse_id))
    }
  }
  
  if (length(all_counts) == 0) {
    cat("No feature counts files found, cannot proceed\n")
    return(NULL)
  }
  
  # Combine all counts
  cat("Combining all counts...\n")
  combined_counts <- do.call(cbind, all_counts)
  
  # Create batch information
  batch_info <- data.frame(
    sample = all_samples,
    gse_id = all_gse_ids,
    batch = factor(all_gse_ids)  # Each GSE is a batch
  )
  
  return(list(
    counts = combined_counts,
    batch_info = batch_info
  ))
}

# Function to perform batch correction
perform_batch_correction <- function(counts, batch_info, output_dir) {
  cat("Starting batch correction...\n")
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Remove genes with zero counts across all samples
  keep <- rowSums(counts) > 0
  counts <- counts[keep, ]
  cat(sprintf("Removed %d genes with zero counts\n", sum(!keep)))
  
  # Check if we have multiple batches
  num_batches <- length(unique(batch_info$batch))
  cat(sprintf("Number of batches: %d\n", num_batches))
  
  if (num_batches > 1) {
    # ComBat-seq batch correction
    cat("Performing ComBat-seq batch correction...\n")
    corrected_counts <- ComBat_seq(counts = as.matrix(counts),
                                  batch = batch_info$batch)
  } else {
    # Skip batch correction if we only have one batch
    cat("Skipping batch correction (only one batch detected)...\n")
    corrected_counts <- as.matrix(counts)
  }
  
  # Save results
  cat("Saving batch-corrected counts...\n")
  results <- list(
    raw_counts = counts,
    corrected_counts = corrected_counts,
    batch_info = batch_info
  )
  
  output_file <- file.path(output_dir, "batch_corrected_counts.rds")
  saveRDS(results, output_file)
  cat(sprintf("Saved batch-corrected counts to: %s\n", output_file))
  
  return(results)
}

# Main workflow
main <- function() {
  # Set up directories
  base_dir <- getwd()
  output_dir <- file.path(base_dir, "results", "batch_correction")
  
  cat("Setting up directories...\n")
  cat(sprintf("Base directory: %s\n", base_dir))
  cat(sprintf("Output directory: %s\n", output_dir))
  
  # Check if batch correction has already been performed
  if (file.exists(file.path(output_dir, "batch_corrected_counts.rds"))) {
    cat("Batch correction has already been performed, skipping\n")
    return(0)
  }
  
  # Collect all feature counts
  all_data <- collect_feature_counts(base_dir)
  
  if (is.null(all_data)) {
    cat("Failed to collect feature counts, cannot proceed\n")
    return(1)
  }
  
  # Perform batch correction
  results <- perform_batch_correction(
    all_data$counts,
    all_data$batch_info,
    output_dir
  )
  
  cat("=== RNA-seq Batch Correction Pipeline Completed Successfully ===\n")
  cat("Batch correction has been performed. Please run normalize_counts.R next.\n")
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 