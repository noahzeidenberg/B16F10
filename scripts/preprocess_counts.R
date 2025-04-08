# Script to preprocess normalized counts to ensure they're all non-negative
# This script will load the normalized counts, replace negative values with zeros,
# and save the preprocessed counts back to the same file

# Load required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c(
  "stringr",
  "data.table"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
  library(pkg, character.only = TRUE)
}

# Set up directories
base_dir <- path.expand("~/scratch/B16F10")
gse_id <- "GSE287957"
gse_dir <- file.path(base_dir, gse_id)
normalization_dir <- file.path(gse_dir, "results", "normalization")

# Load normalized counts
normalized_file <- file.path(normalization_dir, "normalized_counts.rds")
cat(sprintf("Loading normalized counts from: %s\n", normalized_file))

if (!file.exists(normalized_file)) {
  cat(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
  quit(status = 1)
}

normalized_data <- readRDS(normalized_file)

# Extract counts
if (!"normalized_counts" %in% names(normalized_data)) {
  cat("Error: 'normalized_counts' not found in normalized data. Available keys:\n")
  cat(paste(names(normalized_data), collapse = ", "), "\n")
  quit(status = 1)
}

counts <- normalized_data$normalized_counts

# Check for negative values
negative_count <- sum(counts < 0)
cat(sprintf("Found %d negative values in normalized counts\n", negative_count))

if (negative_count > 0) {
  # Replace negative values with zeros
  cat("Replacing negative values with zeros...\n")
  counts[counts < 0] <- 0
  
  # Update the normalized data
  normalized_data$normalized_counts <- counts
  
  # Save the preprocessed counts
  cat(sprintf("Saving preprocessed counts to: %s\n", normalized_file))
  saveRDS(normalized_data, normalized_file)
  
  cat("Preprocessing completed successfully.\n")
} else {
  cat("No negative values found. No preprocessing needed.\n")
} 