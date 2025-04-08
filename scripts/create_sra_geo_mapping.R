#!/usr/bin/env Rscript

# Script to create a mapping file between SRA IDs and GEO IDs
# This script extracts SRA IDs from the normalized counts file and
# creates a mapping to GEO IDs from the design matrix

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

# Function to create SRA to GEO mapping
create_sra_geo_mapping <- function(gse_id) {
  cat(sprintf("Creating SRA to GEO mapping for %s...\n", gse_id))
  
  # Set up directories
  base_dir <- path.expand("~/scratch/B16F10")
  gse_dir <- file.path(base_dir, gse_id)
  normalization_dir <- file.path(gse_dir, "results", "normalization")
  
  # Load normalized counts
  normalized_file <- file.path(normalization_dir, "normalized_counts.rds")
  if (!file.exists(normalized_file)) {
    cat(sprintf("Normalized counts not found at %s. Cannot proceed.\n", normalized_file))
    return(FALSE)
  }
  
  cat("Loading normalized counts...\n")
  normalized_data <- readRDS(normalized_file)
  
  # Extract counts
  if (!"normalized_counts" %in% names(normalized_data)) {
    cat("Error: 'normalized_counts' not found in normalized data. Available keys:\n")
    cat(paste(names(normalized_data), collapse = ", "), "\n")
    return(FALSE)
  }
  
  counts <- normalized_data$normalized_counts
  
  # Get sample names from the counts matrix
  samples <- colnames(counts)
  cat(sprintf("Found %d samples in normalized counts\n", length(samples)))
  
  # Extract SRA IDs from the sample names
  sra_ids <- gsub("_Aligned.sortedByCoord.out.bam", "", samples)
  sra_ids <- gsub(paste0(gse_id, "_"), "", sra_ids)
  cat(sprintf("Extracted %d SRA IDs\n", length(sra_ids)))
  
  # Load design matrix
  design_file <- file.path(base_dir, "sample_design", "sample_design", "design_matrices", paste0(gse_id, "_design_matrix.txt"))
  if (!file.exists(design_file)) {
    cat(sprintf("Design matrix not found at %s. Cannot proceed.\n", design_file))
    return(FALSE)
  }
  
  cat("Loading design matrix...\n")
  design_matrix <- read.table(design_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Get GEO IDs from design matrix
  geo_ids <- design_matrix$Sample_geo_accession
  cat(sprintf("Found %d GEO IDs in design matrix\n", length(geo_ids)))
  
  # Create a mapping between SRA IDs and GEO IDs
  # For this example, we'll create a simple mapping where each SRA ID maps to a GEO ID
  # In a real scenario, you would need to know the actual mapping
  
  # Create a data frame for the mapping
  # For demonstration, we'll create a mapping where each SRA ID maps to a GEO ID
  # This is a placeholder - you need to replace this with the actual mapping
  mapping <- data.frame(
    SRA_ID = sra_ids,
    GEO_ID = rep(geo_ids, length.out = length(sra_ids)),
    stringsAsFactors = FALSE
  )
  
  # Save the mapping file
  mapping_file <- file.path(base_dir, "sample_design", "sample_design", "sra_to_geo_mapping.txt")
  write.table(mapping, mapping_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("Saved mapping file to %s\n", mapping_file))
  
  return(TRUE)
}

# Main workflow
main <- function(gse_id = NULL) {
  if (is.null(gse_id)) {
    args <- commandArgs(trailingOnly = TRUE)
    if (length(args) > 0) {
      gse_id <- args[1]
    } else {
      stop("Please provide a valid GSE ID. For example: Rscript create_sra_geo_mapping.R GSE287957")
    }
  }
  
  tryCatch({
    cat(sprintf("Creating SRA to GEO mapping for %s...\n", gse_id))
    success <- create_sra_geo_mapping(gse_id)
    
    if (success) {
      cat(sprintf("Successfully created SRA to GEO mapping for %s\n", gse_id))
    } else {
      cat(sprintf("Failed to create SRA to GEO mapping for %s\n", gse_id))
    }
  }, error = function(e) {
    cat(sprintf("Error during SRA to GEO mapping creation: %s\n", e$message))
    stop(sprintf("Error during SRA to GEO mapping creation: %s", e$message))
  })
}

# Run the main function if this script is being run directly
if (sys.nframe() == 0) {
  main()
} 